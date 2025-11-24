createCropLayFinalYear1 <- function(Lay, 
                                    potLayTopValid, 
                                    runClusteringInParallel, 
                                    clusterDistance, 
                                    studyAreaHash) {
  # 0) input validation
  if (!is.numeric(clusterDistance) || length(clusterDistance) != 1) {
    stop("`clusterDistance` must be a single numeric value")
  }
  if (clusterDistance < 0) {
    stop("`clusterDistance` must be non-negative")
  }
  
  # 1) crop/reproject Lay to potLayTopValid
  cropLay <- Cache(postProcessTo, Lay, potLayTopValid)
  
  # If nothing inside the mask, return consistent list (old callers now expect list)
  if (nrow(cropLay) == 0L) {
    warning("createCropLayFinalYear1: no input lines fall within potLayTopValid — returning empty lines + original availableArea")
    return(list(
      lines         = Lay[0, ],            # empty SpatVector(lines), preserves schema/CRS
      availableArea = potLayTopValid       # unchanged
    ))
  }
  
  # 2) CRITICAL — replicate old semantics: split lines at potential boundaries and
  #    attach polygon Potential. Using polygon as FIRST arg ensures its attributes
  #    keep their original names (no suffix) while any duplicates from lines get suffixes (.1)
  cropLay <- Cache(terra::intersect, potLayTopValid, cropLay)
  
  # 2a) keep exactly ONE Potential column — the polygon's
  potCols <- grep("^Potential", names(cropLay), value = TRUE)
  if (length(potCols) == 0) {
    stop("createCropLayFinalYear1: no 'Potential' column after intersect — check inputs")
  }
  # keep the bare "Potential" if present; otherwise keep the first and rename it
  keepPot <- if ("Potential" %in% potCols) "Potential" else potCols[1]
  if (!identical(keepPot, "Potential")) names(cropLay)[names(cropLay) == keepPot] <- "Potential"
  dropCols <- setdiff(potCols, "Potential")
  if (length(dropCols)) cropLay <- cropLay[, setdiff(names(cropLay), dropCols)]
  
  # Guard again in case everything got filtered out
  if (nrow(cropLay) == 0L) {
    warning("createCropLayFinalYear1: no line segments remain after intersect — returning empty lines + original availableArea")
    return(list(lines = Lay[0, ], availableArea = potLayTopValid))
  }
  
  # 3) exclude existing lines from available potential area (50 m buffer as before)
  cropLayBuf <- Cache(buffer,    cropLay, width = 50)
  cropLayAg  <- Cache(aggregate, cropLayBuf, dissolve = TRUE)
  finalPotLay <- Cache(erase, potLayTopValid, cropLayAg)
  
  # 4) unique IDs (used later by clustering/simulation)
  cropLay$individualID <- seq_len(nrow(cropLay))
  
  # 5) overlap removal (keep your improved NA-angle guard)
  tic("Time elapsed to identify overlapping features: ")
  overlapMatrix <- data.table::data.table()
  if (nrow(cropLay) > 1L) {
    overlapMatrix <- tryCatch({
      cropLay_sf <- sf::st_as_sf(cropLay)
      relList <- sf::st_relate(cropLay_sf, pattern = "T********", sparse = TRUE)
      if (!length(relList)) {
        data.table::data.table()
      } else {
        counts <- lengths(relList)
        if (!any(counts)) {
          data.table::data.table()
        } else {
          data.table::data.table(
            id.1 = rep.int(seq_along(relList), counts),
            id.2 = unlist(relList, use.names = FALSE)
          )
        }
      }
    }, error = function(e) {
      if (isTRUE(getOption("run_scenario.debug", FALSE))) {
        message("[createCropLayFinalYear1] st_relate failed: ", conditionMessage(e),
                "; falling back to terra::relate for overlap detection.")
      }
      fallback <- tryCatch(
        terra::relate(x = cropLay, relation = "T********", pairs = TRUE),
        error = function(e2) {
          warning("createCropLayFinalYear1: overlap detection failed (",
                  conditionMessage(e2), "); skipping overlap pruning.",
                  immediate. = TRUE)
          data.table::data.table()
        }
      )
      data.table::as.data.table(fallback)
    })
    overlapMatrix <- data.table::as.data.table(overlapMatrix)
    if (nrow(overlapMatrix)) {
      overlapMatrix <- overlapMatrix[id.1 != id.2]
      if (!nrow(overlapMatrix)) {
        overlapMatrix <- overlapMatrix[0, ]
      }
    }
    if (nrow(overlapMatrix)) {
      # normalize (id.1 < id.2) and unique
      overlapMatrix[, c("id.1", "id.2") := .(pmin(id.1, id.2), pmax(id.1, id.2))]
      overlapMatrix <- unique(overlapMatrix, by = c("id.1", "id.2"))
    }
  }
  to_remove <- integer()
  for (i in seq_len(nrow(overlapMatrix))) {
    if (i %in% c(1, 1000) || i %% 10000 == 0)
      message(sprintf("Feature %d of %d", i, nrow(overlapMatrix)))
    a1 <- calculateLineAngle(cropLay[overlapMatrix[i, "id.1"]])
    a2 <- calculateLineAngle(cropLay[overlapMatrix[i, "id.2"]])
    if (is.na(a1) || is.na(a2)) next
    if (round(a1, 4) != round(a2, 4)) next
    l1 <- perim(cropLay[overlapMatrix[i, "id.1"]])
    l2 <- perim(cropLay[overlapMatrix[i, "id.2"]])
    smallest <- if (l1 > l2) overlapMatrix[i, id.2] else overlapMatrix[i, id.1]
    to_remove <- c(to_remove, smallest)
  }
  toc()
  if (length(to_remove)) cropLay <- cropLay[-unique(to_remove), ]
  
  # 5a) Ensure all geometries are valid before clustering. Invalid rings coming
  #      out of the intersect/overlap pruning step can later trigger GEOS/JTS
  #      shell assertions when we union clusters (observed in eccc_parallel_cluster).
  #      We first try terra::makeValid for a cheap repair; if that still leaves
  #      issues fall back to sf/lwgeom and optionally simplify slivers.
  if (any(!terra::is.valid(cropLay))) {
    cropLay <- tryCatch({
      terra::makeValid(cropLay)
    }, error = function(e) cropLay)
  }
  if (any(!terra::is.valid(cropLay))) {
    cropLay_sf <- sf::st_as_sf(cropLay)
    make_valid <- function(x) {
      if (requireNamespace("lwgeom", quietly = TRUE)) {
        return(suppressWarnings(lwgeom::st_make_valid(x)))
      }
      suppressWarnings(sf::st_make_valid(x))
    }
    cropLay_sf <- make_valid(cropLay_sf)
    # When make_valid splits multipart geometries, retain the linework components
    cropLay_sf <- sf::st_cast(cropLay_sf, "LINESTRING", warn = FALSE)
    # Drop empty geometries that can be produced by collection extraction
    cropLay_sf <- cropLay_sf[!sf::st_is_empty(cropLay_sf), , drop = FALSE]
    # If slivers remain invalid, a light simplify can remove self-overlap spikes
    if (any(!sf::st_is_valid(cropLay_sf))) {
      cropLay_sf <- sf::st_simplify(cropLay_sf, dTolerance = 10, preserveTopology = TRUE)
      cropLay_sf <- make_valid(cropLay_sf)
      cropLay_sf <- cropLay_sf[!sf::st_is_empty(cropLay_sf), , drop = FALSE]
    }
    if (!inherits(cropLay_sf, "SpatVector")) {
      cropLay <- terra::vect(cropLay_sf)
    } else {
      cropLay <- cropLay_sf
    }
  }
  if (nrow(cropLay) == 0L) {
    warning("createCropLayFinalYear1: geometry repair left no features — returning empty.")
    return(list(lines = Lay[0, ], availableArea = finalPotLay))
  }
  
  # 6) cluster by Potential as before
  pots <- unique(cropLay$Potential)
  if (length(pots) == 0L) {
    warning("createCropLayFinalYear1: no valid Potential values after preprocessing — returning empty.")
    return(list(lines = Lay[0, ], availableArea = finalPotLay))
  }
  
  cropLayFinal <- do.call(rbind, lapply(pots, function(Pot) {
    Cache(clusterLines,
          cropLay[cropLay$Potential == Pot, ],
          distThreshold  = clusterDistance,
          currPotential  = Pot,
          runInParallel  = runClusteringInParallel,
          totPotential   = length(pots),
          userTags       = studyAreaHash)
  }))
  
  # Drop degenerate clusters (e.g., single-line clusters) that can destabilize downstream refinement
  if ("Pot_Clus" %in% names(cropLayFinal)) {
    cl_counts <- table(cropLayFinal$Pot_Clus)
    bad_clusters <- names(cl_counts)[cl_counts < 2]
    if (length(bad_clusters)) {
      cropLayFinal <- cropLayFinal[!(cropLayFinal$Pot_Clus %in% bad_clusters), ]
      if (nrow(cropLayFinal) == 0L) {
        warning("createCropLayFinalYear1: all clusters dropped as degenerate (<2 lines) — returning empty.")
        return(list(lines = Lay[0, ], availableArea = finalPotLay))
      }
    }
  }
  
  # 7) consistent list return for updated generator
  return(list(
    lines         = cropLayFinal,
    availableArea = finalPotLay
  ))
}
