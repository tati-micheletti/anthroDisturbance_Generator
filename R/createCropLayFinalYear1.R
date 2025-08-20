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
  overlap_matrix <- terra::relate(x = cropLay, relation = "T********", pairs = TRUE)
  overlapMatrix  <- data.table::as.data.table(overlap_matrix)
  overlapMatrix[, self := data.table::fifelse(id.1 == id.2, TRUE, FALSE)]
  overlapMatrix  <- overlapMatrix[self == FALSE][, self := NULL]
  # normalize (id.1 < id.2) and unique
  overlapMatrix  <- unique(overlapMatrix[, .(id.1 = pmin(id.1, id.2),
                                             id.2 = pmax(id.1, id.2))])
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
  
  # 7) consistent list return for updated generator
  return(list(
    lines         = cropLayFinal,
    availableArea = finalPotLay
  ))
}
