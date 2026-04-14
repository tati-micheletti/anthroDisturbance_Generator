createBufferedDisturbances <- function(disturbanceList,
                                       bufferSize,
                                       rasterToMatch,
                                       studyArea,
                                       currentTime,
                                       convertToRaster){
  
  # Keep both a SpatRaster view (for terra ops) and a RasterLayer (for fasterize)
  targetSR <- terra::rast(rasterToMatch)                 ### FIX: canonical SpatRaster view
  if (inherits(rasterToMatch, "SpatRaster")) {
    rasterToMatchR <- raster::raster(rasterToMatch)      # for fasterize template
  } else {
    rasterToMatchR <- rasterToMatch                      # already RasterLayer
  }
  
  studyAreaTotalArea <- terra::expanse(studyArea, transform = FALSE, unit = "km")
  cell_km2 <- prod(terra::res(targetSR)) / 1e6           ### FIX: real cell area
  
  unlistedDL <- unlist(disturbanceList)
  nonPotentialLayers <- names(unlistedDL)[!grepl(x = names(unlistedDL), "potential")]
  
  allBuffered <- lapply(nonPotentialLayers, function(INDEX) {
    lay <- unlistedDL[[INDEX]]
    whichClass <- class(lay)
    tic(paste0("Time elapsed for buffering ", INDEX, ":"))
    
    if (whichClass %in% c("SpatRaster", "RasterLayer")) {
      if (inherits(lay, "RasterLayer")) lay <- terra::rast(lay)
      
      message(paste0("Polygonizing raster ", INDEX))
      if (sum(lay[], na.rm = TRUE) == 0) {
        message(paste0("Layer ", INDEX, " seems empty. Returning NULL."))
        return(NULL)
      }
      lay[lay != 1] <- NA
      lay <- terra::as.polygons(lay, values = TRUE, na.rm = TRUE, digits = 5)
      lay <- terra::subset(x = lay, subset = lay[[names(lay)]] == 1)
    }
    
    message(paste0("Buffering vector ", INDEX))
    bLay <- terra::buffer(lay, width = bufferSize)
    bLay <- terra::aggregate(bLay, dissolve = TRUE)
    
    # early-outs for empties
    if (nrow(bLay) == 0) {
      message(paste0("The layer for ", INDEX, " is empty. Returning NULL"))
      return(NULL)
    }
    
    totalAreaKm2 <- terra::expanse(bLay, transform = FALSE, unit = "km")
    message(crayon::green(paste0("Total area disturbed for ", INDEX, " as vector: ",
                                 round(totalAreaKm2, 3), " km2")))
    message(paste0("This represents ",
                   round(100 * (totalAreaKm2 / studyAreaTotalArea), 4), "% of total area."))
    
    # Rasterize against the RasterLayer template (fasterize requirement)
    bLaySF  <- sf::st_as_sf(bLay)
    bLayRas <- fasterize::fasterize(sf = bLaySF, raster = rasterToMatchR)
    names(bLayRas) <- INDEX
    
    # 1/0 semantics & area message
    distTb <- table(bLayRas[])
    n1     <- if ("1" %in% names(distTb)) as.numeric(distTb[["1"]]) else 0
    totDistKm <- n1 * cell_km2
    message(crayon::green(paste0("Total area disturbed for ", INDEX,
                                 " as raster: ", round(totDistKm, 3), " km2")))
    
    if (convertToRaster) bLayRas else bLay
  })
  
  if (convertToRaster){
    # Drop NULLs
    allBuffered <- Filter(Negate(is.null), allBuffered)
    
    # Nothing? return clean 0/1 mask with template geometry
    if (!length(allBuffered)){
      if (inherits(rasterToMatch, "RasterLayer")) {
        out <- rasterToMatch
        out[!is.na(out[])] <- 0
        return(out)
      } else {
        out <- targetSR
        out[!is.na(out[])] <- 0
        return(out)
      }
    }
    
    # Align every layer to targetSR (SpatRaster view) and force 0/1
    allBuffered <- lapply(allBuffered, function(r) {
      r <- terra::rast(r)  # RasterLayer -> SpatRaster
      if (!terra::compareGeom(r, targetSR, stopOnError = FALSE)) {
        r <- terra::project(r, targetSR)
        r <- terra::resample(r, targetSR, method = "near")
      }
      terra::ifel(r == 1, 1, 0)  # binary
    })
    
    # Stack and merge using a stable reducer (no custom fun)
    st <- do.call(c, allBuffered)
    merged <- if (terra::nlyr(st) == 1L) st[[1]] else terra::app(st, fun = max, na.rm = TRUE)
    
    # Impose template NA mask without evaluating !is.na(targetSR)
    outSR <- tryCatch(
      terra::mask(merged, targetSR),   # preferred (keeps NA where targetSR is NA)
      error = function(e) merged       # fallback if targetSR has "no values"
    )
    
    # Summary (0/1/NA semantics preserved)
    vals <- terra::values(outSR, mat = FALSE)
    n1   <- sum(vals == 1, na.rm = TRUE)
    n0   <- sum(vals == 0, na.rm = TRUE)
    totDistKm <- n1 * cell_km2
    totalKm2  <- (n0 + n1) * cell_km2
    percDist  <- if ((n0 + n1) > 0) round(100 * (n1 / (n0 + n1)), 2) else 0
    message(paste0("Buffered disturbances for all disturbances. Total area disturbed: ",
                   round(totDistKm, 3), " km2 -- ", percDist,
                   "% of the total area (", round(totalKm2, 3), " km2)"))
    
    # Class parity with template
    if (inherits(rasterToMatch, "RasterLayer")) {
      return(raster::raster(outSR))
    } else {
      return(outSR)
    }
    
  } else {
    # (vector branch unchanged)
    if (any(vapply(allBuffered, is.null, FALSE))) {
      allBuffered <- allBuffered[!vapply(allBuffered, is.null, FALSE)]
    }
    disturbanceLayerAll <- do.call(rbind, allBuffered)
    disturbanceLayerAll <- terra::aggregate(disturbanceLayerAll, dissolve = TRUE)
    
    totalDisturbedArea <- terra::expanse(disturbanceLayerAll, transform = FALSE, unit = "km")
    message(paste0("Buffered disturbances for all disturbances. Total area disturbed: ",
                   round(totalDisturbedArea, 2), " km2 -- ",
                   round(100 * (totalDisturbedArea / studyAreaTotalArea), 2),
                   "% of the total area (", round(studyAreaTotalArea, 2), " km2)"))
    return(disturbanceLayerAll)
  }
}
