createBufferedDisturbances <- function(disturbanceList, 
                                       bufferSize,
                                       rasterToMatch,
                                       studyArea,
                                       currentTime,
                                       convertToRaster){
  if (is(rasterToMatch, "SpatRaster"))
    rasterToMatchR <- raster::raster(rasterToMatch) else
      rasterToMatchR <- rasterToMatch
    
    studyAreaTotalArea <- terra::expanse(studyArea, transform = FALSE, unit = "km")
    if (FALSE){
      rtmTotalArea <- as.numeric(table(rasterToMatch[])*prod(res(rasterToMatch)/1000)) # how many pixels * pixel area in km2
      message(paste0("The difference in total area between RasterToMatch and StudyArea is ", 
                     round((1-(studyAreaTotalArea/rtmTotalArea))*100, 2), "%"))
    }
    # Make it generic
    unlistedDL <- unlist(disturbanceList)
    # Exclude any "potential" disturbances!!
    nonPotentialLayers <- names(unlistedDL)[!grepl(x = names(unlistedDL), "potential")]
    allBuffered <- lapply(nonPotentialLayers, function(INDEX){
      # Buffer all shapefiles
      lay <- unlistedDL[[INDEX]]
      whichClass <- class(lay)
      tic(paste0("Time elapsed for buffering ", INDEX, ":"))
      if (whichClass %in% c("SpatRaster", "RasterLayer")){
        if (is(lay, "RasterLayer"))
          lay <- rast(lay)
        message(paste0("Polygonizing raster ", INDEX))
        if (sum(lay[], na.rm = TRUE) == 0){
          message(paste0("Layer ", INDEX, " seems empty. Returning NULL."))
          return(NULL) 
        }
        lay[lay != 1] <- NA
        lay <- as.polygons(lay, values=TRUE, na.rm=TRUE, digits=5)
        # Need to remove the study area polygon (value 0)
        lay <- terra::subset(x = lay, subset = lay[[names(lay)]] == 1)
      }
      message(paste0("Buffering vector ", INDEX))
      # Probably need to aggregate here!
      bLay <- terra::buffer(lay, width = bufferSize)
      bLay <- terra::aggregate(bLay, dissolve = TRUE)
      totalAreaKm2 <- terra::expanse(bLay, transform = FALSE, unit = "km")
      message(crayon::green(paste0("Total area disturbed for ", INDEX, " as vector: ", round(totalAreaKm2, digits = 3), " km2")))
      message(paste0("This represents ", round(100*(totalAreaKm2/studyAreaTotalArea), 4), "% of total area."))
      if (length(bLay) == 0){
        message(paste0("The layer for ", INDEX, " is empty. Returning NULL"))
        return(NULL)
      }
      # Calculating as Raster
        bLaySF <- sf::st_as_sf(bLay)
        bLayRas <- fasterize::fasterize(sf = bLaySF, raster = rasterToMatchR)
        names(bLayRas) <- INDEX
        distTb <- table(bLayRas[])
        totDistKm <- distTb["1"] * 0.0625
        message(crayon::green(paste0("Total area disturbed for ", INDEX,
                                     " as raster: ", totDistKm, " km2")))
        if (convertToRaster){
          return(bLayRas)          
        } else {
          return(bLay)
        }
    })
    if (convertToRaster){
      # Convert all pixels that are disturbed to 1, and non-disturbed to 0
      bufferedAnthropogenicDisturbance500m <- rasterToMatch
      bufferedAnthropogenicDisturbance500m[!is.na(bufferedAnthropogenicDisturbance500m[])] <- 0
      for (i in 1:length(allBuffered)){
        bufferedAnthropogenicDisturbance500m[which(allBuffered[[i]][] == 1)] <- 1
      }
      distTable <- table(bufferedAnthropogenicDisturbance500m[])
      totDistKm <- distTable["1"]*0.0625
      percDist <- round(100*(distTable["1"]/(distTable["0"]+distTable["1"])), 2)
      message(paste0("Buffered disturbances for all disturbances. Total area disturbed: ", 
                     totDistKm, "km2 -- ", percDist, "% of the total area (", 
                     0.0625*(distTable["0"]+distTable["1"]),"km2)"))
      return(bufferedAnthropogenicDisturbance500m)
    } else {
      if (any(unlist(lapply(allBuffered, is.null)))){
        allBufferedNoNull <- allBuffered[-which(sapply(allBuffered, is.null))]# Remove any NULL layers
        disturbanceLayerAll <- do.call(rbind, allBufferedNoNull)
      } else {
        disturbanceLayerAll <- do.call(rbind, allBuffered)
      }
      disturbanceLayerAll <- terra::aggregate(disturbanceLayerAll, dissolve = TRUE)
      totalDisturbedArea <- terra::expanse(disturbanceLayerAll, transform = FALSE, unit = "km")
      message(paste0("Buffered disturbances for all disturbances. Total area disturbed: ", 
                     round(totalDisturbedArea, 2), " km2 -- ", 
                     round(100*(totalDisturbedArea/studyAreaTotalArea), 2), "% of the total area (", 
                     round(studyAreaTotalArea, 2)," km2)"))
      return(disturbanceLayerAll)
    }
}
