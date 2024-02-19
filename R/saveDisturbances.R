saveDisturbances <- function(disturbanceList, 
                             currentTime, 
                             overwrite,
                             runName){
  lapply(names(disturbanceList), function(Sector){
    whichLay <- names(disturbanceList[[Sector]])[!grepl(x = names(disturbanceList[[Sector]]), 
                                                        pattern = "potential")]
    if (class(disturbanceList[[Sector]]) %in% c("SpatVector", 
                                                    "sf",
                                                    "RasterLayer",
                                                    "SpatRaster")){
      message(paste0("Layer: ", Sector, " was likely not generated (i.e., lack of potential ",
                     "disturbance layer). Saving current disturbance layer (i.e., this layer may ",
                     "not vary through time!)"))
      Lay <- disturbanceList[[Sector]]
      if (any(class(Lay) %in% c("SpatVector", "sf"),
              is(Lay, "Spatial"))){
        if (!is(Lay, "SpatVector"))
          Lay <- vect(Lay)
        terra::writeVector(x = Lay,  
                           filename = file.path(Paths[['outputPath']],
                                                paste0("disturbances_", Sector, "_", 
                                                       currentTime, "_",
                                                       runName, ".shp")), 
                           filetype = "ESRI Shapefile", overwrite = overwrite)
      } else {
        if (class(Lay) %in% c("RasterLayer", "SpatRaster")) {
          if (is(Lay, "RasterLayer"))
            Lay <- rast(Lay)
          terra::writeRaster(x = Lay, filename = file.path(Paths[['outputPath']],
                                                           paste0("disturbances_", Sector, "_", 
                                                                  currentTime, "_",
                                                                  runName, ".tif")), 
                             filetype = "GTiff", overwrite = overwrite)
        } else {
          if (is.null(Lay)){
            warning(paste0("The layer for ",Sector," -- ",LAYER, "is NULL. Not saving."), 
                    immediate. = TRUE)
          } else stop(paste0("Objects of class ", class(Lay), "can't be used in the module. ",
                             "Please make sure all spatial objects are of raster, sp, sf or terra format."))
        }
      }
    } else {
    lay <- lapply(whichLay, function(LAYER){
        Lay <- disturbanceList[[Sector]][[LAYER]]
        message(paste0("Saving layer: ", Sector, " -- ", LAYER))
        if (any(class(Lay) %in% c("SpatVector", "sf"),
                is(Lay, "Spatial"))){
          if (!is(Lay, "SpatVector"))
            Lay <- vect(Lay)
          terra::writeVector(x = Lay,  
                             filename = file.path(Paths[['outputPath']],
                                                  paste0("disturbances_", Sector, "_", LAYER, "_", 
                                                         currentTime, "_",
                                                         runName, ".shp")), 
                             filetype = "ESRI Shapefile", overwrite = overwrite)
        } else {
          if (class(Lay) %in% c("RasterLayer", "SpatRaster")) {
            if (is(Lay, "RasterLayer"))
              Lay <- rast(Lay)
            terra::writeRaster(x = Lay, filename = file.path(Paths[['outputPath']],
                                                             paste0("disturbances_", Sector, "_", 
                                                                    LAYER, "_", currentTime, "_",
                                                                    runName, ".tif")), 
                               filetype = "GTiff", overwrite = overwrite)
          } else {
            if (is.null(Lay)){
              warning(paste0("The layer for ",Sector," -- ",LAYER, "is NULL. Not saving."), 
                      immediate. = TRUE)
            } else stop(paste0("Objects of class ", class(Lay), "can't be used in the module. ",
                               "Please make sure all spatial objects are of raster, sp, sf or terra format."))
          }
        }
    })
    }
  })
  message(crayon::green(paste0("All disturbances saved for ", currentTime)))
}