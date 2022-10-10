saveDisturbances <- function(disturbanceList, 
                             currentTime, 
                             overwrite,
                             runName){
  lapply(names(disturbanceList), function(Sector){
    whichLay <- names(disturbanceList[[Sector]])[!grepl(x = names(disturbanceList[[Sector]]), 
                                                        pattern = "potential")]
    lay <- lapply(whichLay, function(LAYER){
      Lay <- disturbanceList[[Sector]][[LAYER]]
      if (any(class(Lay) %in% c("SpatVector", "sf"),
              is(Lay, "Spatial"))){
        Lay <- as(Lay, "Spatial")
        writeOGR(obj = Lay, dsn = Paths[['outputPath']], 
                 layer = paste0("disturbances_", Sector, "_", LAYER, "_", currentTime, "_",
                                runName), 
                 driver = "ESRI Shapefile", overwrite = overwrite)
      } else {
        if (class(Lay) %in% c("RasterLayer", "SpatRaster")) {
          Lay[] <- Lay[]
          writeRaster(x = Lay, filename = file.path(Paths[['outputPath']],
                                                    paste0("disturbances_", Sector, "_", 
                                                           LAYER, "_", currentTime, "_",
                                                           runName)), 
                      format = "GTiff", overwrite = overwrite)
        } else {
          stop(paste0("Objects of class ", class(Lay), "can't be used in the module. ",
                      "Please make sure all spatial objects are of raster, sp, sf or terra format."))
        }
      }
    })
  })
  message(crayon::green(paste0("Disturbances saved for ", currentTime)))
}