saveDisturbances <- function(disturbanceList, 
                             currentTime){
  lapply(names(disturbanceList), function(Sector){
    whichLay <- names(disturbanceList[[Sector]])[!grepl(x = names(disturbanceList[[Sector]]), 
                                                        pattern = "potential")]
    lay <- lapply(whichLay, function(LAYER){
      Lay <- disturbanceList[[Sector]][[LAYER]]
      Lay <- as(Lay, "Spatial")
      writeOGR(obj = Lay, dsn = Paths[['outputPath']], 
               layer = paste0("disturbances_", Sector, "_", LAYER, "_", currentTime), 
               driver = "ESRI Shapefile", overwrite = overwrite) #also you were missing the driver argument
    })
  })
  message(crayon::green(paste0("Disturbances saved for ", currentTime)))
}