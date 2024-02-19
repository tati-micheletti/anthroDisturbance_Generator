calculateSize <- function(disturbanceParameters,
                          disturbanceList,
                          whichToUpdate){

  dP <- disturbanceParameters[whichToUpdate, ]
  updatedDisturbanceParameters <- rbindlist(lapply(1:NROW(dP), function(INDEX){
    # Get the row to be calculated
    sub <- dP[INDEX, ]
    # Get the corresponding layer
    lay <- disturbanceList[[sub[["dataName"]]]][[sub[["disturbanceOrigin"]]]]
    # If the layer was empty, it will return null. That means we don't have data to calculate the 
    # size. In this case, we will use defaults for what we know: 
    if (is.null(lay)){
      # windTurbines: 6070m2 per 2-megawatt of energy. In the NWT the only big turbine is in a
      #               Diamond exploration, with 9.2MW being generated (27,900 m2). Each pixel is
      #               62,500m2, so we will have in average, one big turbine in 2.3 pixels.
      #               There is one turbine being constructed and should be done in 2023 (2021 to be 
      #               easier on the simulations) -- Size: 3.5-megawatt = less than a pixel
      #               There are current no plans for more turbines. Still, we might be able to add
      #               one medium turbine every 10 years. That would mean 1 pixel every 10 years. 
      # There are no other values that need to be inputted. If the lay is not windTurbine, set the
      # size to NULL and basically skip it with a warning
      if (sub[["dataClass"]] == "potentialWindTurbines"){
        sub[, disturbanceSize := 62500]
        message(crayon::yellow(paste0("There is no information on size for ", sub[["dataClass"]],
                                   ". However, this is a potentialWindTurbines. The module will ",
                                   "return a size of 62,500m2 which is equivalent to a medium turbine.",
                                   "If this is wrong, please provide both disturbanceRate and disturbanceSize ",
                                   "in the disturbanceParameters table.")))
      } else {
        message(crayon::red(paste0("There is no information on size for ", sub[["dataClass"]],
                       ". The module will return NULL and this class will not be simulated. ",
                       "If this is wrong, please provide both disturbanceRate and disturbanceSize ",
                       "in the disturbanceParameters table.")))
        sub <- NULL
      }
    } else {
      if (class(lay) %in% c("RasterLayer", "SpatRaster")) {
        message(paste0("The layer ", sub[["dataClass"]], " is not a vector. Trying to convert."))
        if (lay %in% "RasterLayer") lay <- rast(lay)
        lay <- as.polygons(lay, dissolve = FALSE)
      }
      allAreas <- terra::expanse(lay, transform = FALSE, unit = "m")
      sub[, disturbanceSize := paste0("rtnorm(1, ", round(mean(allAreas), 2), ", ", round(sd(allAreas), 2), ", lower = 0)")]
      }
    return(sub)
  }))

  disturbanceParametersUp <- rbind(disturbanceParameters[!whichToUpdate, ], 
                                   updatedDisturbanceParameters, use.names = TRUE)
  return(disturbanceParametersUp)
  
}