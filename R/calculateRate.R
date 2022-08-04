calculateRate <- function(disturbanceParameters,
                          disturbanceList,
                          whichToUpdate,
                          RTM){
  
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
        message(crayon::yellow(paste0("There is no information on rate for ", sub[["dataClass"]],
                                      ". However, this is a potentialWindTurbines. The module will ",
                                      "return a rate of one turbine per 10 years, or 0.1 turbines ",
                                      "per year. The 0.1 pixels correspond to ", round(100*(0.1/sum(RTM[], na.rm = TRUE)), 3), 
                                      "% of the total area, being this the value replacing NA. ",
                                      "If this is wrong, please provide both disturbanceRate and ",
                                      "disturbanceSize in the disturbanceParameters table.")))
        
        sub[, disturbanceRate := 0.1 / sum(RTM[], na.rm = TRUE)] 
        # 0.1 pixels over the total area in pixels, we have the rate
      } else {
        message(crayon::red(paste0("There is no information on rate for ", sub[["dataClass"]],
                                   ". The module will return NULL and this class will not be simulated. ",
                                   "If this is wrong, please provide both disturbanceRate and disturbanceSize ",
                                   "in the disturbanceParameters table.")))
        sub <- NULL
      }
    } else {
      if (sub[["disturbanceType"]] == "Enlarging"){
        # If Enlarging, we can consider the 0.2% of the total area a year as a default
        sub[, disturbanceRate := (0.2/100)]
      } else {
        # If Generating, we can consider the 0.2% of the total area a year as a default
        # Not all data has years, some are in hard formats to convert... here one can 
        # later modify to include specific rates derived from data. For now, we will go 
        # with a default of 0.2% of the current disturbed area, as suggested in 2019 Species at Risk Act 
        # Conservation Agreement for the Conservation of the Boreal Caribou
        sub[, disturbanceRate := (0.2/100)]
        }
    }
    return(sub)
  }))
  
  disturbanceParametersUp <- rbind(disturbanceParameters[!whichToUpdate, ], 
                                   updatedDisturbanceParameters, use.names = TRUE)
  return(disturbanceParametersUp)
  
}