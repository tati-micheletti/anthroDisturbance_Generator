calculateRate <- function(disturbanceParameters,
                          disturbanceDT,
                          disturbanceList,
                          whichToUpdate,
                          RTM,
                          destinationPath,
                          overwriteDisturbanceLayers2015,
                          overwriteDisturbanceLayers2010,
                          studyArea,
                          DisturbanceRate,
                          totalDisturbanceRate,
                          disturbanceRateRelatesToBufferedArea,
                          maskOutLinesFromPolys,
                          aggregateSameDisturbances){
  
  if (all(!is.null(DisturbanceRate), !is.null(totalDisturbanceRate)))
    stop("Both DisturbanceRate and totalDisturbanceRate were provided. Please provide only one,",
         "or none.")
  
  # Subsetting parameters to update
  dP <- disturbanceParameters[whichToUpdate, ]
  classesAvailable <- unique(disturbanceDT[fieldToSearch == "Class", c("classToSearch", "dataName", "dataClass")])
  toLookFor <- classesAvailable[dataClass %in% dP[["disturbanceOrigin"]]]
  
  # Calculating Total study area size
  if (!is(studyArea, "SpatVector"))
    studyAreaV <- terra::vect(studyArea) else studyAreaV <- studyArea
  studyAreaV <- terra::project(x = studyAreaV, y = terra::crs(RTM))
  uniStudyArea <- terra::aggregate(studyAreaV)
  totalstudyAreaVAreaSqKm <- terra::expanse(uniStudyArea, unit = "km", transform = FALSE)
  
  if (is.null(DisturbanceRate)){ # When DisturbanceRate is NOT provided

    warning(paste0("While ECCC 2010 footprint layer covers the whole country, ECCC 2015",
                   " covers only the caribou ranges. This means that the 2015 layer will",
                   " return a biased value for disturbance than it actually is IF the study",
                   " area extends outside of a caribou range. Ideally, either 1) this layer ",
                   " should be replaced by a similar layer which instead covers the ",
                   "whole country, or 2) one should provide the object 'DisturbanceRate' to",
                   " simulate development scenarios."), 
            immediate. = TRUE)
    
    AD_changed_file <- file.path(destinationPath, 
                                 "anthropogenicDisturbance_ECCC_2010_2015.csv")
    if (!file.exists(AD_changed_file)){
      distECCC <- disturbanceInfoFromECCC(studyArea = studyArea, 
                                          RTM = RTM,
                                          totalstudyAreaVAreaSqKm = totalstudyAreaVAreaSqKm,
                                          classesAvailable = classesAvailable,
                                          destinationPath = destinationPath,
                                          bufferedDisturbances = disturbanceRateRelatesToBufferedArea,
                                          maskOutLinesFromPolys = maskOutLinesFromPolys,
                                          aggregateSameDisturbances = aggregateSameDisturbances)
      
      AD_changed <- distECCC[["AD_changed"]]
    } else {
      AD_changed <- data.table::fread(AD_changed_file)
    }
    toFill <- AD_changed[Class %in% toLookFor[["classToSearch"]]]
    toUse <- merge(toFill, toLookFor[, c("classToSearch", "dataClass")], all.x = TRUE, 
                   by.x = "Class", by.y = "classToSearch")
    toUse <- data.table::dcast(toUse, dataClass ~ ., fun.agg = sum, 
                             value.var = c("year2010", "year2015"))
    # Need to recalculate proportions, though!
    toUse[, disturbProportionInArea2010 := year2010/totalstudyAreaVAreaSqKm]
    toUse[, disturbProportionInArea2015 := year2015/totalstudyAreaVAreaSqKm]
    toUse[, totalProportionAreaDisturbed2010 := sum(disturbProportionInArea2010)]
    toUse[, totalProportionAreaDisturbed2015 := sum(disturbProportionInArea2015)]
    
    toUse <- toUse[, c("dataClass", "disturbProportionInArea2010", "disturbProportionInArea2015")]
    toUse[, proportionAreaSqKmChangedPerYear := (disturbProportionInArea2015-disturbProportionInArea2010)/5]
    # If any have negative growth, we ignore
    toUse[proportionAreaSqKmChangedPerYear < 0, proportionAreaSqKmChangedPerYear := 0] 
    #TODO implement reduction of disturbances
  }
    updatedDisturbanceParameters <- rbindlist(lapply(1:NROW(dP), function(INDEX){
      # Get the row to be calculated
      sub <- dP[INDEX, ]
      sub[, disturbanceRate := as.numeric(disturbanceRate)]
      # Get the corresponding layer
      lay <- disturbanceList[[sub[["dataName"]]]][[sub[["disturbanceOrigin"]]]]
      # If the layer was empty, it will return null. That means we don't have data to calculate the 
      # size. In this case, we will use defaults for what we know: 
      message(paste0("Calculating rate for ", sub[["disturbanceType"]], 
                     " of ", sub[["dataClass"]]))
      if (is.null(lay)){
        if (sub[["dataClass"]] == "potentialWindTurbines"){
          
          # Only one wind turbines exists in NWT, there isn't a layer pinpointing where. 
          # Therefore, we researched what would this mean in the NWT (NT1_BCR6: 572448 km2):
          
          # 1 Turbine = 0.0625 km2
          # 0.1 Turbine = 0.00625 km2 
          # areaSqKmChangedPerYear <- 0.00625
          
          # windTurbines: 6070m2 per 2-megawatt of energy. In the NWT the only big turbine is in a
          #               Diamond exploration, with 9.2MW being generated (27,900 m2). Each pixel is
          #               62,500m2 (or 0.0625 km2), so we will have in average, one big turbine in 2.3 pixels.
          #               There is one turbine being constructed and should be done in 2023 (2021 to be 
          #               easier on the simulations) -- Size: 3.5-megawatt = less than a pixel
          #               There are current no plans for more turbines. Still, we might be able to add
          #               one medium turbine every 10 years. That would mean 1 pixel every 10 years 
          #               (or 0.1 pixel every year) over 572448 km2. Each pixel is 0.0625 km2, so  
          #               we have (0.1*0.0625)/572448 = 1.0919e-08 km2 change per km2 of area over 1 year. 
          #               If we multiply this by the study area, we end up with the total area change 
          #               for this area per year. But what we need in the table is the % change per 
          #               year over the area. So we need to divide the resulting value by the area 
          #               and then multiply by 100 to get the % change in relation to the area.
                         
          areaChangePerYearInKm2PerKm2 <- 1.0919e-08 # Which is the same as the proportion of any area
          
          message(crayon::yellow(paste0("There is no information on rate for ", sub[["dataClass"]],
                                        ". However, this is a potentialWindTurbines. The module will ",
                                        "return a rate of one turbine per 10 years, or 0.1 turbines ",
                                        "per year over the NT1_BCR6 area (572448 km2). The 0.1",
                                        "pixels (or 0.625 km2) correspond to ", 
                                        format(areaChangePerYearInKm2PerKm2*100, 
                                               scientific = TRUE), 
                                        "% of the current study area, being this the value replacing NA. ",
                                        "If this is wrong, please provide both disturbanceRate and ",
                                        "disturbanceSize in the disturbanceParameters table or in ",
                                        " the DisturbanceRate")))
          sub[, disturbanceRate := areaChangePerYearInKm2PerKm2*100] #areaSqKmChangedPerYear/totalAreaSqKm
        } else {
          message(crayon::red(paste0("There is no information on rate for ", sub[["dataClass"]],
                                     ". The module will return NULL and this class will not be simulated. ",
                                     "If this is wrong, please provide both disturbanceRate and disturbanceSize ",
                                     "in the disturbanceParameters table.")))
          sub <- NULL
        }
      } else {
        if (is.null(DisturbanceRate)) {
          if (!is.null(totalDisturbanceRate)){            
            warning(paste0("The totalDisturbanceRate was supplied as ", totalDisturbanceRate,
                           ". Using example data derived from the Northwest Territories. If other",
                           " rates are desired, please provide DisturbanceRate"),
                    immediate. = TRUE)
            prop_file <- file.path(destinationPath, 
                                         "proportionTable_ECCC_2010_2015.csv")
            if (!file.exists(prop_file)){
              proportionTable <- distECCC[["proportionTable"]] 
            } else {
              proportionTable <- data.table::fread(prop_file)
            }
            # 1. Determine which disturbances are needed
            neededDistRates <- sub[["disturbanceOrigin"]]
            # 2. Match which of these have a value in proportionTable
            subProportionTable <- proportionTable[dataClass %in% neededDistRates,]
            # 3. Calculate the percentage of the provided totalDisturbanceRate that belongs to each
            # class
            subProportionTable[, calculatedDisturbanceProportion := proportionOfTotalDisturbance*totalDisturbanceRate]
            # 4. Replace in DisturbanceRate the disturbanceRate by calculatedDisturbanceProportion for the
            # ones that are on the table
              sub[disturbanceOrigin == subProportionTable[["dataClass"]],
                                  disturbanceRate := subProportionTable[["calculatedDisturbanceProportion"]]]
          } else {
            # 5. Modify the table, which should be an Input.
            # This process ensures we are not double counting disturbances, and that we can calculate the change 
            # in disturbance per class, which we need to simulate the new disturbances coming.
            toUseSub <- toUse[dataClass %in% sub[["disturbanceOrigin"]], ]
            if (NROW(toUseSub) > 1) {
              warning("Repeated data found in the process. Please debug", immediate. = TRUE)
              browser()
            }
            # Here the disturbance rate is in proportion. Needs to be multiplied for %!
            sub[, disturbanceRate := 100*toUseSub[["proportionAreaSqKmChangedPerYear"]]]
            
            message(paste0("Using yearly disturbance rate for ", sub[["dataName"]],
                           " as ", round(sub[["disturbanceRate"]], 6), 
                           "% of the total area."))
          }
        } else {
          # One can also pass the parameter totalDisturbanceRate. If so, we use the ECCC data to calculate
          # the % each disturbance needs to be to achieve (more or less) the total expected disturbance 
          updatedVal <- DisturbanceRate[dataName == sub[["dataName"]] & 
                                          dataClass == sub[["dataClass"]] & 
                                          disturbanceType == sub[["disturbanceType"]] & 
                                          disturbanceOrigin == sub[["disturbanceOrigin"]], "disturbanceRate"]
          sub[, disturbanceRate := updatedVal]
          # Here the disturbance rate is already in percent!
          message(paste0("Using yearly disturbance rate for ", sub[["dataName"]],
                         " as ", round(sub[["disturbanceRate"]], 4), "% of the total area."))
        }
      }
      return(sub)
    }))
    
    disturbanceParametersUp <- rbind(disturbanceParameters[!whichToUpdate, ], 
                                     updatedDisturbanceParameters, use.names = TRUE)
    message("Final disturbances parameters table: ")
    print(disturbanceParametersUp)
    message(paste0("Total expected yearly rate of new human disturbance (excluding roads and not simulated features): ", 
                   round(sum(disturbanceParametersUp[["disturbanceRate"]], na.rm = TRUE), 4), "% of total area (",
                   round(totalstudyAreaVAreaSqKm, 1), " km2)."))
    return(disturbanceParametersUp)
}