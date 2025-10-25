calculateRate <- function(disturbanceParameters,
                          disturbanceDT,
                          disturbanceList,
                          whichToUpdate,
                          RTM,
                          diffYears = "2010_2015",
                          destinationPath,
                          #overwriteDisturbanceLayersNEW, # currently not in use
                          #overwriteDisturbanceLayersOLD,
                          studyArea,
                          DisturbanceRate,
                          totalDisturbanceRate,
                          disturbanceRateRelatesToBufferedArea,
                          maskOutLinesFromPolys,
                          aggregateSameDisturbances,
                          archiveNEW,
                          targetFileNEW,
                          urlNEW,
                          archiveOLD,
                          targetFileOLD,
                          urlOLD){
  if (length(whichToUpdate) == 0) {
    message("calculateRate: no rows need updating, returning original table.")
    return(disturbanceParameters)
  }
  
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
  if (diffYears == "OLD_2015")
    warning(paste0("While ECCC OLD footprint layer covers the whole country, ECCC 2015",
                   " covers only the caribou ranges. This means that the 2015 layer will",
                   " return a biased value for disturbance than it actually is IF the study",
                   " area extends outside of a caribou range. Ideally, either 1) this layer ",
                   " should be replaced by a similar layer which instead covers the ",
                   "whole country, or 2) one should provide the object 'DisturbanceRate' to",
                   " simulate development scenarios focusing on individual disturbances."), 
            immediate. = TRUE)
    
    AD_changed_file <- file.path(destinationPath, 
                                 paste0("anthropogenicDisturbance_ECCC_", 
                                        diffYears,"_", digest(studyArea),".csv"))
    if (!file.exists(AD_changed_file)){
      distECCC <- disturbanceInfoFromECCC(
        studyArea                     = studyArea, 
        RTM                           = RTM,
        disturbanceList               = disturbanceList,
        totalstudyAreaVAreaSqKm       = totalstudyAreaVAreaSqKm,
        classesAvailable              = classesAvailable,
        destinationPath               = destinationPath,
        bufferedDisturbances          = disturbanceRateRelatesToBufferedArea,
        maskOutLinesFromPolys         = maskOutLinesFromPolys,
        aggregateSameDisturbances     = aggregateSameDisturbances,
        #overwriteDisturbanceLayersNEW = overwriteDisturbanceLayersNEW,
        #overwriteDisturbanceLayersOLD = overwriteDisturbanceLayersOLD,
        archiveNEW                    = archiveNEW,
        diffYears                     = diffYears,
        targetFileNEW                 = targetFileNEW,
        urlNEW                        = urlNEW,
        archiveOLD                    = archiveOLD,
        targetFileOLD                 = targetFileOLD,
        urlOLD                        = urlOLD
      )
      
      AD_changed <- distECCC[["AD_changed"]]
    } else {
      AD_changed <- data.table::fread(AD_changed_file)
    }
    if (nrow(AD_changed[Class %in% toLookFor$classToSearch]) == 0) {
      warning("No historical ECCC disturbance for these classes → setting rate to zero")
      # set disturbanceRate = 0 for all rows we intended to update
      disturbanceParameters[whichToUpdate, disturbanceRate := 0]
      return(disturbanceParameters)
    }
    
    toFill <- AD_changed[Class %in% toLookFor[["classToSearch"]]]
    toUse <- merge(toFill, toLookFor[, c("classToSearch", "dataClass")], all.x = TRUE, 
                   by.x = "Class", by.y = "classToSearch")
    toUse <- data.table::dcast(toUse, dataClass ~ ., fun.agg = sum, 
                               value.var = c("yearOLD", "yearNEW"))
    # Need to recalculate proportions, though!
    toUse[, disturbProportionInAreaOLD := yearOLD/totalstudyAreaVAreaSqKm]
    toUse[, disturbProportionInAreaNEW := yearNEW/totalstudyAreaVAreaSqKm]
    toUse[, totalProportionAreaDisturbedOLD := sum(disturbProportionInAreaOLD)]
    toUse[, totalProportionAreaDisturbedNEW := sum(disturbProportionInAreaNEW)]
    
    toUse <- toUse[, c("dataClass", "disturbProportionInAreaOLD", "disturbProportionInAreaNEW")]
    toUse[, proportionAreaSqKmChangedPerYear := (disturbProportionInAreaNEW-disturbProportionInAreaOLD)/5]
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
      if (any(is.null(lay), length(lay) == 0)){
        # If no current origin layer exists in the AOI, but the user supplied
        # DisturbanceRate, proceed using that rate (generation will rely on potential).
        if (!is.null(DisturbanceRate)) {
          message(crayon::yellow(paste0(
            "No current '", sub[["disturbanceOrigin"]], "' layer in AOI; proceeding with supplied DisturbanceRate for ",
            sub[["dataClass"]], ".")))
          updatedVal <- DisturbanceRate[dataName == sub[["dataName"]] & 
                                          dataClass == sub[["dataClass"]] & 
                                          disturbanceType == sub[["disturbanceType"]] & 
                                          disturbanceOrigin == sub[["disturbanceOrigin"]], "disturbanceRate"]
          sub[, disturbanceRate := updatedVal]
        } else {
          message(crayon::red(paste0("There is no potential for ", sub[["dataClass"]],
                                     " in the study Area. The module will return NULL and this",
                                     "class will not be simulated. ",
                                     "If this is wrong, please provide a layer with this information ")))
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
                                   paste0("proportionTable_ECCC_", diffYears, "_", 
                                          digest(studyArea),".csv"))
            if (!file.exists(prop_file)) {
              if (!exists("distECCC")) {
                distECCC <- disturbanceInfoFromECCC(
                  studyArea                     = studyArea, 
                  RTM                           = RTM,
                  disturbanceList               = disturbanceList,
                  totalstudyAreaVAreaSqKm       = totalstudyAreaVAreaSqKm,
                  classesAvailable              = classesAvailable,
                  destinationPath               = destinationPath,
                  bufferedDisturbances          = disturbanceRateRelatesToBufferedArea,
                  maskOutLinesFromPolys         = maskOutLinesFromPolys,
                  aggregateSameDisturbances     = aggregateSameDisturbances,
                  #overwriteDisturbanceLayersNEW = overwriteDisturbanceLayersNEW,
                  #overwriteDisturbanceLayersOLD = overwriteDisturbanceLayersOLD,
                  archiveNEW                    = archiveNEW,
                  diffYears                     = diffYears,
                  targetFileNEW                 = targetFileNEW,
                  urlNEW                        = urlNEW,
                  archiveOLD                    = archiveOLD,
                  targetFileOLD                 = targetFileOLD,
                  urlOLD                        = urlOLD
                )
              }
              proportionTable <- distECCC[["proportionTable"]]
            } else {
              proportionTable <- data.table::fread(prop_file)
            }
            # 1. Determine which disturbances are needed
            neededDistRates <- sub[["disturbanceOrigin"]]
            # 2. Match which of these have a value in proportionTable
            subProportionTable <- proportionTable[dataClass %in% neededDistRates,]
            # 2b. If the disturbance exists but was negative/absent, keep the row with zero rate.
            if (nrow(subProportionTable) == 0) {
              message(paste0("The disturbance ", sub[["dataName"]], " of the sector ",
                             sub[["disturbanceOrigin"]], " was virtually zero or negative in the area.",
                             " Setting rate to 0 and keeping row for scheduling/diagnostics."))
              sub[, disturbanceRate := 0]
              return(sub)
            }
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
            toUseSub <- toUse[dataClass %in% sub[["disturbanceOrigin"]]]
            if (NROW(toUseSub) == 0L) {
              warning("No ECCC-derived rate for ", sub$dataName, "/", sub$dataClass, 
                      "; setting rate to zero", immediate. = TRUE)
              sub[, disturbanceRate := 0]
            } else {
              if (NROW(toUseSub) > 1L) {
                warning("Repeated ECCC rows for ", sub$dataName, "/", sub$dataClass, "; using first")
                toUseSub <- toUseSub[1]
              }
              # Pull the single value, sanitize it
              val <- toUseSub[["proportionAreaSqKmChangedPerYear"]]
              if (length(val) == 0L || is.na(val)) val <- 0
              # Already clipped negatives to 0 earlier, but be defensive:
              val <- max(0, val)
              sub[, disturbanceRate := 100 * val]
              
              message(paste0("Using yearly disturbance rate for ", sub[["dataName"]],
                             " as ", round(sub[["disturbanceRate"]], 6),
                             "% of the total area."))
            }
          }
        } else {
          # One can also pass the parameter totalDisturbanceRate. If so, we use the ECCC data to calculate
          # the % each disturbance needs to be to achieve (more or less) the total expected disturbance 
          subDR <- DisturbanceRate[
            dataName == sub$dataName &
              dataClass == sub$dataClass &
              disturbanceType == sub$disturbanceType &
              disturbanceOrigin == sub$disturbanceOrigin
          ]
          if (nrow(subDR) == 0) {
            warning("No matching row in DisturbanceRate for ",
                    sub$dataName, "/", sub$dataClass, "; setting rate to zero",
                    immediate. = TRUE)
            sub[, disturbanceRate := 0]
            return(sub)
          }
          if (nrow(subDR) > 1) {
            warning("Multiple matches in DisturbanceRate for ", 
                    sub$dataName, "/", sub$dataClass, "; using first")
          }
          updatedVal <- subDR[1, disturbanceRate]
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
