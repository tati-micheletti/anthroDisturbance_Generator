generateDisturbancesShp <- function(disturbanceParameters,
                                 disturbanceList,
                                 rasterToMatch,
                                 studyArea,
                                 fires,
                                 currentTime,
                                 firstTime,
                                 growthStepGenerating,
                                 growthStepEnlargingPolys,
                                 growthStepEnlargingLines,
                                 currentDisturbanceLayer,
                                 connectingBlockSize,
                                 disturbanceRateRelatesToBufferedArea,
                                 outputsFolder,
                                 seismicLineGrids,
                                 checkDisturbancesForBuffer = FALSE,
                                 runName,
                                 useRoadsPackage,
                                 siteSelectionAsDistributing,
                                 probabilityDisturbance,
                                 maskWaterAndMountainsFromLines,
                                 featuresToAvoid,
                                 altitudeCut,
                                 clusterDistance,
                                 distanceNewLinesFactor,
                                 runClusteringInParallel,
                                 useClusterMethod,
                                 refinedStructure){
  
  # Extracting layers from previous ones
  # Total study area
  rasterToMatchR <- raster::raster(rasterToMatch)
  probabilityDisturbance <- if (is.null(probabilityDisturbance)) list() else probabilityDisturbance
  studyAreaHash <- digest(studyArea)
  
  if (!is(studyArea, "SpatVector"))
    studyArea <- terra::vect(studyArea)
  
  studyArea <- terra::project(x = studyArea, y = terra::crs(rasterToMatch))
  uniStudyArea <- terra::aggregate(studyArea)
  totalstudyAreaVAreaSqKm <- terra::expanse(uniStudyArea, unit = "km", transform = FALSE)
  totalstudyAreaVAreaSqm <- terra::expanse(uniStudyArea, unit = "m", transform = FALSE)
  totNPix <- sum(rasterToMatch[], na.rm = TRUE)
  calculatedPixelSizem2 <- totalstudyAreaVAreaSqm/totNPix # in m2
  
  # FIRST: Enlarging
  #############################################
  message(crayon::white("Enlarging disturbances..."))
  dPar <- disturbanceParameters[disturbanceType  %in% "Enlarging", ]
  if (nrow(dPar) == 0){
    # Means we don't have any Enlarging (i.e., there is no settlement there!)
    message(crayon::white("No Enlarging disturbances in the study area, returning NULL..."))
    Enlarged <- NULL
  } else {
    whichSector <- dPar[["dataName"]]
    Enlarged <- lapply(1:nrow(dPar), function(ROW) {
      Sector <- dPar[ROW, dataName]
      whichOrigin <- dPar[ROW, disturbanceOrigin]
      updatedL <- lapply(whichOrigin, function(ORIGIN) {
        Lay <- disturbanceList[[Sector]][[ORIGIN]]
        if (is.null(Lay)) {
          stop(
            paste0(
              "Layer of dataName = ",
              Sector,
              " and disturbanceOrigin = ",
              ORIGIN,
              " needs to exist in disturbanceList"
            )
          )
        } # Bug catch for mismatches between disturbanceOrigin and the layers
        # available
        dParOri <- dPar[dataName == Sector & disturbanceOrigin == ORIGIN,]
        # Select the growth rate. For Enlarging we simply buffer the current layer at the percent
        # mentioned in disturbanceRate.
        # Detect current size: terra::expanse(Lay)
        # Calculate the total area to expand it: (currentSize*dParOri[["disturbanceRate"]])
        # NOTE: The value from the table should be in percent because we divide it by 100 to get numeric
        # We need to assume a form for polygons and for lines.
        # Polygons: CIRCLE
        # Lines: RECTANGLES
        # Then we need the formula to calculate the buffer based on the increase in
        # area size (parameter from table) using a normal
        originalForm <- geomtype(Lay)
        if (!disturbanceRateRelatesToBufferedArea){
          if (originalForm %in% c("lines", "points")){
            # Check if we have the resolution
            if ("resolutionVector" %in% dPar){
              if (originalForm == "lines")
                RES <- dPar[index, resolutionVector] else 
                  RES <- dPar[index, resolutionVector]/2
            } else RES <- NULL
            if (is.null(RES)){
              message(crayon::red(paste0("resolutionVector in disturbanceParameters",
                                         " table was not supplied for a lines file vector. Will ",
                                         "default to 15m total (7.5m in each direction).",
                                         " If this is not correct, please provide the ",
                                         "resolution used for developing the layer ", ORIGIN,
                                         " (", Sector,")")))
              if (originalForm == "lines")
                RES <- 15 else 
                  RES <- 7.5
                
            }
            # NOTE: The buffer value is applied on all sides of the shapefile / line
            # This means that a 15m buffer ends up being a 30 m increase in the line
            # in all directions. If the original resolution was 15, then we need to divide it by 2,
            # so we get in the end a TOTAL increase in of 15 m.
            Lay <- terra::buffer(x = Lay, width = RES)
          }
          currArea <- terra::expanse(Lay, unit = "m", transform = FALSE)
        } else { # If the disturbanceRate Relates To Buffered Area (i.e., caribou-wise)
          # It doesn't matter if poly or lines or points, if buffered, is buffered to 500m
          LayBuff <- terra::buffer(Lay, width = 500) # Need to aggregate to avoid double counting!
          LayBuff <- terra::aggregate(x = LayBuff, dissolve = TRUE)
          currArea <- terra::expanse(LayBuff, unit = "m", transform = FALSE) # NEEDS TO BE METERS. USED BELOW
          if (disturbanceRateRelatesToBufferedArea) 
            message(paste0("Buffered (500m) area for ", Sector, 
                           " -- ", ORIGIN, ": ", round(currArea/1000000, 2), " Km2"))
          message(paste0("Percentage of current area: ", round(100*((currArea/1000000)/totalstudyAreaVAreaSqKm), 3), "%."))
        }
        # The dParOri[["disturbanceRate"]] (or disturbRate) is (by default) the calculated difference between
        # the 2010 and 2015 disturbance divided by the total amount of years, over the entire area
        # Therefore, it refers to the total study area size. 
        disturbRate <- dParOri[["disturbanceRate"]]
        #Convert Rate into numeric (i.e. 0.2% = 0.002)
        Rate <- disturbRate/100
        # Disturbance area to add to current:
        if (Rate == 0){
          message(paste0("Rate of disturbance for ", ORIGIN, " is 0. This is either a disturbance ",
                         "that did not (or should not) change, or a disturbance that was reduced through ",
                         "time in the study area. Disturbance reduction has not yet been implemented.",
                         " Returning the layer unchanged."))
          return(Lay)
        }
        disturbAreaSqKm <- Rate*totalstudyAreaVAreaSqKm*dParOri[["disturbanceInterval"]]
        # Total increase in area to be distributed across all disturbances:
        disturbAreaSqM <- disturbAreaSqKm * 1000000 # expected EXTRA disturbed area in sqm 
        # (i.e., how much 0.2% of disturbance over the entire area actually represents)
        # Now calculate how much I actually expect the area to be disturbed (i.e., total future 
        # disturbed area = current disturbed area + extra disturbed area)
        expectedTotalDisturbedArea <- currArea + disturbAreaSqM
        # We will do the buffer iteratively as there are no great solutions for getting the correct 
        # value by calculations because the polygon's forms is not known.
        # While totalDisturbedAreaAchieved is below expectedTotalDisturbedArea, we will keep increasing the 
        # buffer value and recalculating it. Past trials are first below:
        
        totalDisturbedAreaAchieved <- 0
        expandTo <- 0
        iter <- 1
        tictoc::tic("Total elapsed time for calculating disturbance percentage: ")
        while (totalDisturbedAreaAchieved < expectedTotalDisturbedArea) {
          if (iter %in% 1:10){ 
            message(paste0("Calculating total enlarging disturbance size for ", Sector, " for ",  
                           ORIGIN, " (Year ", currentTime,"; iteration ", iter,")"))
          }
          if (iter %% 50 == 0) {
            message(paste0("Calculating total enlarging disturbance size for ", Sector, " for ",  
                           ORIGIN, " (Year ", currentTime,"; iteration ", iter,")"))
          }
          if (originalForm %in% c("lines", "points")){
            expandTo <- expandTo + growthStepEnlargingLines
          } else {
            expandTo <- expandTo + growthStepEnlargingPolys 
          }
          # Expand the area
          LayUpdated <- terra::buffer(x = Lay, width = expandTo)
          # Calculate new total disturbed area
          if (disturbanceRateRelatesToBufferedArea){
            LayUpdatedBuff <- terra::buffer(x = LayUpdated, width = 500)
            futureAreaAgg <- terra::aggregate(LayUpdatedBuff)
            futureArea <- terra::expanse(futureAreaAgg, unit = "m", transform = FALSE)
          } else {
            futureArea <- terra::expanse(LayUpdated, unit = "m", transform = FALSE)
          }
          totalDisturbedAreaAchieved <- sum(futureArea)
          iter <- iter + 1
        }
        tictoc::toc()
        
        message(paste0("Percentage of disturbed future area after buffer: ", 
                       round(100*(totalDisturbedAreaAchieved/totalstudyAreaVAreaSqm), 3), "%."))
        
        cat(crayon::yellow(paste0("Difference between expected and achieved change for ",
                                  crayon::red(Sector), " -- ", crayon::red(ORIGIN), ": ",
                                  crayon::red(format(100*round((totalDisturbedAreaAchieved - 
                                                                  expectedTotalDisturbedArea)/expectedTotalDisturbedArea, 4), 
                                                     scientific = FALSE), " % (ideal value = 0)."),
                                  "\nDisturbance achieved: ", round(totalDisturbedAreaAchieved/1000000,3),
                                  " km2 -- Disturbance expected: ", round(expectedTotalDisturbedArea/1000000,3), " km2")))
        cat(paste0(Sector, 
                   " ", 
                   ORIGIN, 
                   " ",
                   currentTime,
                   " ",
                   format(100*round((totalDisturbedAreaAchieved - expectedTotalDisturbedArea)/expectedTotalDisturbedArea, 8), 
                          scientific = FALSE)),
            file = file.path(Paths[["outputPath"]], paste0("PercentageDisturbances_", currentTime, 
                                                           "_", runName, ".txt")),
            append = TRUE, sep = "\n")
        
        return(LayUpdated)
      })
      names(updatedL) <- whichOrigin
      return(updatedL)
    })
    names(Enlarged) <- whichSector
  }
  #############################################
  
  
  # SECOND: Generating
  #############################################
  message(crayon::white("Generating disturbances..."))
    # First, mask the current disturbances on the potential layer so we don't choose them again
  # Get the resolution for lines
  RES <- unique(disturbanceParameters[["resolutionVector"]])
  if (length(RES) > 1)
    stop(paste0("Different resolutions have not yet been implemented.",
                "Please modify the code to allow for it"))
  if (is.null(currentDisturbanceLayer)){
    message(paste0("currentDisturbanceLayer is NULL. This is likely the first year of the",
                   " simulation. Creating current disturbed polygons..."))
    # Current layer is null, which indicates this might be the first year of the simulation. I
    # will need to create it then, based on the existing disturbances.
    # Get all current disturbances
    unDL <- unlist(disturbanceList)
    allCurrentDistLayNames <- names(unDL)[!grepl(x = names(unDL), pattern = "potential")]
    currDist <- unDL[names(unDL) %in% allCurrentDistLayNames]
    # Combine all layers, both polygons and lines, separtely
    linesLays <- currDist[sapply(currDist, geomtype) == "lines"]
    polyLays <- currDist[sapply(currDist, geomtype) == "polygons"]
    linesLays <- unlist(lapply(linesLays, function(X){
      terra::buffer(x = X, width = RES)
      # Here we need to buffer as this layer will be used to avoid specific locations (i.e., already 
      # disturbed!) to be selected for new disturbances.
    }))
    linesAndPolys <- do.call(what = c, list(polyLays, linesLays))
    names(linesAndPolys) <- NULL
    allLays <- do.call(rbind, linesAndPolys)
  } else {
    if (is(currentDisturbanceLayer, "RasterLayer"))
      allLaysRas <- terra::rast(currentDisturbanceLayer) else
        if (is(currentDisturbanceLayer, "SpatRaster")) 
          allLaysRas <- currentDisturbanceLayer else {
            curDistVcs <- unlist(currentDisturbanceLayer)
            curDistVcsAll <- lapply(names(curDistVcs), function(eachVectNm){
              eachVect <- curDistVcs[[eachVectNm]]
              if (length(eachVect) == 0) return(NULL)
              message(paste0("Buffering and/or merging polygons for ", eachVectNm))
              if (geomtype(eachVect) != "polygons"){
                if (geomtype(eachVect) == "points"){
                  buffVect <- terra::buffer(eachVect, width = RES/2)
                } else {
                  buffVect <- terra::buffer(eachVect, width = RES)
                }
              } else {
                buffVect <- terra::aggregate(eachVect, dissolve = TRUE)
              }
              return(buffVect)
            })
            names(curDistVcsAll) <- NULL
            allLays <- do.call(rbind, curDistVcsAll)
          }
        # Create allLaysRas, which is the disturbance raster with all disturbances rasterized, with
        # value == 1 for the disturbed pixels. Non-disturbed pixels are NA
  }

  dPar <- disturbanceParameters[disturbanceType  %in% "Generating", ]
  if (nrow(dPar) == 0){
    # Means we don't have any Generated (i.e., there is nothing disturbed there!)
    message(crayon::white("No Generated disturbances in the study area, returning NULL..."))
    Generated <- NULL
  } else {
    whichSector <- dPar[["dataName"]]
    Generated <- lapply(1:nrow(dPar), function(ROW) {
      Sector <- dPar[ROW, dataName]
      whichOrigin <- dPar[ROW, disturbanceOrigin]
      updatedL <- lapply(whichOrigin, function(ORIGIN) {
        dParOri <- dPar[dataName == Sector & disturbanceOrigin == ORIGIN,]
        Lay <- disturbanceList[[Sector]][[ORIGIN]] # If there are no disturbances, Lay will be NULL
        if (is.null(Lay)) {
          message(crayon::red(paste0("Layer of dataName = ", Sector, " and disturbanceOrigin = ",
                                     ORIGIN, " doesn't exist in disturbanceList. Potentially this disturbance ",
                                     "hasn't happened in the area yet. The function will try to created it based ",
                                     "on the information from the disturbanceParameters and potential layer.")))
        }
        # Select the growth rate. For Generating we get the original layer (Potential) and generate 
        # X polygons based on total size (represented by the rate in function of currently 
        # existing structures) / average size of polygons (using the function to randomly select the 
        # number and sizes)
        
        # 1. Get the potential layer
        potLay <- disturbanceList[[Sector]][[dParOri[["dataClass"]]]]
        if (is.null(potLay)){
          # No potential layer, but disturbances happened.
          message(crayon::white(paste0("No potential for disturbances Sector ", Sector,
                                       " in the study area, returning NULL...")))
          newDisturbs <- NULL
        } else {
          # Make sure we have only available places (i.e., remove fire places for forestry and 
          # remove the existing disturbances for all!)
          
          if (Sector == "forestry"){
            # Update with fire layer if forestry
            
            message(paste0("Generating disturbance for forestry. Updating potential layer for ",
                           "occurred fires and currently productive forest for year ", currentTime))
            
            # First: Select only productive forests
            # Previous pixel strategy for Generating Disturbances
            # 2. Fasterize it
            potLaySF <- sf::st_as_sf(x = potLay)
            potField <- dParOri[["potentialField"]]
            
            if (any(is.na(potField), 
                    potField == "",
                    is.null(potField))){ 
              # If NA, it doesn't matter, but need to create a 
              # field so not the whole thing becomes one big 1 map
              potField <- "Potential"
              potLaySF$Potential <- 1
            }
            potLaySF <- subset(potLaySF, potLaySF$ORIGIN < (currentTime - 50))
            potLayF <- fasterize::fasterize(sf = st_collection_extract(potLaySF, "POLYGON"),
                                            raster = rasterToMatchR, field = potField)
            # Second: remove fires
            if (!is.null(fires))
              potLayF[fires[] == 1] <- NA
            # Convert the fire raster to polygons
            potLayT <- rast(potLayF)
            potLay <- terra::as.polygons(potLayT)
            # Third: give preference to cutblocks that are closer to current cutblocks
            # Using terra is quicker!
            # # WITH NWT DATA, THE FOLLOWING LINES ARE USELESS AS ALL FOREST IS CLOSE ENOUGH
            # # WE CAN IMPROVE THAT FOR OTHER AREAS, HOWEVER, BY USING BUFFERING METHODS. 
            # Something in these lines but for polygons...
            # potLayFt <- terra::rast(potLayF)
            # tictoc::tic("Distance raster time elapsed: ")
            # distRas <- terra::distance(potLayFt) # Distance raster time elapsed: : 8348.948 sec elapsed
            # toc()
            # distRas2 <- (distRas-maxValue(distRas))*-1
            # distRas2[is.na(potLayF)] <- NA
          }
          
          
          # We need to aggregate the potential layer to make sure we have all possible polygons with the 
          # same values with the same probability of developing the disturbance
          aggby <- "Potential"
          if (!"Potential" %in% names(potLay)){ # If we don't have Potential, we need to create it
            if (length(names(potLay)) == 1){
              names(potLay) <- "Potential"
            } else {
              potLay[["Potential"]] <- 1
            }
          }
          if (NROW(potLay) > 1)
            potLay <- aggregate(potLay, by = aggby, dissolve = TRUE, count = FALSE)
          if (length(potLay) == 0){# In case the cropped area doens't have anything
            message(paste0("The potential area for ", Sector, " class ", ORIGIN, " is NULL.",
                           " Likely cropped out from studyArea. Returning NULL."))
            
            return(NULL)
          }
          
          # Get the disturbance rate
          disturbRate <- dParOri[["disturbanceRate"]]
          #Convert Rate into numeric (i.e. 0.2% = 0.002)
          Rate <- disturbRate/100
          # Here we multiply the rate by the total number of pixels by the interval we are using to 
          # generate the disturbances (as rate is passed as yearly rate!), to know how many pixels we 
          # are expected to disturb in the given interval for this study area, given the rate applied
          if (Rate == 0){
            message(paste0("Rate of disturbance for ", ORIGIN, " is 0. This is either a disturbance ",
                           "that did not (or should not) change, or a disturbance that was reduced through ",
                           "time in the study area. Disturbance reduction has not yet been implemented.",
                           " Returning the layer unchanged."))
            
            newDistLay <- subset(potLay, subset = NA) # Create a template to rbind with the new disturbances
            return(newDistLay)
          }
          expectedNewDisturbAreaSqKm <- Rate*totalstudyAreaVAreaSqKm*dParOri[["disturbanceInterval"]]
          # Total increase in area to be distributed across all disturbances:
          expectedNewDisturbAreaSqM <- expectedNewDisturbAreaSqKm * 1000000 # expected EXTRA disturbed area in sqm 
          # (i.e., how much 0.2% of disturbance over the entire area actually represents)
          
          # Before I do the iterations I wanna know what is the currently disturbed area for this specific disturbance type
          if (all(!is.null(Lay),
                  disturbanceRateRelatesToBufferedArea)){
            if (is(Lay, "RasterLayer"))
              Lay <- rast(Lay)
            if (is(Lay, "SpatRaster")){ # Need to convert to polygon for area
              Lay[Lay == 0] <- NA # Otherwise buffers weird places!
              Lay <- terra::as.polygons(Lay)
            }
            LayBuff <- terra::buffer(Lay, width = 500) # Need to aggregate to avoid double counting!
            LayBuff <- terra::aggregate(x = LayBuff, dissolve = TRUE)
            currArea <- terra::expanse(LayBuff, unit = "m", transform = FALSE) # NEEDS TO BE METERS. USED BELOW
            if (length(currArea) == 0){
              message(paste0("No disturbance for ", Sector, 
                             " -- ", ORIGIN, ": currentArea = 0 Km2"))
            } else {
              message(paste0("Buffered (500m) area for ", Sector, 
                             " -- ", ORIGIN, ": ", round(currArea/1000000, 2), " Km2"))
              message(paste0("Percentage of current area: ", round(100*((currArea/1000000)/totalstudyAreaVAreaSqKm), 3), "%."))
            }
          }
          # Make sure we select the areas to disturb in a way we have enough space to disturb the necessary
          # sizes
          ITR <- 1
          totalAreaAvailable <- 0
          # # Now check if there is enough space for the expected disturbance!
          # # 1.2. We need to chose the best places:
          while (totalAreaAvailable < expectedNewDisturbAreaSqM){
            # Get the possible places for the disturbances. Make sure there is enough place! Best potential 
            # are the places with max values
            valuesAvailable <- sort(as.numeric(unlist(potLay[[names(potLay)]])))
            if (ORIGIN %in% siteSelectionAsDistributing){
              rowsToChoose <- 1:length(valuesAvailable)
            } else {
              rowsToChoose <- (length(valuesAvailable)-(ITR-1)):length(valuesAvailable)
            }
            potLayTop <- potLay[rowsToChoose]
            # Now exclude where disturbance already exists
            # Use intersect to see if they overlap.
            if (!ORIGIN == "seismicLines") # For seismic lines, intersection happens below
              message(paste0("Intersecting existing disturbances with potential for development for ",
                             Sector, " -- ", ORIGIN))
            if (Sector == "forestry"){
              croppedAllLays <- reproducible::cropInputs(allLays, potLayTop)
              potLayTopValid <- terra::erase(potLayTop, croppedAllLays) # BEST AREA!
            } else
              if (!ORIGIN == "seismicLines") { # If "seismicLines", we don't need to erase as they can overlap
                doIntersect <- terra::intersect(allLays, potLayTop)
                if (nrow(doIntersect) > 0){
                  # If so, erase.
                  croppedAllLays <- reproducible::cropInputs(allLays, potLayTop)
                  potLayTopValid <- terra::erase(potLayTop, croppedAllLays) # BEST AREA!
                } else {        
                  # If not, move on
                  potLayTopValid <- potLayTop
                }
              } else potLayTopValid <- potLayTop
            
            # Now check if there is enough space for the expected disturbance!
            totalAreaAvailable <- tryCatch(terra::expanse(terra::aggregate(potLayTopValid, dissolve = TRUE), 
                                                          transform = FALSE, unit = "m"), error = function(e) browser())
            ITR <- ITR + 1
          }
          
          # If probabilityDisturbance is NOT provided, calculate
          if (all(ORIGIN %in% siteSelectionAsDistributing,
                  is.null(probabilityDisturbance[[ORIGIN]]))){ # NOTE: Potentiall yery time consuming!
            message(paste0("probabilityDisturbance for ", ORIGIN, " is NULL. Calculating from data..."))
            # 1. Extract the total area of each polygon type
            potLayTopValid$areaPerPoly <- terra::expanse(potLayTopValid, unit = "m", transform = FALSE)
            # 2. Sum all to calculate the total area of all polys
            totAreaPolys <- sum(potLayTopValid$areaPerPoly)
            # 3. Get the total disturbance area for each polygon type       
            areaDistPerPoly <- terra::intersect(Lay, potLayTopValid) # This is super time demanding. 
            # Allow to pass a proportion. This would also avoid completely selecting all polygons if passed! 
            # 4. Calculate the total percentage of the disturbance per polygon type
            buffSeisL <- terra::buffer(areaDistPerPoly, width = 3) # Need to aggregate to avoid double counting!
            areaDistPerPoly$area <- terra::expanse(buffSeisL, unit = "m", transform = FALSE)
            areaDT <- as.data.table(areaDistPerPoly[, c("Potential", "area")])
            totAreaDT <- areaDT[, disturbedArea := sum(area), by = Potential]
            totAreaDT[, area := NULL]
            totAreaDT <- unique(totAreaDT)
            totAreaDT[, probPoly := disturbedArea/sum(totAreaDT$disturbedArea)]
            # 5. This percentage is = to the probability a disturbance will occur in that polygon --> Need to pass this value!
            probabilityDisturbance[[ORIGIN]] <- totAreaDT[, c("Potential", "probPoly")]
          } else {
            # If provided, test that probabilityDisturbance matches the Potential
            potValsPassed <- unique(sort(probabilityDisturbance[[ORIGIN]][["Potential"]]))
            potValsLay <- unique(sort(potLayTopValid$Potential))
            passTest <- all(potValsLay %in% potValsPassed)
            if (!passTest){
              if (!is.null(Lay)){
                # Cases where we just don't have the probability for the area will 
                # have Lay as NULL.
                # We will not have potential values in this case because we don't 
                # have the data.These should be fine. 
                
                # We could calculate the probability for the area, but it is not happening for some wicked reason.
                warning(paste0("The probabilityDisturbance was provided, but do not match the expected Potential values for ",
                               ORIGIN, ".\n",
                               "Potential values passed: ", paste(potValsPassed, collapse = ", "), "\n",
                               "Potential values in layer: ", paste(potValsLay, collapse = ", "), 
                               ". This may happen in small study areas and should not be",
                               "cause of concern but may help identify errors."), 
                        immediate. = TRUE)
              }
            }
          }
          # For seismic Lines we need to do some area processing of the before we can 
          # generate the disturnances. And we don't want to repeat this every time inside
          # a while loop, as it doesn't change, so we do it outside
          if (ORIGIN == "seismicLines"){
            if (firstTime){
              message("First time generating sesmicLines, proceeding with clustering...")
              cropLayFinal <- Cache(createCropLayFinalYear1, 
                                    Lay = Lay, 
                                    potLayTopValid = potLayTopValid, 
                                    runClusteringInParallel = runClusteringInParallel, 
                                    clusterDistance = clusterDistance, 
                                    studyAreaHash = studyAreaHash)
              terra::writeVector(x = cropLayFinal, 
                                 filename = file.path(outputsFolder, 
                                                      paste0("seismicLinesYear", 
                                                             currentTime, 
                                                             "_",
                                                             studyAreaHash,
                                                             ".shp")),
                                 overwrite = TRUE)
            } else {
              cropLayFinal <- Lay
            }
          }
          
          # 1. Make the iteration, and while the area is not achieved, continue 
          areaChosenTotal <- 0
          IT <- 1
          newDisturbs <- NULL
          alreadyReduced <- FALSE
          
          while (areaChosenTotal < expectedNewDisturbAreaSqM){
            if (IT %in% 1:10){
              message(paste0("Calculating total generated disturbance size for ", Sector, " for ", 
                             ORIGIN, " (Year ", currentTime,"; iteration ", IT, ", ", 
                             round(100*(round(areaChosenTotal, 0)/round(expectedNewDisturbAreaSqM, 0)), 2),"% achieved)"))
            }
            if (IT %% 10 == 0){
              message(paste0("Calculating total generated disturbance size for ", Sector, " for ",  
                             ORIGIN, " (Year ", currentTime,"; iteration ", IT, ", ", 
                             round(100*(round(areaChosenTotal, 0)/round(expectedNewDisturbAreaSqM, 0)), 2),"% achieved)"))
            }
            # 1.1. Get each disturbance's size. If the expected disturbance is smaller than the 
            # normal size, we only disturb what is expected
            Size <- round(eval(parse(text = dParOri[["disturbanceSize"]])), 0)
            # If the Size is larger than expectedNewDisturbAreaSqM, we might first make a probability of 
            # the disturbance happening. 
            if (Size > expectedNewDisturbAreaSqM){
              p <- rbinom(1, size = 1, prob = (expectedNewDisturbAreaSqM/Size))
              disturbanceHappening <- if (p == 1) TRUE else FALSE
            } else disturbanceHappening <- TRUE
            # If the disturbance doesn't happen:
            if (!disturbanceHappening){
              message(paste0("Rate of disturbance for ", ORIGIN, " is very small and the probability of ",
                             "this disturbance happening returned 0.",
                             " Returning layer without disturbances."))
              break
            }
            # 1.3. Once the best places are chosen, we place the disturbance in a new layer, which will 
            # have to be updated until we leave the while loop
            # Make an inner buffer with the half of the size of the Size to make sure we have the full 
            # new disturbance within the area designated! 
            # UPDATE: I can't do this as it is crashing RStudio. I will just hope that all fall within 
            # the area. Probably crashing because the area is smaller than I am asking it to make the inner buffer
            # potLayTopValidIn <- terra::buffer(potLayTopValid, width = -(Size/2))
            # potLayTopValidIn <- terra::makeValid(potLayTopValidIn)
            if (ORIGIN == "seismicLines"){
              if (useClusterMethod){
                
                # NOTES: Seismic lines, because we don't connect these, can reach close the specified information
                # on total amount of disturbance. However, the lines created are RANDOM and NOT exactly as the 
                # existing ones, so we still need to do it iteratively.
                
                # Here I need to think about a way to assess how much each cluster represents in the total needed area
                # expectedNewDisturbAreaSqM so I can randomly choose, with higher probability in better areas 
                # clusters to be duplicated. This duplication will then get VERY close to the expected area,
                # likely slightly below (as not many lines overlap) --> The smaller the buffer 
                # (clusterDistance and distanceNewLinesFactor) new lines will be closer 
                
                if (geomtype(cropLayFinal) != "lines"){
                  print("cropLayFinal needs to be lines and needs to have all previous lines")
                  browser()
                }
                
                # 0. Calculate the area of each line buffered by 3m (min width in the field)
                cropLayFinal$buff3mAreaM2 <- expanse(buffer(x = cropLayFinal, width = 3))
                
                # 1. Add a column in how much each cluster represents in the total in m2 -- not by individual line!
                cropLayFinalDT <- as.data.table(as.data.frame(cropLayFinal))
                if (!"Pot_Clus" %in% names(cropLayFinalDT)){
                  message("Pot_Clus not found in cropLayFinalDT. Debug")
                  browser()
                }
                cropLayFinalDT[, sumBuff3mAreaM2 := sum(buff3mAreaM2), by = "Pot_Clus"] 
                # Need to do by potential as for each potential, cluster numbers are repeated
                totalBuff3mArea <- sum(cropLayFinalDT$buff3mArea)
                cropLayFinalDT[, PercBuff3mAreaOfTotalM2 := 100*(sumBuff3mAreaM2/totalBuff3mArea),  by = "Pot_Clus"]
                if (sum(unique(cropLayFinalDT[, c("Pot_Clus","PercBuff3mAreaOfTotalM2")]$PercBuff3mAreaOfTotalM2)) > 100.001){ 
                  message(paste0("Total contribution of clusters in total area is higher than 100%.",
                                 "Something may be wrong. Entering debug mode."))
                  browser()
                } # TODO test
                # 2. Choose randomly clusters (with different probabilities) that sum to the total expected new. 
                # sumBuff3mAreaM2 --> by cluster
                # PercBuff3mAreaOfTotalM2 --> representation if each cluster over the total
                # Potential --> Represents the highest potential for being chosen.
                # --> Draw for all potentials, the probability a cluster within these will be chosen (higher potential, higher chances)
                # Normalize probabilities to sum to 1
                probabilities <- unique(cropLayFinalDT$Potential) / sum(unique(cropLayFinalDT$Potential))
                if (length(unique(cropLayFinalDT$Potential)) == 1){
                  sampledClusters <- rep(unique(cropLayFinalDT$Potential), times = growthStepEnlargingLines)### <~~~~~~~~~~~~ Changed here from 1000 to 1 to try increasing iterations number in Seismic lines, currently overdoing it.
                } else {
                  sampledClusters <- sample(unique(cropLayFinalDT$Potential), 
                                            size = growthStepEnlargingLines, ### <~~~~~~~~~~~~ Changed here from 1000 to 1 to try increasing iterations number in Seismic lines, currently overdoing it. 
                                            replace = TRUE, 
                                            prob = probabilities)
                }
                selectedClusters <- NULL
                for (uniqueSampClus in unique(sampledClusters)){
                  toChoseFrom <- unique(cropLayFinalDT[Potential == uniqueSampClus, Pot_Clus])
                  howManyINeed <- sum(sampledClusters == uniqueSampClus)
                  if (length(toChoseFrom) == 1){
                    sampledOnes <- rep(toChoseFrom, times = howManyINeed)
                  } else {
                    sampledOnes <- sample(toChoseFrom,
                                          replace = TRUE,
                                          size = howManyINeed)
                  }
                  selectedClusters <- c(selectedClusters,sampledOnes)
                }              
                #TODO Here I can parallelize using future!
                newLines <- simulateLines(Lines = cropLayFinal[cropLayFinal$Pot_Clus %in% selectedClusters, ],
                                          distThreshold = clusterDistance,
                                          distNewLinesFact = distanceNewLinesFactor,
                                          refinedStructure = refinedStructure)
                # newDist NEEDS to be buffered, first by 3m and then (happening below in the code) by 
                # 500m, if 
                # disturbanceRateRelatesToBufferedArea 
                newDist <- terra::buffer(newLines, width = 3) # THIS IS WHAT I NEED IN THE END HERE!
              } else {
                distApropExp <- areaChosenTotal/expectedNewDisturbAreaSqM
                if (all(!alreadyReduced, 
                        distApropExp > 0.85)){
                  seismicLineGridsO <- seismicLineGrids
                  seismicLineGrids <- round(((seismicLineGrids*IT*(1-distApropExp))/distApropExp)/4, 0)
                  message(paste0("Total expected disturbance almost achieved, reducing the number of grids from ", 
                                 seismicLineGridsO, " to ", seismicLineGrids))
                  alreadyReduced <- TRUE
                }
                rs <- matrix(runif(2*seismicLineGrids, 50, 100), ncol = 2, byrow = TRUE) # For the raster below
                if (ORIGIN %in% siteSelectionAsDistributing){
                  # 6. Select from the polygon finalPotLay the one which will have the disturbance by applying the probability
                  polyToChoose <- sample(probabilityDisturbance[[ORIGIN]][["Potential"]],
                                         size = 1,
                                         prob = probabilityDisturbance[[ORIGIN]][["probPoly"]])
                  currFinalPotLay <- finalPotLay[finalPotLay$Potential == polyToChoose]
                  centerPoint <- terra::spatSample(currFinalPotLay, size = seismicLineGrids, method = "random")
                } else {
                  centerPoint <- terra::spatSample(finalPotLay, size = seismicLineGrids, method = "random")              
                }
                while (length(centerPoint) != seismicLineGrids){
                  # If center point has less rows than size, it means that the sampling likely got all 
                  # the area possible if replace = FALSE, the default. So the missing ones need to be 
                  # added to the center point in the next best available area.
                  howManyMissing <- seismicLineGrids-length(centerPoint)
                  if (howManyMissing < 0){
                    message(paste0("spatSample has sampled ", length(centerPoint), " while the expected ",
                                   "number of points was ", seismicLineGrids, ". Entering browser mode."))
                    browser()
                  } 
                  message(paste0("spatSample has sampled ", length(centerPoint), " while the expected ",
                                 "number of points was ", seismicLineGrids, ". Choosing more points from ",
                                 "next best area..."))
                  valsToExclude <- as.numeric(unique(potLayTopValid[["Potential"]]))
                  nextBestValue <- max(setdiff(valuesAvailable, valsToExclude))
                  rowsToChoose <- which(valuesAvailable == nextBestValue)
                  potLayTopValid <- potLay[rowsToChoose]
                  # 1. UPDATING THE LAYER: 
                  cropLay <- postProcessTo(Lay, potLayTopValid)
                  cropLayBuf <- buffer(cropLay, width = 50)
                  cropLayAg <- aggregate(cropLayBuf, dissolve = TRUE)
                  finalPotLay <- erase(potLayTopValid, cropLayAg)
                  centerPointToAdd <- terra::spatSample(finalPotLay, size = howManyMissing, method = "random")
                  centerPoint <- rbind(centerPoint, centerPointToAdd)
                }
                # 5. Draw the new grid based on the total length expected
                # 5.1. Draw a square based on the centerPoint, where the distance from point to the lines 
                # is the diagonal of a square of the lineLenght you want.
                
                print("Currently not using, but should test! Something is likely not working")
                # browser() # HERE is where Mean and SD is used from cropLay
                lineLength <- rtnorm(length(centerPoint), Mean, Sd, lower = 0)
                # 5.2. Make a square polygon with the center point and the distance
                gridReady <- lapply(1:seismicLineGrids, function(ROW){
                  pnt <- vect(matrix(c(xmin(centerPoint[ROW,])-(lineLength[ROW]/2), ymin(centerPoint[ROW,]),
                                       xmin(centerPoint[ROW,])+(lineLength[ROW]/2), ymin(centerPoint[ROW,]),
                                       xmin(centerPoint[ROW,]), ymin(centerPoint[ROW,])+(lineLength[ROW]/2),
                                       xmin(centerPoint[ROW,]), ymin(centerPoint[ROW,])-(lineLength[ROW]/2)), 
                                     ncol = 2, byrow = TRUE), type="points", crs=crs(centerPoint[ROW,]))
                  polyArea <- as.polygons(pnt, extent=TRUE)
                  # 5.3. Create a raster with the desired distance between the points
                  polRas <- tryCatch(rast(polyArea, resolution = rs[ROW,]), error = function(e) browser())
                  # To points: This is now a vector
                  gridPoints <- suppressWarnings(as.points(polRas, na.rm = TRUE))
                  # Now I need to find the first and last dots to connect
                  # Rows: 
                  # NOTE: We can tilt by simply choosing different end points!!
                  rowsOfPointsL <- lapply(0:(dim(polRas)[2]-1), function(rowIndex){
                    thePair <- c(1, (ncell(polRas)-(dim(polRas)[2]-1)))+rowIndex
                    # create the line based on the points by extracting the points based on thePair 
                    theLine <- as.lines(gridPoints[thePair, ])
                    return(theLine)
                  })
                  rowsOfPoints <- do.call(rbind, rowsOfPointsL)
                  # Cols:
                  colsOfPointsL <- lapply(0:(dim(polRas)[1]-1), function(colIndex){
                    thePair <- c(1, dim(polRas)[2])+(colIndex*dim(polRas)[2])
                    # create the line based on the points by extracting the points based on thePair 
                    theLine <- as.lines(gridPoints[thePair, ])
                    return(theLine)
                  })
                  colsOfPoints <- do.call(rbind, colsOfPointsL)
                  # Now we do 0 to two random crossing lines
                  howManyCross <- sample(x = seq(0, 2), size = 1)
                  if (howManyCross > 0){
                    # Find all points that are in each edge
                    up <- c(1:dim(polRas)[2])
                    bottom <- c((ncell(polRas)-(dim(polRas)[2]-1)):ncell(polRas))
                    left <- NULL
                    for (i in 0:(dim(polRas)[1]-1)){
                      left <- c(left, 1+(i*dim(polRas)[2]))
                    }
                    right <- left + (dim(polRas)[2]-1)
                    allDirs <- c("up", "bottom", "left", "right")
                    crossLines <- lapply(1:howManyCross, function(tms){
                      dir1 <- sample(x = allDirs, 1)
                      dir2 <- sample(setdiff(allDirs, dir1), 1)
                      thePair <- c(sample(get(dir1), 1), sample(get(dir2), 1))
                      # create the line based on the points by extracting the points based on thePair 
                      theLine <- as.lines(gridPoints[thePair, ])
                      return(theLine)
                    })
                    crossedLines <- do.call(rbind, crossLines)
                  } else {
                    crossedLines <- howManyCross <- NULL
                  }
                  gridReady <- do.call(rbind, list(rowsOfPoints, colsOfPoints, crossedLines))# Is the newly created grid
                  return(gridReady)
                })
                gridReadyB <- do.call(rbind, gridReady)
                # newDist NEEDS to be buffered, first by 3m and then by by 500m, if disturbanceRateRelatesToBufferedArea 
                # Although dParOri[["resolutionVector"]] has vector resolution, seismic lines are much slimmer.
                # In the past, they used to be placed 300–500 m apart, and 5–10 m wide. 
                # Nowadays, the average is 3m wide (2-4, usually not more than 5.5m) and about 50–100m apart 
                # (Dabros et al., 2018 - DOI: 10.1139/er-2017-0080)
                newDist <- terra::buffer(gridReadyB, width = 3)
              }
            } else {
              # Get a point within the layer:
              if (ORIGIN %in% siteSelectionAsDistributing){
                polyToChoose <- sample(probabilityDisturbance[[ORIGIN]][["Potential"]],
                                       size = 1,
                                       prob = probabilityDisturbance[[ORIGIN]][["probPoly"]])
                currFinalPotLay <- potLayTopValid[potLayTopValid$Potential == polyToChoose]
                centerPoint <- terra::spatSample(currFinalPotLay, size = seismicLineGrids, method = "random")             
              } else {
                centerPoint <- terra::spatSample(potLayTopValid, size = 1, method = "random")
              }
              while (nrow(centerPoint) < 1){
                centerPoint <- terra::spatSample(potLayTopValid, size = 1, method = "random")
              }
              # Here we create a random study area based on the center point for the area where 
              # the disturbance is more likely to happen.
              newDist <- SpaDES.tools::randomStudyArea(center = centerPoint, size = Size)
              
              ### CATCHING PROBLEMS ###
              if (is.na(newDist)) browser()
              #########################
              
              names(newDist) <- names(centerPoint)
              # 1.4. If we relate to buffered area, we need to make an inner buffer, as it is included in 
              # areaChosenTotal.
              if (disturbanceRateRelatesToBufferedArea){
                # Then, if disturbanceRateRelatesToBufferedArea, we buffer this point to identify what 
                # if the minimum area we can have. Then, we also calculate the area of the new disturbance 
                # and we check if the area of the new disturbance is smaller than the minimum area of the 
                # buffered point. 
                minDistBuff <- terra::buffer(centerPoint, width = 500)
                areaMinDB <- terra::expanse(minDistBuff)
                areaNewDist <- terra::expanse(newDist)
                # If the new disturbance is smaller than its buffered to 500m version, means we can't reduce it, so
                # we use just a point instead, so when it is buffered, it has the similar size expected.
                if (areaNewDist < areaMinDB){
                  newDist <- centerPoint
                } else {
                  # Otherwise, we buffer innwards and that's our disturbance
                  newDist <- terra::buffer(newDist, width = -500)
                  if (any(length(newDist) == 0, # In this case, the disturbance was almost exactly 
                          # 500m but just a wee bigger, so the result is 0 but it passed the other tests
                          !terra::is.valid(newDist),
                          any(is.na(as.vector(ext(newDist)))))){
                    newDist <- centerPoint
                  }
                }
              }
            }
            
            if (disturbanceRateRelatesToBufferedArea){
              newDistBuff <- terra::buffer(newDist, width = 500) 
              if (nrow(newDistBuff) > 1){
                newDistBuff <- terra::aggregate(newDistBuff)
              }
              areaChosenTotal <- areaChosenTotal + terra::expanse(newDistBuff, unit = "m", transform = FALSE)
            } else { 
              areaChosenTotal <- areaChosenTotal + terra::expanse(newDist, unit = "m", transform = FALSE)
            }
            if (geomtype(newDist) %in% c("points", "lines")){
              if (IT == 2)
                message("Types of geometry differ, buffering new disturbance to a minimum value")
              wid <- 0.0000000003  
              newDistBuf <- terra::buffer(newDist, width = wid)
              while (any(is.na(as.vector(ext(newDistBuf))))){
                message(paste0("Minimum buffering size (", wid,") failed. Increasing buffering size..."))
                wid <- wid + wid
                newDistBuf <- terra::buffer(newDist, width = wid)
              }
              newDist <- newDistBuf
            }
            
            if (ORIGIN == "seismicLines"){
              # Crop again because are being simulated out of SA
              newDist <- reproducible::postProcess(newLines, studyArea)
            }
            
            if (IT == 1){
              newDisturbs <- newDist
            } else {
              newDisturbs <- rbind(newDisturbs, newDist)
            }
            IT <- IT + 1
          }
          
          
          message(paste0("Percentage of disturbed future area after buffer: ", 
                         round(100*(areaChosenTotal/totalstudyAreaVAreaSqm), 4), "%."))
          
          cat(crayon::yellow(paste0("Difference between expected and achieved change for ",
                                    crayon::red(Sector), " -- ", crayon::red(ORIGIN), ": ",
                                    crayon::red(format(100*round((areaChosenTotal - 
                                                                    expectedNewDisturbAreaSqM)/expectedNewDisturbAreaSqM, 4), 
                                                       scientific = FALSE), " % (ideal value = 0)."),
                                    "\nDisturbance achieved: ", round(areaChosenTotal/1000000,3),
                                    " km2 -- Disturbance expected: ", round(expectedNewDisturbAreaSqM/1000000,3), " km2")))
          
          cat(paste0(Sector, 
                     " ", 
                     ORIGIN, 
                     " ",
                     currentTime,
                     " ",
                     format(100*round((areaChosenTotal - expectedNewDisturbAreaSqM)/expectedNewDisturbAreaSqM, 8), 
                            scientific = FALSE)),
              file = file.path(Paths[["outputPath"]], paste0("PercentageDisturbances_", currentTime, 
                                                             "_", runName, ".txt")),
              append = TRUE, sep = "\n")
        }
        return(newDisturbs)
      })
      names(updatedL) <- whichOrigin
      return(updatedL)
    })
    names(Generated) <- whichSector
  }
  #############################################
  
  # THIRD: CONNECTING (Use generated and table)
  #############################################
  message(crayon::white("Connecting disturbances..."))
  dPar <- disturbanceParameters[disturbanceType  %in% "Connecting", ]
  if (nrow(dPar) == 0){
    # Means we don't have any Connecting (i.e., there is nothing new!)
    message(crayon::white("No Connecting disturbances in the study area, returning NULL..."))
    Connected <- NULL
  } else {
    Connected <- lapply(1:NROW(dPar), function(index) {
      print(paste0("Connected Index ", index))
      Sector <- dPar[index, dataName]
      disturbanceOrigin <- dPar[index, disturbanceOrigin]
      disturbanceEnd <- dPar[index, disturbanceEnd]
      # 1. Get the info on which layer to connect to which layer
      oriLay <- Generated[[Sector]][[disturbanceOrigin]]
      if (is.null(oriLay)) {
        # Find alternative dataName for roads
        foundSectors <- unique(dPar[disturbanceOrigin == dPar[index, disturbanceOrigin], dataName])
        newSector <- foundSectors[foundSectors != "roads"]
        if (length(newSector) != 0)
          oriLay <- Generated[[newSector]][[disturbanceOrigin]]
        if (is.null(oriLay)){
          # This might happen for mining for example, as it is the only 
          # one that has connecting from mines to roads, but no other 
          # layer in disturbanceParameters[disturbanceType  %in% "Connecting", ]
          oriLay <- Generated[[disturbanceOrigin]][[disturbanceOrigin]]
          if (is.null(oriLay)){
            if (disturbanceOrigin == "cutblocks") {
              oriLay <- Generated[["forestry"]][[disturbanceOrigin]]
            }
            if (is.null(oriLay)){
              warning(paste0("The disturbance origin layer (", Sector," for ", disturbanceOrigin,
                             ") couldn't be found in the disturbanceList. It is possible that the",
                             " disturbance is not available in the study area. Returning NULL."), 
                      immediate. = TRUE)
              return(NULL)
            }
          } # Bug catch for mismatches between disturbanceOrigin and the layers
        }
      }
      # available
      # 3. Get the end layer (from )
      endLay <- disturbanceList[[Sector]][[disturbanceEnd]]
      # If we are dealing with 'roads' we have a special case where the endLay will be NULL. Then
      # we can try to mitigate by checking if we are expecting the road layer, and check in the 
      # disturbanceList for it
      if (is.null(endLay)) {
        # Find alternative dataName for roads
        foundSectors <- unique(dPar[disturbanceOrigin == dPar[index, disturbanceOrigin], dataName])
        newSector <- foundSectors[foundSectors != "roads"]
        if (length(newSector) != 0) # This might happen for mining for example, as it is the only 
          # one that has connecting from mines to roads, but no other 
          # layer in disturbanceParameters[disturbanceType  %in% "Connecting", ]
          endLay <- Generated[[newSector]][[disturbanceOrigin]]
        if (is.null(endLay)){
          endLay <- Generated[[disturbanceOrigin]][[disturbanceOrigin]]
          if (is.null(endLay)){
            warning(paste0("The disturbance end layer (", Sector," for ", disturbanceOrigin,
                           ") couldn't be found in the disturbanceList.  It is possible that the",
                           " disturbance is not available in the study area. Returning NULL."), 
                    immediate. = TRUE)
            return(NULL)
          }
        } # Bug catch for mismatches between disturbanceOrigin and the layers
      }
      if (is(oriLay, "RasterLayer"))
        oriLay <- rast(oriLay) # Using terra is faster
      # 4. Merge the layers and then connect the features
      if (is(oriLay, "SpatRaster")){
        # 4.1. Convert the raster to poly
        # Need to convert all that is zero into NA otherwise 2 polygons are created (one for 
        # 0's and one for 1's):
        oriLay[oriLay != 1] <- NA
        oriLayVect <- terra::as.polygons(oriLay)
      } else oriLayVect <- oriLay
      if (any(NROW(oriLayVect) == 0)){
        # If the disturbance doesn't exist, we need to return NULL
        return(NULL)
      }
      # Remove any polygons that are NA
      # Identify which rows have NA in first column and remove them
      whichRowsNA <- which(is.na(oriLayVect[,1]))
      RowsNotNA <- !is.na(oriLayVect[,1])
      if (length(whichRowsNA) > 0){
        message(paste0("Origin layer has disturbances that are NA. Removing ", 
                       length(whichRowsNA),
                       " ", geomtype(oriLayVect), "..."))
        oriLayVect <- subset(oriLayVect, subset = RowsNotNA)
      }
      # Here I should go over each individual disturbance and connect it iteratively.
      # This might eliminate the problem with so many lines at the same time very close to 
      # each other. I probably need to:
      # 1. Make a for loop over the number of lines, 
      # 2. select just one feature at a time
      # oriLayVectSingle <- st_as_sf(oriLayVect) # For sf and st_connect. 
      # 3. use the st_connect or terra::nearest on it # Changed to terra's native nearest. Much faster!
      classEndLay <- na.omit(unique(endLay[["Class"]])) 
      if (all(!is.null(connectingBlockSize),
              NROW(oriLayVect) > connectingBlockSize)){
        message(paste0("connectingBlockSize is not NULL and connecting layer has ", 
                       NROW(oriLayVect), " rows. Applying blocking technique to speed up",
                       "disturbance generation type connecting. ",
                       " If too many lines are connecting from the same place, decrease the",
                       "connectingBlockSize or ",
                       "set it to NULL, which will improve final result, but needs considerable",
                       "more time to run."))
        blockFullList <- 1:NROW(oriLayVect)
        amountBlocks <- ceiling(NROW(oriLayVect)/connectingBlockSize)
        blockList <- list()
        for (i in 1:amountBlocks) {
          if (length(blockFullList) <= connectingBlockSize){
            blockList[[i]] <- blockFullList
          } else {
            blockList[[i]] <- sample(x = blockFullList, 
                                     size = connectingBlockSize, replace = FALSE)
            # Update the list
            blockFullList <- setdiff(blockFullList, blockList[[i]])
          }
        }
        names(blockList) <- paste0("block_", 1:amountBlocks)
        for (i in names(blockList)){
          whichVecs <- blockList[[i]]
          message(paste0("Connecting ", Sector," for year ", currentTime,
                         ": ", i ," of ", length(blockList),
                         " (", 100*round(length(whichVecs)/NROW(oriLayVect),3),"%)"))
          connectedOne <- terra::nearest(oriLayVect[whichVecs, ], endLay, 
                                         pairs = FALSE, 
                                         centroids = TRUE, 
                                         lines = TRUE)
          # 4. Update the endLayer with the new connection
          if (!is(connectedOne, "SpatVector")){
            {
              cExt <- ext(connected)
              if (any(is.nan(cExt[1:4]))){
                print("Point44.2")
                browser()
              }
            }
            connectedOne <- terra::vect(connectedOne)
            {
              cExt <- ext(connected)
              if (any(is.nan(cExt[1:4]))){
                print("Point44.1")
                browser()
              }
            }
            
          }
          if (exists("connected")){
            {
              cExt <- ext(connected)
              if (any(is.nan(cExt[1:4]))){
                print("Point4.00")
                browser()
              }
            }
            connected <- rbind(connected, connectedOne)
            {
              cExt <- ext(connected)
              if (any(is.nan(cExt[1:4]))){
                print("Point4.0")
                browser()
              }
            }
            
          } else {
            {
              cExt <- ext(connected)
              if (any(is.nan(cExt[1:4]))){
                print("Point4.21")
                browser()
              }
            }
            connected <- connectedOne
            {
              cExt <- ext(connected)
              if (any(is.nan(cExt[1:4]))){
                print("Point4.31")
                browser()
              }
            }
          }
          endLay <- rbind(endLay, connectedOne)
        }          
      } else {
        if (any(useRoadsPackage, maskWaterAndMountainsFromLines)){
          Require("geodata")
          DEMraw <- elevation_30s(country = "CAN", path = Paths[["inputPath"]])
          DEM <- postProcessTo(from = DEMraw, to = rasterToMatch, cropTo = studyArea, 
                               projectTo = rasterToMatch, maskTo = studyArea, 
                               writeTo = file.path(Paths[["inputPath"]], 
                                                   paste0("DEM_", 
                                                          reproducible::.robustDigest(studyArea))))
        }
        if (useRoadsPackage){
          # Need to create `connected`
          # But roads package is SUPER slow... 
          tic("Time Elapsed for Road Building: ")
          connected <- projectRoads(landings = st_as_sf(oriLayVect), # points (works with lines?) where to join to the roads (landings)
                                    weightRaster = DEM,     # raster with zeros on roads and "weight" (DEM?)
                                    roads = st_as_sf(endLay),  #an existing road network as lines
                                    plotRoads = TRUE,
                                    roadsInWeight = FALSE)
          toc()
        } else {
          if (maskWaterAndMountainsFromLines){
            if (is.null(featuresToAvoid)){
              message(paste0("Features to avoid (i.e., lakes, rivers, wetland and mountain tops higher than ",
                             altitudeCut,"m) is NULL.",
                             "The module will create those."))
              # Alternative to roads package:
              # # TRY IMPLEMENTING A VERSION OF THE LINES WHICH AVOID CERTAIN POLYGONS?
              # # THEN I CAN MAKE RIVERS, LAKES, WETLAND AND HIGH MOUNTAIN DISAPPEAR FROM THE MAP, 
              # AND THUS BE AVOIDED BY LINES BUILDING
              waterraw <- landcover("water", path = Paths[["inputPath"]])
              water <- postProcessTo(from = waterraw, to = rasterToMatch, cropTo = studyArea, 
                                     projectTo = rasterToMatch, maskTo = studyArea, 
                                     writeTo = file.path(Paths[["inputPath"]], 
                                                         paste0("water_", 
                                                                reproducible::.robustDigest(studyArea))))
              wetlandsraw <- landcover("wetland", path = Paths[["inputPath"]])
              wet <- postProcessTo(from = wetlandsraw, to = rasterToMatch, cropTo = studyArea, 
                                   projectTo = rasterToMatch, maskTo = studyArea, 
                                   writeTo = file.path(Paths[["inputPath"]], 
                                                       paste0("wet_", 
                                                              reproducible::.robustDigest(studyArea))))
              # Now we put all layers together and exclude the pixels for building the roads (but without modifying the final maps!)
              DEM[DEM[] < altitudeCut] <- 0
              DEM[DEM[] >= altitudeCut] <- 1
              DEM[wet[] > 0.8] <- 1
              DEM[water[] > 0.8] <- 1
              featuresToAvoid <- copy(DEM)
              featuresToAvoid[featuresToAvoid == 0] <- NA
              featuresToAvoid[is.na(rasterToMatch)] <- 1
            }
          }
          # Here comes what is below...
          for (i in 1:NROW(oriLayVect)){
            if(i%%100==0)
              message(paste0("Connecting ", Sector," for year ", currentTime,
                             ": ", i ," of ", NROW(oriLayVect),
                             " (", 100*round(i/NROW(oriLayVect),4),"%)"))
            endLayCropped <- terra::crop(endLay, oriLayVect[i, ])  
            if (!NROW(endLayCropped) == 0){ # If you still have overlap with the crop, make sure they really overlap
              intersectedOriEnd <- terra::intersect(endLayCropped, oriLayVect[i, ])
              if (!NROW(intersectedOriEnd) == 0){
                message(paste0("Connecting ", Sector," for year ", currentTime,
                               ": ", i ," of ", NROW(oriLayVect),
                               " OVERLAPS with final layer, moving to next polygon"))
                next
              }
            }
            if (any(is.na(as.vector(ext(oriLayVect[i, ]))))){
              print("oriLayVect extent is NA. Fix didn't work. Debug.")
              browser()
            } 
            if (maskWaterAndMountainsFromLines){
              Require::Require("spaths")
              px <- extract(rasterToMatch, ext(oriLayVect[i, ]), cells = TRUE)
              oRas <- copy(rasterToMatch)
              oRas[] <- 0
              oRas[px[["cell"]]] <- 1
              cEndLay <- terra::nearest(oriLayVect[i, ], endLay, 
                                        pairs = FALSE, 
                                        centroids = TRUE, 
                                        lines = FALSE)
              message(paste0("Finding shortest path for feature ", i, " of ", NROW(oriLayVect),
                             ": ", 100*round(i/NROW(oriLayVect), 4),"% of ", 
                             Sector," (",disturbanceOrigin,
                             ") for year ", currentTime))
              connectedOne <- tryCatch({
                spaths::shortest_paths(rst = oRas, origins = oriLayVect[i, ],
                                       destinations = cEndLay,
                                       update_rst = as.polygons(featuresToAvoid), 
                                       output = "lines", show_progress = FALSE)
              }, error = function(e){
                print("Caught error in shortest path.")
                browser()
              })
              if (any(is.nan(ext(connectedOne[connectedOne$layer > 0,])[1:4]))){
                print(paste0("NaN being created likely because features to avoid",
                             "are within all possible solutions. Trying to fix it..."))
                connectedOne <- terra::aggregate(connectedOne)
              } else {
                connectedOne <- connectedOne[connectedOne$layer > 0,]
              }
              
              # message("***** END Finding shortest paths *******")
            } else {
              message(paste0("Connecting feature ", i, " of ", NROW(oriLayVect),
                             ": ", 100*round(i/NROW(oriLayVect), 4),"% of ", 
                             Sector," (",disturbanceOrigin,
                             ") for year ", currentTime))
              connectedOne <- terra::nearest(oriLayVect[i, ], endLay, 
                                             pairs = FALSE, 
                                             centroids = TRUE, 
                                             lines = TRUE)
            }
            # 4. Update the endLayer with the new connection
            if (!is(connectedOne, "SpatVector")){
              {
                cExt <- ext(connectedOne)
                if (any(is.nan(cExt[1:4]))){
                  print("Point2.1")
                  browser()
                }
              }
              connectedOne <- terra::vect(connectedOne)
              {
                cExt <- ext(connectedOne)
                if (any(is.nan(cExt[1:4]))){
                  print("Point2.2")
                  browser()
                }
              }
              
            }
            if (exists("connected")){
              {
                cExt <- ext(connected)
                if (any(is.nan(cExt[1:4]))){
                  print("Point0000")
                  browser()
                }
              }
              connected <- rbind(connected, connectedOne)
              {
                cExt <- ext(connectedOne)
                if (any(is.nan(cExt[1:4]))){
                  print("Point1111")
                  browser()
                }
              }
              {
                cExt <- ext(connected)
                if (any(is.nan(cExt[1:4]))){
                  print("Point2222")
                  browser()
                }
              }
              
            } else {
              connected <- connectedOne
            }
            endLay <- rbind(endLay, connectedOne)
          }
          message(paste0("For loop for ", Sector," -- ", disturbanceEnd, " finished."))
        }
        
      }
      if (!exists("connected")){ # It's likely because the one or few created disturbances were 
        # overlapping with their final connection. This means, connected doesn't exist and we need to 
        # return(NULL)
        return(NULL)
      }
      connected[["Class"]] <- na.omit(classEndLay)
      Lay <- list(connected)
      names(Lay) <- disturbanceEnd
      if (disturbanceRateRelatesToBufferedArea){
        if (length(disturbanceEnd)>1){
          print("Size > 1. Debug")
          browser()
        } 
        LayBuff <- terra::buffer(Lay[[disturbanceEnd]], width = 500) # Need to aggregate to avoid double counting!
        LayBuff <- terra::aggregate(x = LayBuff, dissolve = TRUE)
        currArea <- terra::expanse(LayBuff, unit = "m", transform = FALSE) # NEEDS TO BE METERS. USED BELOW
        message(paste0("Buffered (500m) area for ", Sector, 
                       " -- ", disturbanceEnd, ": ", round(currArea/1000000, 2), " Km2"))
        message(paste0("Percentage of current area: ", round(100*((currArea/1000000)/totalstudyAreaVAreaSqKm), 3), "%."))
      }
      
      return(Lay)
    })
    names(Connected) <- dPar[["dataName"]]
    # Mash together what has the same name
    # Which Sector has layers that are the same?
    whichSectorHasMoreThanOne <- unique(names(Connected)[duplicated(names(Connected))])  
    whichSectorHasOne <- setdiff(names(Connected), whichSectorHasMoreThanOne)
    
    if (length(whichSectorHasMoreThanOne) > 1){
      stop(paste0("There are at least two sectors with more than one layers to merge, which is",
                  "not been implemented. Please adapt the code in generateDisturbanceShp.R to ",
                  "allow for it."))
    }
    toMerge <- unlist(Connected[which(names(Connected) %in% whichSectorHasMoreThanOne)], 
                      recursive = TRUE)
    names(toMerge) <- NULL
    if (length(toMerge) == 0){
      warning(paste0("There is one sector with multiple layers to merge (i.e.,",
                     whichSectorHasMoreThanOne,"), but they are all NULL.",
                     "This has not been extensively tested and might fail!"),immediate. = TRUE)
      merged <- list(NULL)
    } else {
      merged <- list(do.call(rbind, toMerge)) 
    }
    names(merged) <- names(Connected[whichSectorHasMoreThanOne]) # Second level
    merged <- list(merged) # First level
    names(merged) <- whichSectorHasMoreThanOne # First level
    
    # Replace the merged class in the Connected object
    updatedToMerge <- Connected[names(Connected) %in% whichSectorHasOne]
    Connected <- c(updatedToMerge, merged)
  }

  # LAST: Updating and returning
  #############################################
  
  # Put all updated layers in a list to return
  individuaLayers <- c(Enlarged, Generated, Connected)
  
  if (is.null(individuaLayers)){
    stop("The study area has no potential for disturbances. Please enlarge or change the area.")
  }
  # Now merge all disturbances (raster format) to avoid creating new disturbances where it has 
  # already been disturbed
  currentDisturbance <- copy(individuaLayers)
  
    curDistRas <- lapply(1:length(currentDisturbance), function(index1){

    SECTOR <- names(currentDisturbance)[index1]
    if (is.null(currentDisturbance[[index1]])){ 
      # If NULL on the outter side, return NULL
      curDistRas <- NULL
    } else {
      maxLeng <- length(currentDisturbance[[index1]])
    curDistRas <- lapply(1:maxLeng, function(index2){
      Class <- names(currentDisturbance[[index1]])[index2]
      ras <- currentDisturbance[[index1]][[index2]]
      isNULLras <- is.null(ras)
      isMAX0ras <- if (is(ras, "SpatVector"))
        length(ras) == 0 else
          max(ras[], na.rm = TRUE) == 0
      if (any(isNULLras, isMAX0ras)){
        message(paste0(Class, " of ", SECTOR, " is NULL. Returning NULL"))
        return(NULL)
      }
      if (class(ras) %in% c("RasterLayer", "SpatRaster")){
        message(paste0("Converting ", Class, " of ", SECTOR, " from raster to vector..."))
        if (!is(ras, "SpatRaster"))
          ras <- rast(ras)
        ras[ras != 1] <- NA
        currentDisturbanceLay <- terra::as.polygons(ras)
        return(currentDisturbanceLay)
      } else {
        if (!is(ras, "SpatVector")){
          message(paste0(Class, " of ", SECTOR, " is not SpatVector, converting to it..."))
          ras <- rast(ras)
        } else {
          message(paste0(Class, " of ", SECTOR, " is either already a SpatVector or NULL, skipping conversion..."))
        }
        return(ras)
      }
    })
    if (!is.null(unlist(curDistRas))){
      if (length(curDistRas) == length(names(currentDisturbance[[index1]]))){
        names(curDistRas) <- names(currentDisturbance[[index1]])
      } else {
        warning("Rasters are not NULL, but the number of rasters doesn't match the number of names. Please debug.", 
                immediate. = TRUE)
        browser()
      }
    } else {
      message(paste0("The current disturbance raster for ", SECTOR, " is NULL. Returning without internal names."))
    }
    }
    return(curDistRas)
  })
  if (length(curDistRas) == length(names(currentDisturbance))){
    names(curDistRas) <- names(currentDisturbance)
  } else {
    warning("The number of rasters doesn't match the number of names. Please debug.", 
            immediate. = TRUE)
    browser()
  }
    # Cleanup the NULLS!
    curDistRas <- cleanupList(curDistRas, 
                              outter = TRUE, 
                              inner = TRUE, 
                              cleanEmpty = TRUE, 
                              nullEmpty = FALSE)

  ########################### FINAL LAYERS ###########################
  
  # individualLayers: mixed rasters and vectors, coming directly from the disturbances generated
  
  # curDistRas: all vectors, unbuffered new disturbances, no old ones here except for seismic lines (if enlarging) and settlements
  #             While seismic lines are not lines anymore (they got buffered to match real width), we still have lines
  #             belonging to roads and pipelines.     
  
  ########################### FINAL LAYERS ###########################
  
  # Now I want to test it with 500m buffer --> For knowledge. The layers below are not 
  # returned!!
  
  if (all(disturbanceRateRelatesToBufferedArea,
          checkDisturbancesForBuffer)){ 
    #TODO I don't think this is working. New disturbances appear as much smaller than the original ones. 
    # Maybe some disturbances (i.e., seismic lines) in the new layer are not including the existing "old"  
    # layer? I am pretty sure that this is the problem... We should fix this!
    curDistVcs <- unlist(curDistRas)
    
    curDistVcsAll <- lapply(names(curDistVcs), function(eachVectNm){
      eachVect <- curDistVcs[[eachVectNm]]
      if (length(eachVect) == 0) return(NULL)
      message(paste0("Buffering and/or merging polygons for ", eachVectNm))
      buffVect <- terra::buffer(eachVect, width = 500)
      buffVect <- terra::aggregate(buffVect, dissolve = TRUE)
      return(buffVect)
    })
    newDisturbanceLayers <- do.call(rbind, curDistVcsAll)
    newDisturbanceLayers <- terra::aggregate(newDisturbanceLayers, dissolve = TRUE)
    
    # Now I get the previous layers and do the same
    oldDisturbanceLayers <- createBufferedDisturbances(disturbanceList = disturbanceList, 
                                                       bufferSize = 500,
                                                       rasterToMatch = rasterToMatch,
                                                       studyArea = studyArea,
                                                       currentTime = currentTime,
                                                       convertToRaster = FALSE)
    
    newDisturbanceLayers$totAreaKm2 <- terra::expanse(newDisturbanceLayers, 
                                                      unit = "km", transform = FALSE) 
    oldDisturbanceLayers$totAreaKm2 <- terra::expanse(oldDisturbanceLayers, 
                                                      unit = "km", transform = FALSE) 
    
    message(paste0("Buffered 500m (polygons) old disturbance percent of the area (",currentTime,"): ", 
                   round(100*(oldDisturbanceLayers$totAreaKm2/totalstudyAreaVAreaSqKm), 3), "%."))
    
    message(paste0("Buffered 500m (polygons) new disturbance percent of the area (",currentTime,"): ", 
                   round(100*(newDisturbanceLayers$totAreaKm2/totalstudyAreaVAreaSqKm), 3), "%."))
    
    message(paste0("This means a total distubance growth rate of ", 
                   (round(100*(newDisturbanceLayers$totAreaKm2/totalstudyAreaVAreaSqKm), 3)-round(100*(oldDisturbanceLayers$totAreaKm2/totalstudyAreaVAreaSqKm), 3))/unique(disturbanceParameters[["disturbanceInterval"]]),
                   " per year."))
  }
  
  ### RETURN!
  if (firstTime){
    seismicLinesFirstYear <- vect(file.path(outputsFolder, 
                                                paste0("seismicLinesYear", 
                                                       currentTime, 
                                                       "_",
                                                       studyAreaHash,
                                                       ".shp")))
  } else {
    seismicLinesFirstYear <- NULL
  }

  individuaLayers <- cleanupList(individuaLayers, 
                              outter = TRUE, 
                              inner = TRUE, 
                              cleanEmpty = TRUE, 
                              nullEmpty = FALSE)
  
  return(list(individuaLayers = individuaLayers, 
       currentDisturbanceLayer = curDistRas,
       seismicLinesFirstYear = seismicLinesFirstYear))
}
