generateDisturbances <- function(disturbanceParameters,
                                 disturbanceList,
                                 rasterToMatch,
                                 studyArea,
                                 fires,
                                 currentTime,
                                 growthStepEnlargingPolys,
                                 growthStepEnlargingLines,
                                 currentDisturbanceLayer,
                                 runName){

  # Extracting layers from previous ones
  # Total study area
  studyArea <- vect(studyArea)
  studyArea <- project(x = studyArea, y = terra::crs(rasterToMatch))
  uniStudyArea <- terra::aggregate(studyArea)
  totalstudyAreaVAreaSqKm <- terra::expanse(uniStudyArea, unit = "km")
  totalstudyAreaVAreaSqm <- terra::expanse(uniStudyArea, unit = "m")
  
  # First, mask the current disturbances on the potential layer so we don't choose them again
  if (is.null(currentDisturbanceLayer)){
    message(paste0("currentDisturbanceLayer is NULL. This is likely the first year of the",
                   " simulation. Creating current disturbed raster..."))
    # Current layer is null, which indicates this might be the first year of the simulation. I
    # will need to create it then, based on the existing disturbances.
    # Get all current disturbances
    unDL <- unlist(disturbanceList)
    allCurrentDistLayNames <- names(unDL)[!grepl(x = names(unDL), pattern = "potential")]
    currDist <- unDL[names(unDL) %in% allCurrentDistLayNames]
    # Combine all layers, both polygons and lines, separtely
    linesLays <- currDist[sapply(currDist, geomtype) == "lines"]
    polyLays <- currDist[sapply(currDist, geomtype) == "polygons"]
    # Buffer lines at the resolution of the rasterToMatch to convert them to polys. 
    # This is not the ideal resolution, but as we are NOT deriving the buffered disturbances 
    # from this layer but just using it to avoid choosing pixels already chosen, then it should 
    # be ok at 250m.  
    # NOTE: We devide the res by 2 so the total around is exactly the resolution of a pixel,
    # as the buffer uses the width on all sides. 
    
    #TODO Think this through on diagonals: Are these also captured?
    # If not, we would need actually the half of the diagonal as buffer size, not side. 
    # ~~~ IS THIS REALLY NEEDED ? ~~~ 
    # (sqrt(2)*res(rasterToMatch)[1])/2
    
    linesLays <- unlist(lapply(linesLays, function(X){
      terra::buffer(x = X, width = res(rasterToMatch)[1]/2) 
      # Here it needs to be the raster res, not the real one because we are not dealing with 
      # real area to be calculated, but just the buffering to avoid specific pixels to be 
      # selected 
    }))
    linesAndPolys <- do.call(what = c, list(polyLays, linesLays))
    names(linesAndPolys) <- NULL
    allLays <- do.call(rbind, linesAndPolys)
    # Now fasterize all polys
    allLaysSF <- sf::st_as_sf(x = allLays)
    allLaysSF$Polys <- 1
    polysField <- "Polys"
    allLaysRas <- fasterize::fasterize(sf = st_collection_extract(allLaysSF, "POLYGON"),
                                       raster = rasterToMatch, 
                                       field = polysField)
  } else {
    allLaysRas <- currentDisturbanceLayer
    # Create allLaysRas, which is the disturbance raster with all disturbances rasterized, with
    # value == 1 for the disturbed pixels and non-disturbed pixels are NA
  }
  
  # FIRST: Enlarging
  #############################################
  message(crayon::white("Enlarging disturbances..."))
  dPar <- disturbanceParameters[disturbanceType  %in% "Enlarging", ]
  whichSector <- dPar[["dataName"]]
  Enlarged <- lapply(whichSector, function(Sector) {
    whichOrigin <- dPar[dataName == Sector, disturbanceOrigin]
    if (length(whichOrigin) != 1) {
      print("Repeated origins and sector. Debug")
      browser()
    } # Bug Catch
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
      # NOTE: The value from the table should already be in percent!!
      # We need to assume a form for polygons and for lines.
      # Polygons: CIRCLE
      # Lines: RECTANGLES
      # Then we need the formula to calculate the buffer based on the increase in
      # area size (parameter from table) using a normal
      originalForm <- geomtype(Lay)
      if (originalForm %in% c("lines", "points")){
        # Check if we have the resolution
        if ("resolutionVector" %in% dPar){
          RES <- dPar[index, resolutionVector]/2
        } else RES <- NULL
        if (is.null(RES)){
          message(crayon::red(paste0("resolutionVector in disturbanceParameters",
                                     " table was not supplied for a lines file vector. Will ",
                                     "default to 15m total (7.5m in each direction).",
                                     "If this is not correct, please provide the ",
                                     "resolution used for developing the layer ", ORIGIN,
                                     " (", Sector,")")))
          RES <- 7.5
        }
        # NOTE: The buffer value is applied on all sides of the shapefile / line
        # This means that a 15m buffer ends up being a 30 m increase in the line
        # in all directions. If the original resolution was 15, then we need to divide it by 2,
        # so we get in the end a TOTAL increase in of 15 m.
        Lay <- terra::buffer(x = Lay, RES)
      }
        currArea <- terra::expanse(Lay, unit = "m")
        # The dParOri[["disturbanceRate"]] (or disturbRate) is (by default) the calculated difference between
        # the 2010 and 2015 disturbance divided by the total amount of years, over the entire area
        # Therefore, it refers to the total study area size. 
        disturbRate <- dParOri[["disturbanceRate"]]
        # Disturbance area to add to current:
        disturbAreaSqKm <- disturbRate*totalstudyAreaVAreaSqKm
        # Total increase in area to be distributed across all disturbances:
        disturbAreaSqM <- disturbAreaSqKm * 1000000
        # What is the % to growth for each disturbance:
        totalCurrDistArea <- sum(currArea)
        totalPercGrowth <- disturbAreaSqM/totalCurrDistArea
        
        # We will do the buffer iteratively as there are no great solutions for getting the correct 
        # value by calculations because the polygon's forms is not known.
        # While totalPercGrowthAchieved is below totalPercGrowth, we will keep increasing the 
        # buffer value and recalculating it. Past trials are first below:
        
        totalPercGrowthAchieved <- 0
        expandTo <- 0
        iter <- 1
        tictoc::tic("Total elapsed time for calculating disturbance percentage: ")
        while (totalPercGrowthAchieved < totalPercGrowth) {
          if (iter == 1)
            message(paste0("Calculating disturbance percentage for ",
                           ORIGIN, " (iteration ", iter,")"))
          if (iter %% 100 == 0) {
            message(paste0("Calculating disturbance percentage for ",
                         ORIGIN, " (iteration ", iter,")"))
          }
          if (originalForm %in% c("lines", "points")){
            expandTo <- expandTo + growthStepEnlargingLines
          } else {
            expandTo <- expandTo + growthStepEnlargingPolys 
          }
          LayUpdated <- terra::buffer(x = Lay, width = expandTo)
          totalCurrDistArea <- sum(currArea)
          futureArea <- terra::expanse(LayUpdated, unit = "m")
          totalFutureDistArea <- sum(futureArea)
          totalPercGrowthAchieved <- (totalFutureDistArea-totalCurrDistArea)/totalCurrDistArea
          iter <- iter + 1 
        }
        tictoc::toc()
        
        cat(crayon::yellow(paste0("Percentage difference between mean disturbance rate",
                                  " and mean change in total area for sector ",
                                  crayon::red(Sector), " disturbance type ", crayon::red(ORIGIN), ": \n",
                                  crayon::red(format(100*round((totalPercGrowthAchieved - 
                                                                    totalPercGrowth)/totalPercGrowth, 4), 
                                                       scientific = FALSE), " %."),
                                  "\nThis value can help interpret how close the increased areas ",
                                  "are to the expected increases. The ideal value is 0.\n")))
        
        cat(paste0(Sector, 
                   " ", 
                   ORIGIN, 
                   " ",
                   currentTime,
                   " ",
                   format(100*round((totalPercGrowthAchieved - totalPercGrowth)/totalPercGrowth, 8), 
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
  
  # SECOND: Generating
  #############################################
  message(crayon::white("Generating disturbances..."))
  dPar <- disturbanceParameters[disturbanceType  %in% "Generating", ]
  whichSector <- dPar[["dataName"]]
  Generated <- lapply(whichSector, function(Sector) {
    whichOrigin <- dPar[dataName == Sector, disturbanceOrigin]
    if (length(whichOrigin) != 1) {
      print("Repeated origins and sector. Debug")
      browser()
    } # Bug Catch
    updatedL <- lapply(whichOrigin, function(ORIGIN) {
      dParOri <- dPar[dataName == Sector & disturbanceOrigin == ORIGIN,]
      Lay <- disturbanceList[[Sector]][[ORIGIN]]
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
      
      # 2. Fasterize it
      potLaySF <- sf::st_as_sf(x = potLay)
      potField <- dParOri[["potentialField"]]
      
      if (any(is.na(potField), 
              potField == "",
              is.null(potField))){ # If NA, it doesn't matter, but need to create a 
        # field so not the whole thing becomes one big 1 map
        potLaySF$Potential <- 1
        potField <- "Potential"
        message(crayon::yellow(paste0("Potential field not declared for ", ORIGIN, " (",
                                      Sector, "). Considering all polygons as likely to develop ",
                                      "generated disturbance.")))
      }
      
      # Update with fire layer if forestry
      if (Sector == "forestry"){

        message(paste0("Generating disturbance for forestry. Updating potential layer for ",
                       "occurred fires and currently productive forest"))
        
        # First: Select only productive forests
        # potLaySF <- potLaySF[potLaySF$ORIGIN < (currentTime - 50), ] # For some weird reason, this doesn't work anymore... sigh.
        potLaySF <- subset(potLaySF, potLaySF$ORIGIN < (currentTime - 50))
        potLayF <- fasterize::fasterize(sf = st_collection_extract(potLaySF, "POLYGON"),
                                        raster = rasterToMatch, field = potField)
        # Second: remove fires
        if (!is.null(fires))
          potLayF[fires[] == 1] <- NA
        # Third: give preference to cutblocks that are closer to current cutblocks
        # Using terra is quicker!
        # # WITH NWT DATA, THE FOLLOWING LINES ARE USELESS AS ALL FOREST IS CLOSE ENOUGH
        # # WE CAN IMPROVE THAT, HOWEVER, BY USING BUFFERING METHODS
        # potLayFt <- terra::rast(potLayF)
        # tictoc::tic("Distance raster time elapsed: ")
        # distRas <- terra::distance(potLayFt) # Distance raster time elapsed: : 8348.948 sec elapsed
        # toc()
        # distRas2 <- (distRas-maxValue(distRas))*-1
        # distRas2[is.na(potLayF)] <- NA
      } else {
        potLayF <- fasterize::fasterize(sf = st_collection_extract(potLaySF, "POLYGON"),
                                        raster = rasterToMatch, field = potField)
      }
      
      # 3. Check how many pixels we need to choose, and in which sizes based on the table
      Rate <- dParOri[["disturbanceRate"]]
      totNPix <- sum(rasterToMatch[], na.rm = TRUE)
      # Here we multiply the rate by the total number of pixels by the interval we are using to 
      # generate the disturbances (as rate is passed as yearly rate!), to know how many pixels we 
      # are expected to disturb in the given interval for this study area, given the rate applied
      expectedDistPixels <- Rate*totNPix*dParOri[["disturbanceInterval"]]
      if (expectedDistPixels < 1){
        # If the number of expected pixels for that given year is smaller than 1, then we need to
        # apply a probability of that disturbance happening. 
        message(crayon::red("The number of expected disturbed pixels for step ", 
                            currentTime, " is less than 1. The function will apply a probability ",
                            "of this disturbance happening, with the given size expected."))
        totPixToChoose <- rbinom(1, 1, expectedDistPixels)
      } else {
        nPixChosenTotal <- 0
        chosenDistrib <- numeric()
        IT <- 1
        while (round(nPixChosenTotal, 0) < round(expectedDistPixels, 0)){
          if (IT == 1){
            print(paste0("Calculating total disturbance size for ", 
                         ORIGIN, " (iteration ", IT, ", ", 
                         round(100*(round(nPixChosenTotal, 0)/round(expectedDistPixels, 0)), 2),"% achieved)"))
          }
          if (IT %% 10 == 0){
            print(paste0("Calculating total disturbance size for ", 
                         ORIGIN, " (iteration ", IT, ", ", 
                         round(100*(round(nPixChosenTotal, 0)/round(expectedDistPixels, 0)), 2),"% achieved)"))
          }
          # This is the size in meter square that each disturbance should have
          Size <- round(eval(parse(text = dParOri[["disturbanceSize"]])), 0)
          # This is the size as the number of pixels to be chosen
          # This below doesn't work well as the resolution doesn't account for the geodesic form of 
          # the planet:
          # nPixChosen <- round(sqrt(Size)/res(rasterToMatch)[1], 0)
           
          # This solution accounts more or less for this!
          # totNPix is the number of pixels that correspont to totalstudyAreaVAreaSqKm
          # As the calculation based on pixel resolution is not good, we can try a better approximation
          # by using the study area shapefile calculated with terra. It is not the best as southern pixels
          # are in reality smaller than the northern ones, but it is a good approximation for now: 
          calculatedPixelSizem2 <- totalstudyAreaVAreaSqm/totNPix # in m2
          calculatedPixelSizeLength <- sqrt(calculatedPixelSizem2)
          # Here I am getting the size of the 
          nPixChosen <- round(sqrt(Size)/calculatedPixelSizeLength, 0)
          nPixChosenTotal <- nPixChosenTotal + nPixChosen
          chosenDistrib <- c(chosenDistrib, nPixChosen)
          IT <- IT + 1
        }
        message(paste0("Disturbance size successfully calculated for ", ORIGIN))
        totPixToChoose <- chosenDistrib[chosenDistrib > 0] # Need to clean up as sometimes we
        # get zero pixels if the size of each disturbance is not too big
      }
      # 4. Select which pixels are the most suitable and then get the number of neighbors based on 
      # the needed sizes
      # Can't forget to take pixels that already have disturbance out of the potential availability!
      # Otherwise it might choose pixels to disturb in existing disturbed ones
      # Second, bring the raster with masked disturbances into a table to choose the pixels
      potLayF[allLaysRas[] == 1] <- NA
      potentialDT <- na.omit(data.table(pixelID = 1:ncell(potLayF),
                                vals = terra::values(potLayF)))
      # Third, subset the best pixels
      bestPotential <- potentialDT[vals == max(unique(potentialDT[["vals"]])), pixelID]
      # Then, randomly, choose the pixels. Here need to pay attention to choose one per "development"
      # and only then use the adjacent to get to the correct sizes.
      # length(totPixToChoose) --> is the number of pixels minus the ones that don't need adjacent. 
      # The totPixToChoose is then the sizes each of these need to have with adjacent cells
      # and sum(totPixToChoose) is the total number of pixels where disturbance will be generated, which
      # should match very closely to expectedDistPixels
      allSelected <- FALSE
      if (sum(totPixToChoose) > length(bestPotential)){
        # This means that the best potential is not enough. So we need to increase to the 2 best 
        # potentials.
        counter <- 1
        while (sum(totPixToChoose) > length(bestPotential)){
          classesToInclude <- 1:(counter+1)
          # subset the best pixels
          bestPotential <- potentialDT[vals %in% sort(decreasing = TRUE, 
                                                    unique(potentialDT[["vals"]]))[classesToInclude], 
                                       pixelID]
          counter <- counter+1
          if (all(classesToInclude %in% unique(potentialDT[["vals"]]))){
            message(crayon::red(paste0("All pixels available for ", ORIGIN, " have been selected. ",
                                       "This distubance will remain static until the end of the ",
                                       "simulation.")))
            pixNewDist <- bestPotential
            allSelected <- TRUE
            break
          }
        }
      } else {
        pixNewDist <- sample(x = bestPotential, 
                             size = length(totPixToChoose), 
                             replace = FALSE)
      }
      # We then need to identify the adjacent cells, if any of the values of totPixToChoose is > 1
      if (all(any(totPixToChoose > 1),
              !allSelected)){
        # Create whichPixelsChosen using the pixNewDist and adjacent cells
        pixDontNeedAdj <- pixNewDist[totPixToChoose == 1]
        pixNeedAdj <- pixNewDist[totPixToChoose > 1]
        howManyAdj <- totPixToChoose[totPixToChoose > 1]
        # Create an adjacent matrix
        # The distance depends on how many pixels I need
        whichAdj <- unique(howManyAdj)
        nbMatrices <- lapply(whichAdj, function(X){
          m <- focalWeight(x = rasterToMatch, d = X * res(rasterToMatch)[1], 
                      type = c('circle'))
          m[m > 0] <- 1
          return(m)
        })
        names(nbMatrices) <- paste0("NB_", whichAdj)
        print(ORIGIN)  
        allPixAdj <- sapply(whichAdj, function(totAdj){
          pix <- pixNeedAdj[which(howManyAdj == totAdj)]
          adjN <- totAdj - 1  # We need to exclude the focal cell
          nb <- nbMatrices[[paste0("NB_", totAdj)]]
          # A neighborhood matrix identifies the cells around each cell that are considered adjacent. 
          # The matrix should have one, and only one, cell with value 0 (the focal cell); 
          # at least one cell with value 1 (the adjacent cell(s)); All other cells are not considered 
          # adjacent and ignored.
          nb[nb == 0] <- -1
          nb[ceiling(ncol(nb)/2), ceiling(nrow(nb)/2)] <- 0
          Adj <- as.data.table(raster::adjacent(x = rasterToMatch, cells = pix,
                                  directions = nb, pairs = TRUE,
                                  include = FALSE, id = TRUE))
          # Now we sample the number of adjacents needed from the whole table, by ID (each pix that 
          # needs those adjacents)
          chosen <- Adj[,.SD[sample(.N, adjN)], by = "id"]
          return(c(unique(chosen[["from"]]), unique(chosen[["to"]])))
        })
        allPixAdj <- unique(as.numeric(unlist(allPixAdj)))
        whichPixelsChosen <- c(allPixAdj, pixDontNeedAdj)
      } else {
        whichPixelsChosen <- pixNewDist
      }
      # Now, we make a new disturbance layer. We need to keep track of the 
      # old one because we need to connect the new stuff. We can use RTM as a template.
      newDistLay <- copy(rasterToMatch)
      newDistLay[] <- newDistLay[]
      names(newDistLay) <- ORIGIN
      newDistLay[newDistLay == 1] <- 0
      newDistLay[whichPixelsChosen] <- 1
      cat(crayon::yellow(paste0("Percentage difference between mean disturbance rate",
                                " and mean change in total area for sector ",
                                crayon::red(Sector), " disturbance type ", crayon::red(ORIGIN), ": \n",
                                crayon::red(format(100*round((length(whichPixelsChosen) - 
                                                                expectedDistPixels)/expectedDistPixels, 4), 
                                                   scientific = FALSE), " %."),
                                "\nThis value can help interpret how close the increased areas ",
                                "are to the expected increases. The ideal value is 0.\n")))
      cat(paste0(Sector, 
                 " ", 
                 ORIGIN, 
                 " ",
                 currentTime,
                 " ",
                 format(100*round((length(whichPixelsChosen) - 
                                     expectedDistPixels)/expectedDistPixels, 8), 
                        scientific = FALSE)),
          file = file.path(Paths[["outputPath"]], paste0("PercentageDisturbances_", currentTime, 
                                                         "_", runName, ".txt")),
          append = TRUE, sep = "\n")
      
      
      return(newDistLay)
    })
    names(updatedL) <- whichOrigin
    return(updatedL)
  })
  names(Generated) <- whichSector

  # THIRD: CONNECTING (Use generated and table)
  #############################################
  message(crayon::white("Connecting disturbances..."))
  dPar <- disturbanceParameters[disturbanceType  %in% "Connecting", ]
  Connected <- lapply(1:NROW(dPar), function(index) {
    Sector <- dPar[index, dataName]
    disturbanceOrigin <- dPar[index, disturbanceOrigin]
    disturbanceEnd <- dPar[index, disturbanceEnd]
      # 1. Get the info on which layer to connect to which layer
      # 2. Get the origin layer (from Generated)
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
                stop(paste0("The disturbance origin layer (", Sector," for ", disturbanceOrigin,
                          ") couldn't be found in the disturbanceList. Please debug"))
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
              stop(paste0("The disturbance end layer (", Sector," for ", disturbanceOrigin,
                          ") couldn't be found in the disturbanceList. Please debug"))
            }
          } # Bug catch for mismatches between disturbanceOrigin and the layers
        }
          oriLay <- rast(oriLay) # Using terra is faster
        # 4. Merge the layers and then connect the features
        # 4.1. Convert the raster to poly
        # Need to convert all that is zero into NA otherwise 2 polygons are created (one for 
        # 0's and one for 1's):
        oriLay[oriLay != 1] <- NA
        oriLayVect <- as.polygons(oriLay, dissolve = FALSE)
        # Here I should go over each individual disturbance and connect it iteratively.
        # This might eliminate the problem with so many lines at the same time very close to 
        # each other. I probably need to:
        # 1. Make a for loop over the number of lines, 
        # 2. select just one feature at a time
        oriLayVectSingle <- st_as_sf(oriLayVect)
        # 3. use the st_connect on it
        classEndLay <- na.omit(unique(endLay[["Class"]]))
        for (i in 1:NROW(oriLayVectSingle)){
          connectedOne <- nngeo::st_connect(oriLayVectSingle[i, ], 
                                         sf::st_as_sf(endLay), 
                                         ids = NULL, 
                                         progress = TRUE)
          # 4. Update the endLayer with the new connection
          if (exists("connected")){
            connected <- rbind(connected, vect(connectedOne))
          } else {
            connected <- vect(connectedOne)
          }
          endLay <- rbind(endLay, vect(connectedOne))
        }
        connected[["Class"]] <- classEndLay
        Lay <- list(connected)
        names(Lay) <- disturbanceEnd
        return(Lay)
  })
  names(Connected) <- dPar[["dataName"]]
  # Mash together what has the same name
  # Which Sector has layers that are the same?

  nms <- names(unlist(lapply(Connected, names)))
  whichSectorHasOne <- unique(nms[ave(seq_along(nms), nms, FUN = length)==1])
  # To control for other cases:
  whichSectorHasMoreThanOne <- unique(nms[ave(seq_along(nms), nms, FUN = length)>1])
  if (length(whichSectorHasMoreThanOne) > 1){
    stop(paste0("There are at least two sectors with more than one layers to merge, which is currently not supported. ",
                "Please re-code (line 466 of generateDisturbances.R) to allow for this case"))
  }
  toMerge <- unlist(Connected)[which(!names(Connected) %in% whichSectorHasOne)]
  names(toMerge) <- NULL
  merged <- list(do.call(rbind, toMerge))
  nm <- unique(names(Connected[[whichSectorHasMoreThanOne]]))
  names(merged) <- nm # Second level 
  
  merged <- list(merged)
  names(merged) <- whichSectorHasMoreThanOne # First level

  # Replace the merged class in the Connected object
  updatedToMerge <- Connected[names(Connected) %in% whichSectorHasOne]
  Connected <- c(updatedToMerge, merged)
  
  # LAST: Updating and returning
  #############################################

  # Put all updated layers in a list to return
  individuaLayers <- c(Enlarged, Generated, Connected)
  
  # Now merge all disturbances (raster format) to avoid creating new disturbances where it has 
  # already been disturbed

  currentDisturbance <- individuaLayers
  curDistRas <- lapply(1:length(currentDisturbance), function(index1){
    SECTOR <- names(currentDisturbance)[index1]
    curDistRas <- lapply(1:length(currentDisturbance[[index1]]), function(index2){
      Class <- names(currentDisturbance[[index1]])[index2]
      ras <- currentDisturbance[[index1]][[index2]]
      if (!class(ras) %in% c("RasterLayer", "SpatRaster")){
        message(paste0("Converting ", Class, " of ", SECTOR, " to a raster..."))
        # Now we fasterize so this is used to exclude pixels that have already been disturbed
        # As some of these are lines, we will need to buffer them so they become polygons at
        # the original resolution they had in the dataset
        if (!is(ras, "SpatVector"))
          ras <- vect(ras)
        if  (terra::geomtype(ras) %in% c("points", "lines")) {
          if ("resolutionVector" %in% dPar){
            RES <- unique(dPar[disturbanceEnd == Class, resolutionVector])/2
          } else RES <- NULL
          if (any(is.na(RES),
                  is.null(RES))){
            message(crayon::red(paste0("resolutionVector in disturbanceParameters",
                                       " table was not supplied for a lines file vector. Will ",
                                       "default to 15m total (7.5m in each direction).",
                                       "If this is not correct, please provide the ",
                                       "resolution used for developing the layer ", Class,
                                       " (", SECTOR,")")))
            RES <- 7.5
          }
          currDistT <- terra::buffer(x = ras, width = RES)
        }  else {
          currDistT <- ras
        }
        currDistSF <- sf::st_as_sf(x = currDistT)
        currDistSF$disturbance <- 1
        fld <- "disturbance"
        currentDisturbanceLay <- suppressWarnings(fasterize::fasterize(sf = st_collection_extract(currDistSF, "POLYGON"),
                                                      raster = rasterToMatch, field = fld))
        # NOTE: 
        # Seems that the fasterize is not picking up roads from currDistSF, which makes sense considering they are 
        # much smaller than the resolution... However, this layer is just so we exclude the pixels that 
        # have been disturbed already. This means that lines should be fine, as we have space for more 
        # stuff to come in in a 250mx250m pixel that has 7.5 or 15m disturbed. For boo's models, we will 
        # use the full layers and buffer before we fasterize (PopGrowth) and do the calculations of decay 
        # in the polygon or at a higher resolution. So for now, this solution is fine.
        
        return(currentDisturbanceLay)
      } else {
        message(paste0(Class, " of ", SECTOR, " is already a raster, skipping conversion..."))
        return(ras)
      }
    })
    names(curDistRas) <- names(currentDisturbance[[index1]])
    return(curDistRas)
  })
  names(curDistRas) <- names(currentDisturbance)
  
  stk <- raster::stack(unlist(curDistRas))
  stkDT <- data.table(pixelIndex = 1:ncell(stk),
                      getValues(stk))
  stkDT[, currDisturbance := rowSums(.SD, na.rm = TRUE), .SDcols = names(stk)]
  
  newDisturbanceLayers <- setValues(x = raster(rasterToMatch), 
                                    values = stkDT[["currDisturbance"]])
  # Cleanup the zeros where is actually NA
  newDisturbanceLayers[is.na(rasterToMatch)] <- NA
  
  # Some disturbances happened in the same pixel (0.001% of the area), 
  # but this is such a small problem I will not deal with right now.   
  newDisturbanceLayers[newDisturbanceLayers[] > 1] <- 1
  
  if (!is.null(currentDisturbanceLayer)) {
    newDisturbanceLayers[currentDisturbanceLayer == 1] <- 1
  }
  
  ### RETURN!
  list(individuaLayers = individuaLayers, currentDisturbanceLayer = newDisturbanceLayers)
}
