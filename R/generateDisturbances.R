generateDisturbances <- function(disturbanceParameters,
                                 disturbanceList,
                                 rasterToMatch,
                                 currentTime,
                                 currentDisturbanceLayer){

  # FIRST: Enlarging
  #############################################
  message("Enlarging disturbances...")
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
      if (geomtype(Lay) == "polygons") {
        currArea <- terra::expanse(Lay)
        rCurrA <-  sqrt(currArea / pi) # Using a circle formula
        desiredArea <- currArea + (currArea * dParOri[["disturbanceRate"]])
        rDesiredArea <- sqrt(desiredArea / pi) # Using a circle formula
        expandTo <- rDesiredArea - rCurrA
        LayUpdated <- terra::buffer(x = Lay, width = expandTo)
        # To see how far we are from the form assumption we need to do
        print(paste0("Difference between disturbance rate and change in total area for sector ",
            Sector, " disturbance type ",ORIGIN))
        print(format(summary(dParOri[["disturbanceRate"]] - ((terra::expanse(LayUpdated) - currArea) /
                                                        currArea)), scientific = FALSE))
      } else {
        if (geomtype(Lay) == "lines"){
          currLength <- perim(Lay)
          # Check if we have the resolution
          RES <- dParOri[["resolutionVector"]]
          if (is.null(RES)){
            message(crayon::red(paste0("resolutionVector in disturbanceParameters",
                                       " table was not supplied for a lines file vector. Will ",
                                       "default to 15m. If this is not correct, please provide the ",
                                       "resolution used for developing the layer ", ORIGIN,
                                       " (", Sector,")")))
            RES <- 15
          }
          Lay <- terra::buffer(x = Lay, RES)
          currArea <- terra::expanse(Lay)
          # NOTE :: Noticed that using the circle formula we have a very small 
          #         difference in terms of average area gained and expected proportional gain
          #         So I will not worry now about updating to rectangle formula!  
          rCurrA <- sqrt(currArea / pi) # Using a circle formula
          desiredArea <- currArea + (currArea * dParOri[["disturbanceRate"]])
          rDesiredArea <- sqrt(desiredArea / pi) # Using a circle formula
          expandTo <- rDesiredArea - rCurrA
          LayUpdated <- terra::buffer(x = Lay, width = expandTo)
          # To see how far we are from the form assumption we need to do
          print(paste0("Difference between disturbance rate and change in total area for sector ",
                       Sector, " disturbance type ",ORIGIN))
          print(format(summary(dParOri[["disturbanceRate"]] - ((terra::expanse(LayUpdated) - currArea) /
                                                          currArea)), scientific = FALSE))
        } else {
          if (geomtype(Lay) == "points"){
            stop("Disturbances as points are still not implemented")
            # Need to add buffer around points. Maybe this can be done with lines?
          } else {
            stop("Enlarging vectors can for now only be polygons or lines")
          }
        }
      }
      return(LayUpdated)
    })
    names(updatedL) <- whichOrigin
    return(updatedL)
  })
  names(Enlarged) <- whichSector
  
  # SECOND: Generating
  #############################################
  message("Generating disturbances...")
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
      # 
      # 1. Get the potential layer
      potLay <- disturbanceList[[Sector]][[dParOri[["dataClass"]]]]
      # 2. Fasterize it
      potLaySF <- sf::st_as_sf(x = potLay)
      potField <- dParOri[["potentialField"]]
      
      if (any(is.na(potField), potField == "")){ # If NA, it doesn't matter, but need to create a field so not the whole thing
        # Becomes one big 1 map
        potLaySF$Potential <- 1
        potField <- "Potential"
        message(crayon::yellow(paste0("Potential field not declared for ", ORIGIN, " (",
                                      Sector, "). Considering all polygons as likely to develop ",
                                      "generated disturbance.")))
      }
      potLayF <- fasterize::fasterize(sf = st_collection_extract(potLaySF, "POLYGON"),
                                      raster = rasterToMatch, field = potField)

      # 3. Check how many pixels we need to choose, and in which sizes based on the table
      Rate <- dParOri[["disturbanceRate"]]
      totAreaPix <- sum(rasterToMatch[], na.rm = TRUE)
      # Here we multiply the rate by the total area and the interval to know how many pixels we 
      # are expected to disturb in the given interval for this study area, given the rate applied
      expectedDistPixels <- Rate*totAreaPix*dParOri[["disturbanceInterval"]]
      # newDisturbance <- oldDisturbance + (disturbanceRate * oldDisturbance) 
      if (expectedDistPixels < 1){
        # If the number of expected pixels for that given year is smaller than 1, then we need to
        # apply a probability of that disturbance happening. 
        message(crayon::red("The number of expected disturbed pixels for step ", 
                            currentTime, " is less than 1. The function will apply a probability ",
                            "of this disturbance happening, with the given size expected."))
        Size <- round(eval(parse(text = dParOri[["disturbanceSize"]])), 0)
        SizePix <- sqrt(Size)/res(rasterToMatch)[1]
        probHappening <- expectedDistPixels/SizePix
        totPixToChoose <- if (rbinom(1, 1, probHappening) == 1) SizePix else 0
      } else {
        totSizeChosen <- 0
        chosenDistrib <- numeric()
        IT <- 1
        while (round(totSizeChosen, 0) < round(expectedDistPixels, 0)){
          if (IT %% 1000 == 0){
            print(paste0("Calculating total disturbance size for ", 
                         ORIGIN, " (iteration ", IT, ", ", 
                         round(100*(round(totSizeChosen, 0)/round(expectedDistPixels, 0)), 2),"% achieved)"))
          }
          Size <- round(eval(parse(text = dParOri[["disturbanceSize"]])), 0)
          SizePix <- round(sqrt(Size)/res(rasterToMatch)[1], 0) # Convert the original res to pixels
          totSizeChosen <- totSizeChosen + SizePix
          chosenDistrib <- c(chosenDistrib, SizePix)
          IT <- IT + 1
        }
        message(paste0("Disturbance size successfully calculated for ", ORIGIN))
        totPixToChoose <- chosenDistrib
        totPixToChoose <- totPixToChoose[totPixToChoose > 0] # Need to clean up as sometimes we
      }
      # 4. Select which pixels are the most suitable and then get the number of neighbors based on 
      # the needed sizes
      # Can't forget to take pixels that already have disturbance out of the potential availability!
      # Otherwise it might choose pixels to disturb in existing disturbed ones
      # First, mask the current disturbances on the potential layer so we don't choose them again
      if (is.null(currentDisturbanceLayer)){
        message("currentDisturbanceLayer is NULL. This is likely the first year of the simulation")
        # Current layer is null, which indicates this might be the first year of the simulation. I
        # will need to create it then, based on the existing disturbances.
        # Get all current disturbances
        unDL <- unlist(disturbanceList)
        allCurrentDistLayNames <- names(unDL)[!grepl(x = names(unDL), pattern = "potential")]
        currDist <- unDL[names(unDL) %in% allCurrentDistLayNames]
        # Combine all layers, both polygons and lines, separtely
        linesLays <- currDist[sapply(currDist, geomtype) == "lines"]
        polyLays <- currDist[sapply(currDist, geomtype) == "polygons"]
        # Buffer lines at to convert them to polys to the data resolution of rasterToMatch 
        # This is not the ideal resolution, but if we are NOT deriving the buffered disturbancces 
        # from this layer but just using it to avoid choosing pixels already chosen, then it should 
        # be ok at 250m.  
        # NOTE: We devide the res by 2 so the total around is exactly the resolution of a pixel,
        # as the buffer uses the width on all sides. However, we need to ensure diagonals are also captured. 
        # So we need actually the half of the diagonal, not side. ~~~ Probably not needed... 
        # (sqrt(2)*res(rasterToMatch)[1])/2
        linesLays <- unlist(lapply(linesLays, function(X){
          terra::buffer(x = X, width = res(rasterToMatch)[1]/2)
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
        print("currentDisturbanceLayer is NOT null, second year ")
        browser()
        # Create allLaysRas
        
      }
      # Second, bring the raster with masked disturbances into a table to choose the pixels
      potLayF[allLaysRas == 1] <- NA
      
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
        allPixAdj <- as.numeric(sapply(whichAdj, function(totAdj){
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
        }))
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
      return(newDistLay)
    })
    names(updatedL) <- whichOrigin
    return(updatedL)
  })
  names(Generated) <- whichSector
  
  # THIRD: CONNECTING (Use generated and table)
  #############################################
  message("Connecting disturbances...")
  dPar <- disturbanceParameters[disturbanceType  %in% "Connecting", ]
  whichSector <- dPar[["dataName"]]
  Connected <- lapply(whichSector, function(Sector) {
    whichOrigin <- dPar[dataName == Sector, disturbanceOrigin]
    updatedL <- lapply(whichOrigin, function(ORIGIN) {
      Lay <- Generated[[Sector]][[ORIGIN]]
      #   whichOrigin2 <- dPar[dataName == Sector & disturbanceOrigin == ORIGIN, disturbanceOrigin]
      # if (length(whichOrigin) != 1) {
      #   print("Repeated origins and sector. Debug")
      #   browser()
      # } # Bug Catch
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
      # 1. Get the info on which layer to connect to which layer
      dParOri <- dPar[dataName == Sector & disturbanceOrigin == ORIGIN,]
      # 2. Get the origin layer (from Generated)
      oriLay <- Generated[[Sector]][[dParOri[["disturbanceOrigin"]]]]
      # 3. Get the end layer (from )
      endLay <- disturbanceList[[Sector]][[dParOri[["disturbanceEnd"]]]]
      # 4. Merge the layers and then connect the features
      # 4.1. Convert the raster to poly
      oriLaySF <- sf::st_as_sf(stars::st_as_stars(oriLay), as_points = FALSE, merge = TRUE)
      connected <- nngeo::st_connect(oriLaySF, sf::st_as_sf(endLay), ids = NULL, progress = TRUE)
      RES <- dParOri[["resolutionVector"]]
      if (any(is.na(RES),
              is.null(RES))){
        message(crayon::red(paste0("resolutionVector in disturbanceParameters",
                                   " table was not supplied for a lines file vector. Will ",
                                   "default to 15m. If this is not correct, please provide the ",
                                   "resolution used for developing the layer ", ORIGIN,
                                   " (", Sector,")")))
        RES <- 15
      }
      connectedLay <- terra::buffer(x = vect(connected), width = RES)
      return(connectedLay)
    })
    names(updatedL) <- whichOrigin
    return(updatedL)
  })
  names(Connected) <- whichSector
  
  # AT LAST: Update all layers so we can replace them in the list in the next function
  #############################################
  print("Update all layers so we can replace them in the list in the next function")
  browser()
  
  ### RETURN!
  # list(individuaLayers, currentDisturbanceLayer, forestry Potential layer -- which has been updated!! It needs to come back to this function)
}