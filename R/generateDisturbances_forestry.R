generateDisturbances_forestry <- function(disturbanceParameters,
                                 disturbanceList,
                                 fires,
                                 rasterToMatch,
                                 currentTime,
                                 currentDisturbanceLayer){
  
  # Forestry only has Generating and connecting

  # SECOND: Generating
  #############################################
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
      
      if (ORIGIN == "cutblocks"){
        # If we are dealing with cutblocks, we need to create a specific potential for the layers
        # based on the tree ages, crown closure and height at a certain age.
        # This is not exactly generic as there is a specific threshold.
        
        # TODO In future versions, I can pass a parameter saying how the age is classified
        # For now: 'Productive' forest stands were defined as those with a minimum Site Index (SI_50) 
        # value of 8 (which means stands have a height of at least 8 m at 50 years) and a minimum 
        # crown closure (CROWNCL) of 30%. (James Hodson GNWT's Definition). Also, trees need to be 
        # at least 50 years (based on ORIGIN) old.
        
        # On the original data, there is an error in one of the features (0.0009%): ORIGIN is defined as 1067
        # It is likely 1967. I will correct
        potLaySF[potLaySF$ORIGIN %in% c(1067), "ORIGIN"] <- 1967
        # There are also some features (0.28%) that have ORIGIN as 0 and no info on SI_50. So I will 
        # exclude these, as I have no way to know if they are productive or not.
        potLaySF <- potLaySF[potLaySF$ORIGIN != 0, ]
        
        # Applying James' rules (which I believe have been already applied!)
        potLaySF <- potLaySF[potLaySF$SI_50 >= 8 & potLaySF$CROWNCL >= 30, ]

        ######################################################################## TESTING!!
        fires <- raster::raster("rstCurrentBurn_2025_year2025.tif")
        ######################################################################## TESTING!!
        
        if (!is.null(fires)){
          
          # Also for forestry, if we are dealing with cutblocks, we need to update 
          # the layer to remove the fires. It is not a problem for any other disturbance, but with 
          # fire, there are no trees and therefore, no forestry.
          
          # As this layer was created in 2011 and has all the ages until 2010, I can use it
          # to update the age of the trees at each step, together with the fire layer. I just need to 
          # update the layer and feedback it to the function so I can keep track of the following years.
          
          # 1. Layover the fire layer on top of the potential layer. 
          browser() # HERE!!!!!!!!!!!!!!!!!!!!!
          fires
          # 2. Update the years of trees (ORIGIN) there based on the current fires
          # 3. Think if the rest of the code is fine!
          
        } else {
          # If fires is NULL, we are not passing data/simulated ones. The layer won't update, it will be
          # as there have been no fires.
          
          message(crayon::red(paste0("No fires were passed to the function. The layer won't be updated for ",
                                     "fires; it will be as there have been no fires during the simulation.",
                                     "This will, however, only affect forestry.")))
        }
      }
      
      # 3. Check how many pixels we need to choose, and in which sizes based on the table
      Rate <- dParOri[["disturbanceRate"]]
      totAreaPix <- sum(rasterToMatch[], na.rm = TRUE)
      # Here we multiply the rate by the total area and the interval to know how many pixels we 
      # are expected to disturb in the given interval for this study area, given the rate applied
      expectedDistPixels <- Rate*totAreaPix*dParOri[["disturbanceInterval"]] 
      
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
        totPixToChoose <- chosenDistrib
        totPixToChoose <- totPixToChoose[totPixToChoose > 0] # Need to clean up as sometimes we
      }
      # 4. Select which pixels are the most suitable and then get the number of neighbors based on 
      # the needed sizes
      # Can't forget to take pixels that already have disturbance out of the potential availability!
      # Otherwise it might choose pixels to disturb in existing disturbed ones
      # First, mask the current disturbances on the potential layer so we don't choose them again
      if (is.null(currentDisturbanceLayer)){
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
        # as the buffer uses the width on all sides.
        linesLays <- unlist(lapply(linesLays, function(X){
          terra::buffer(x = X, width = (res(rasterToMatch)[1])/2)
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
  dPar <- disturbanceParameters[disturbanceType  %in% "Connecting", ]
  whichSector <- dPar[["dataName"]]
  Connected <- lapply(whichSector, function(Sector) {
    whichOrigin <- dPar[dataName == Sector, disturbanceOrigin]
    if (length(whichOrigin) != 1) {
      print("Repeated origins and sector. Debug")
      browser()
    } # Bug Catch
    updatedL <- lapply(whichOrigin, function(ORIGIN) {
      Lay <- Generated[[Sector]][[ORIGIN]]
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
      browser()
      oriLaySF <- sf::st_as_sf(stars::st_as_stars(oriLay), as_points = FALSE, merge = TRUE)
      connected <- nngeo::st_connect(oriLaySF, sf::st_as_sf(endLay), ids = NULL, progress = TRUE)
      RES <- dParOri[["resolutionVector"]]
      if (any(is.na(RED),
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