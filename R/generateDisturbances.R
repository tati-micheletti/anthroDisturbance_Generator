generateDisturbances <- function(disturbanceParameters,
                                 disturbanceList,
                                 rasterToMatch,
                                 studyArea,
                                 fires,
                                 currentTime,
                                 growthStepGenerating,
                                 growthStepEnlargingPolys,
                                 growthStepEnlargingLines,
                                 currentDisturbanceLayer,
                                 connectingBlockSize,
                                 disturbanceRateRelatesToBufferedArea,
                                 outputsFolder,
                                 runName,
                                 checkDisturbancesForBuffer = FALSE){

  # Extracting layers from previous ones
  # Total study area
  rasterToMatchR <- raster::raster(rasterToMatch)
  growthStepEnlargingLines <- if (is.null(growthStepEnlargingLines) || is.na(growthStepEnlargingLines)) 1 else growthStepEnlargingLines

  safeMergeSpatVectors <- function(items,
                                   context = "spatial vectors",
                                   bufferWidth = NULL,
                                   targetCRS = NULL) {
    if (is.null(items)) return(NULL)
    if (!is.list(items)) items <- list(items)
    items <- Filter(function(x) {
      inherits(x, "SpatVector") &&
        tryCatch(nrow(x) > 0L, error = function(...) FALSE)
    }, items)
    if (!length(items)) {
      if (isTRUE(getOption("run_scenario.debug", FALSE))) {
        message("[generateDisturbances] no ", context, " to merge; returning NULL")
      }
      return(NULL)
    }
    baseWidth <- bufferWidth
    if (is.null(baseWidth) || is.na(baseWidth) || baseWidth <= 0) {
      baseWidth <- tryCatch(res(rasterToMatchR)[1], error = function(...) NA_real_)
    }
    if (is.na(baseWidth) || baseWidth <= 0) baseWidth <- 15
    convertToPolygon <- function(vec) {
      if (is.null(vec) || !inherits(vec, "SpatVector")) return(vec)
      geomType <- tryCatch(terra::geomtype(vec), error = function(...) NA_character_)
      if (all(!is.na(geomType) & geomType == "polygons")) return(vec)
      width <- baseWidth
      if (any(geomType == "points", na.rm = TRUE)) {
        width <- width/2
      }
      if (isTRUE(getOption("run_scenario.debug", FALSE))) {
        message("[generateDisturbances] buffering vector (geom=", paste(geomType, collapse = "/"),
                ") by ", signif(width, 4), " prior to merge in ", context)
      }
      vec <- terra::buffer(vec, width = width)
      terra::aggregate(vec, dissolve = TRUE)
    }
    items <- lapply(items, function(obj) {
      geomType <- tryCatch(terra::geomtype(obj), error = function(...) NA_character_)
      if (all(!is.na(geomType) & geomType == "polygons")) return(obj)
      convertToPolygon(obj)
    })
    if (length(items) == 1L) {
      merged <- items[[1]]
    } else {
      merged <- Reduce(function(acc, x) {
      if (is.null(acc) || !inherits(acc, "SpatVector") ||
          tryCatch(nrow(acc) == 0L, error = function(...) TRUE)) {
        return(x)
      }
      if (is.null(x) || !inherits(x, "SpatVector") ||
          tryCatch(nrow(x) == 0L, error = function(...) TRUE)) {
        return(acc)
      }
      acc <- convertToPolygon(acc)
      x <- convertToPolygon(x)
      tryCatch(
        do.call(rbind, list(acc, x)),
        error = function(err) {
          geomAcc <- tryCatch(terra::geomtype(acc), error = function(...) NA_character_)
          geomX <- tryCatch(terra::geomtype(x), error = function(...) NA_character_)
          warning(paste0("[generateDisturbances] failed to merge vectors (geom acc=",
                         paste(geomAcc, collapse = "/"), ", geom x=",
                         paste(geomX, collapse = "/"), ") for ", context,
                         ": ", conditionMessage(err), ". Returning accumulator."),
                  call. = FALSE, immediate. = TRUE)
          acc
        }
      )
    }, items[-1], init = items[[1]])
    }
  if (!is.null(targetCRS) && inherits(merged, "SpatVector")) {
    merged <- tryCatch(terra::project(merged, targetCRS), error = function(e) merged)
  }
  merged
}
  expanse_vec <- function(x, unit = "m", transform = FALSE) {
    vals <- tryCatch(terra::expanse(x, unit = unit, transform = transform),
                     error = function(...) NA_real_)
    if (is.data.frame(vals)) {
      vals <- if ("area" %in% names(vals)) vals[["area"]] else unlist(vals, use.names = FALSE)
    } else if (is.list(vals)) {
      vals <- unlist(vals, use.names = FALSE)
    }
    suppressWarnings(as.numeric(vals))
  }
  sumExpanse <- function(x, unit = "m", transform = FALSE) {
    vals <- expanse_vec(x, unit = unit, transform = transform)
    if (!length(vals)) return(0)
    sum(vals, na.rm = TRUE)
  }
  
  if (!is(studyArea, "SpatVector"))
    studyArea <- terra::vect(studyArea)
    
  studyArea <- terra::project(x = studyArea, y = terra::crs(rasterToMatch))
  uniStudyArea <- terra::aggregate(studyArea)
  totalstudyAreaVAreaSqKm <- sumExpanse(uniStudyArea, unit = "km", transform = FALSE)
  totalstudyAreaVAreaSqm <- sumExpanse(uniStudyArea, unit = "m", transform = FALSE)
  totNPix <- sum(rasterToMatch[], na.rm = TRUE)
  calculatedPixelSizem2 <- totalstudyAreaVAreaSqm/totNPix # in m2
  
  # First, mask the current disturbances on the potential layer so we don't choose them again
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
    # Buffer lines at the resolution of the rasterToMatch to convert them to polys.
    # This is not the ideal resolution, but as we are NOT deriving the buffered disturbances
    # from this layer but just using it to avoid choosing pixels already chosen, then it should
    # be ok at 250m.
    # NOTE: We devide the res by 2 so the total around is exactly the resolution of a pixel,
    # as the buffer uses the width on all sides.

    if (length(linesLays)) {
      linesLays <- lapply(linesLays, function(X){
        terra::buffer(x = X, width = res(rasterToMatch)[1]/2)
        # Here we need to buffer as this layer will be used to avoid specific pixels (i.e., already 
        # disturbed!) to be selected for new disturbances.
      })
    }
    if (length(linesLays)) {
      linesLays <- Filter(function(x) {
        inherits(x, "SpatVector") &&
          tryCatch(nrow(x) > 0L, error = function(...) FALSE)
      }, linesLays)
    }
    if (length(polyLays)) {
      polyLays <- Filter(function(x) {
        inherits(x, "SpatVector") &&
          tryCatch(nrow(x) > 0L, error = function(...) FALSE)
      }, polyLays)
    }
    linesAndPolys <- c(polyLays, linesLays)
    if (!length(linesAndPolys)) {
      message("No current disturbances found in study area; skipping generation for this year.")
      return(currentDisturbanceLayer)
    }
    if (isTRUE(getOption("run_scenario.debug", FALSE))) {
      message("[generateDisturbances] merging linesAndPolys length=", length(linesAndPolys))
      message("[generateDisturbances] linesAndPolys classes=", paste(vapply(linesAndPolys, function(x) {
        if (is.null(x)) return("NULL")
        paste(class(x), collapse = ":")
      }, character(1)), collapse=" | "))
    }
    mergeWidth <- tryCatch(res(rasterToMatchR)[1], error = function(...) NA_real_)
    allLays <- safeMergeSpatVectors(linesAndPolys,
                                    context = "current disturbance layers",
                                    bufferWidth = mergeWidth,
                                    targetCRS = terra::crs(rasterToMatchR))
  } else {
    if (is(currentDisturbanceLayer, "RasterLayer"))
      allLaysRas <- terra::rast(currentDisturbanceLayer) else
        if (is(currentDisturbanceLayer, "SpatRaster")) 
          allLaysRas <- currentDisturbanceLayer else {
            curDistVcs <- unlist(currentDisturbanceLayer)
            RES <- unique(disturbanceParameters[["resolutionVector"]])
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
            curDistVcsAll <- Filter(Negate(is.null), curDistVcsAll)
            mergeWidth <- tryCatch(res(rasterToMatchR)[1], error = function(...) NA_real_)
            allLays <- safeMergeSpatVectors(curDistVcsAll,
                                            context = "existing disturbance vectors",
                                            bufferWidth = mergeWidth,
                                            targetCRS = terra::crs(rasterToMatchR))
          }
    # Create allLaysRas, which is the disturbance raster with all disturbances rasterized, with
    # value == 1 for the disturbed pixels. Non-disturbed pixels are NA
  }
  # Now fasterize all polys
  allLaysRas <- NULL
  if (!is.null(allLays) && inherits(allLays, "SpatVector") &&
      tryCatch(nrow(allLays) > 0L, error = function(...) FALSE)) {
    allLaysSF <- sf::st_as_sf(x = allLays)
    if (nrow(allLaysSF) > 0L) {
      geomTypes <- sf::st_geometry_type(allLaysSF, by_geometry = TRUE)
      polyIdx <- geomTypes %in% c("POLYGON", "MULTIPOLYGON")
      if (!any(polyIdx, na.rm = TRUE)) {
        if (isTRUE(getOption("run_scenario.debug", FALSE))) {
          message("[generateDisturbances] no polygon geometries available for fasterize; skipping raster conversion")
        }
      } else {
        polySF <- allLaysSF[polyIdx, , drop = FALSE]
        if (nrow(polySF) > 0L) {
          emptyIdx <- sf::st_is_empty(polySF)
          if (any(emptyIdx, na.rm = TRUE)) {
            if (isTRUE(getOption("run_scenario.debug", FALSE))) {
              message("[generateDisturbances] dropping ", sum(emptyIdx, na.rm = TRUE),
                      " empty geometries before fasterize")
            }
            polySF <- polySF[!emptyIdx, , drop = FALSE]
          }
        }
        if (nrow(polySF) > 0L) {
          validIdx <- sf::st_is_valid(polySF)
          if (any(!validIdx, na.rm = TRUE)) {
            if (requireNamespace("lwgeom", quietly = TRUE)) {
              polySF <- lwgeom::st_make_valid(polySF)
            } else {
              warnMsg <- "[generateDisturbances] invalid geometries detected but 'lwgeom' not available; dropping invalid features before fasterize."
              warning(warnMsg, call. = FALSE, immediate. = TRUE)
              polySF <- polySF[validIdx, , drop = FALSE]
            }
          }
        }
        if (nrow(polySF) > 0L) {
          polySF$Polys <- 1
          polysField <- "Polys"
          message("Rasterizing disturbances...")
          rasterToMatchR <- rasterToMatchR
          polyVec <- terra::vect(polySF)
          allLaysRas <- terra::rasterize(polyVec,
                                         terra::rast(rasterToMatchR), 
                                         field = polysField,
                                         background = NA_real_,
                                         touches = TRUE)
        }
      }
    }
  }
  

  # FIRST: Enlarging
  #############################################
  message(crayon::white("Enlarging disturbances..."))
  dPar <- disturbanceParameters[disturbanceType  %in% "Enlarging", ]
  whichSector <- dPar[["dataName"]]
  Enlarged <- lapply(whichSector, function(Sector) {
    whichOriginRaw <- dPar[dataName == Sector, disturbanceOrigin]
    whichOrigin <- unique(whichOriginRaw)
    if (!length(whichOrigin)) {
      warning("[generateDisturbances] No disturbanceOrigin rows found for ", Sector,
              " in Enlarging; skipping.", call. = FALSE)
      return(NULL)
    }
    if (length(whichOriginRaw) != length(whichOrigin)) {
      message(crayon::yellow(paste0(
        "[generateDisturbances] Duplicate origins for ", Sector, "; using unique values: ",
        paste(whichOrigin, collapse = ",")
      )))
    }
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
        currArea <- sumExpanse(Lay, unit = "m", transform = FALSE)
      } else { # If the disturbanceRate Relates To Buffered Area (i.e., caribou-wise)
        # It doesn't matter if poly or lines or points, if buffered, is buffered to 500m
          LayBuff <- terra::buffer(Lay, width = 500) # Need to aggregate to avoid double counting!
          LayBuff <- terra::aggregate(x = LayBuff, dissolve = TRUE)
          currArea <- sumExpanse(LayBuff, unit = "m", transform = FALSE) # NEEDS TO BE METERS. USED BELOW
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
            futureArea <- sumExpanse(futureAreaAgg, unit = "m", transform = FALSE)
          } else {
            futureArea <- sumExpanse(LayUpdated, unit = "m", transform = FALSE)
          }
          totalDisturbedAreaAchieved <- futureArea
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
  
  # SECOND: Generating
  #############################################
  message(crayon::white("Generating disturbances..."))
  dPar <- disturbanceParameters[disturbanceType  %in% "Generating", ]
  whichSector <- dPar[["dataName"]]
  Generated <- lapply(whichSector, function(Sector) {
    whichOriginRaw <- dPar[dataName == Sector, disturbanceOrigin]
    whichOrigin <- unique(whichOriginRaw)
    if (!length(whichOrigin)) {
      warning("[generateDisturbances] No disturbanceOrigin rows found for ", Sector,
              " in Generating; skipping.", call. = FALSE)
      return(NULL)
    }
    if (length(whichOriginRaw) != length(whichOrigin)) {
      message(crayon::yellow(paste0(
        "[generateDisturbances] Duplicate origins for ", Sector, "; using unique values: ",
        paste(whichOrigin, collapse = ",")
      )))
    }
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
      if (length(potLay) == 0){# In case the cropped area doens't have anything
        message(paste0("The potential area for ", Sector, " class ", ORIGIN, " is NULL.",
                       " Likely cropped out from studyArea. Returning NULL."))
        return(NULL)
      } 
      # Previous pixel strategy for Generating Disturbances
      # 2. Fasterize it
      potLaySF <- sf::st_as_sf(x = potLay)
      potField <- dParOri[["potentialField"]]
      
      if (any(is.na(potField), 
              potField == "",
              is.null(potField))){ 
        # If NA, it doesn't matter, but need to create a 
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
                       "occurred fires and currently productive forest for year ", currentTime))
        
        # First: Select only productive forests
        # potLaySF <- potLaySF[potLaySF$ORIGIN < (currentTime - 50), ] # For some weird reason, this doesn't work anymore... sigh.
        potLaySF <- subset(potLaySF, potLaySF$ORIGIN < (currentTime - 50))
        potLayF <- fasterize::fasterize(sf = st_collection_extract(potLaySF, "POLYGON"),
                                        raster = rasterToMatchR, field = potField)
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
                                        raster = rasterToMatchR, field = potField)
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
        newDistLay <- copy(rasterToMatch)
        names(newDistLay) <- ORIGIN
        newDistLay[newDistLay == 1] <- 0
        return(newDistLay)
      }
      # Ideally, I would do pixels here. However, because a lot of the disturbances are smaller than 
      # 1 pixel (6.25 ha), doing with vectors is better.
      
      expectedDistPixels <- Rate*totNPix*dParOri[["disturbanceInterval"]] # The disturbance expected here 
      # is already considering the buffer if disturbanceRateRelatesToBufferedArea, or the not buffered if
      # the parameter is FALSE
      
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
        currArea <- sumExpanse(LayBuff, unit = "m", transform = FALSE) # NEEDS TO BE METERS. USED BELOW
        if (length(currArea) == 0){
          message(paste0("No disturbance for ", Sector, 
                         " -- ", ORIGIN, ": currentArea = 0 Km2"))
        } else {
          message(paste0("Buffered (500m) area for ", Sector, 
                         " -- ", ORIGIN, ": ", round(currArea/1000000, 2), " Km2"))
          message(paste0("Percentage of current area: ", round(100*((currArea/1000000)/totalstudyAreaVAreaSqKm), 3), "%."))
        }
      }
      
      if (expectedDistPixels < 1){
        # If the number of expected pixels for that given year is smaller than 1, then we need to
        # apply a probability of that disturbance happening. 
        message(crayon::red("The number of expected disturbed pixels for ", ORIGIN," for step ", 
                            currentTime, " is less than 1. The function will apply a probability ",
                            "of this disturbance happening, with the given expected size."))
        expectedDistPixelsProb <- if (disturbanceRateRelatesToBufferedArea) (1/8)*expectedDistPixels else 
          expectedDistPixels # As the disturbance needs to be buffered, we need to reduce the probability by 1/8 (total buffer of 1 pixel)
        totPixToChoose <- rbinom(1, 1, expectedDistPixelsProb)
        if (totPixToChoose == 0){
          message(paste0("Rate of disturbance for ", ORIGIN, " is very small and the probability of ",
                         "this disturbance happening returned 0.",
                         " Returning layer without disturbances."))
          newDistLay <- copy(rasterToMatch)
          names(newDistLay) <- ORIGIN
          newDistLay[newDistLay == 1] <- 0
          return(newDistLay)
        }
        nPixChosenTotal <- totPixToChoose
      } else {
        nPixChosenTotal <- 0
        chosenDistrib <- numeric()
        IT <- 1
        while (round(nPixChosenTotal, 0) < round(expectedDistPixels, 0)){
          if (IT %in% 1:10){
            message(paste0("Calculating total generated disturbance size for ", Sector, " for ", 
                           ORIGIN, " (Year ", currentTime,"; iteration ", IT, ", ", 
                           round(100*(round(nPixChosenTotal, 0)/round(expectedDistPixels, 0)), 2),"% achieved)"))
          }
          if (IT %% 10 == 0){
            message(paste0("Calculating total generated disturbance size for ", Sector, " for ",  
                           ORIGIN, " (Year ", currentTime,"; iteration ", IT, ", ", 
                           round(100*(round(nPixChosenTotal, 0)/round(expectedDistPixels, 0)), 2),"% achieved)"))
          }
          # This is the size in meter square that each disturbance should have
          expr <- dParOri[["disturbanceSize"]]
          # Ensure truncated-normal sampling uses fully-qualified msm::rtnorm
          expr <- gsub("(?<!msm::)rtnorm\\(", "msm::rtnorm(", expr, perl = TRUE)
          expr <- gsub("msm::msm::", "msm::", expr, fixed = TRUE)
          Size <- round(eval(parse(text = expr)), 0)
          # Here I generate one disturbance at a time. 
          # This is the size as the number of pixels to be chosen IN m2!
          # totNPix is the number of pixels that correspont to totalstudyAreaVAreaSqKm
          # As the calculation based on pixel resolution is not good, we can try a better approximation
          # by using the study area shapefile calculated with terra. It is not the best as southern pixels
          # are in reality smaller than the northern ones, but it is a good approximation for now: 
          calculatedPixelSizeLength <- sqrt(calculatedPixelSizem2)
          # sqrt(Size) needs to be at least the same as calculatedPixelSizeLength (or better, 
          # Size at least calculatedPixelSizeLength^2 == at least 1 pixel), otherwise we have 0 pixels selected
          nPixChosen <- round(sqrt(Size)/calculatedPixelSizeLength, 0)
          if (nPixChosen == 0)
            nPixChosen <- 1
          # TODO Improve/simplify this. I might already here have either 1 or more.
          if (Size < calculatedPixelSizem2){
            # This means that, if I have a 0, I need to force it to 1. This is because even if we had 
            # everything as vector, we would still have one pixel when rasterizing, even if the 
            # disturbance is really small.
            if (disturbanceRateRelatesToBufferedArea){ # Size smaller than a pixel and buffered!
              # browser() # Implement the speeding up here
              # if (expectedDistPixels > (nPixChosen * growthStepGenerating)) # Trying to make things run faster when they can
              #   nPixChosen <- nPixChosen * growthStepGenerating
              
              # 1. Create a raster with the same resolution than RTM
              nrowscols <- 100*nPixChosen
              # Make sure the raster is more or less in the same region as we are,
              # so use the ymax and xmax of rasterToMatch to guide the extent
              tmpRas <- rast(ncols = nrowscols,
                             nrows = nrowscols,
                             res = terra::res(rasterToMatch),
                             extent = c(xmax(rasterToMatch)-(res(rasterToMatch)*nrowscols)[1], 
                                        xmax(rasterToMatch),
                                        ymax(rasterToMatch)-(res(rasterToMatch)*nrowscols)[1], 
                                        ymax(rasterToMatch)), # xmin, xmax, ymin, ymax
                             crs = terra::crs(rasterToMatch),
                             vals = NA)
              # 2. Get the pixToChoose in the middle of the raster
              middlePoint <- getRasterMiddlePoint(tmpRas)
              # This means that we have to randomly choose the nearest 8 pixel neighbors
              adjCels <- c(middlePoint,
                           sample(x = terra::adjacent(tmpRas, cells = middlePoint, 
                                                      directions = 8, include = FALSE), 
                                  size = (nPixChosen - 1),
                                  replace = FALSE))
              # Place the chosen pixels in the raster
              tmpRas[adjCels] <- 1
              # 3. Buffer it using terra buffer by 500m (background = NA)
              bufferedTmpRas <- terra::buffer(tmpRas,
                                              width = 500, background = -1)
              bufferedTmpRas[bufferedTmpRas == -1] <- NA
              # 4. Get how many pixels became 1's, this is my 
              nPixChosenBuff <- sum(bufferedTmpRas[], na.rm = TRUE)
              #TODO Implement growthStepGenerating when buffering...
              # if (expectedDistPixels > (nPixChosenBuff * growthStepGenerating)) # Trying to make things run faster when they can
              #   nPixChosen <- nPixChosen * growthStepGenerating  
              nPixChosenTotal <- nPixChosenTotal + nPixChosen + nPixChosenBuff
              
            } else { # Size smaller than a pixel, and NOT buffered
              # Similarly, Size needs to be at least 1 pixel for non-disturbanceRateRelatesToBufferedArea.
              if (expectedDistPixels > (nPixChosen * growthStepGenerating)) # Trying to make things run faster when they can
                nPixChosen <- nPixChosen * growthStepGenerating
              nPixChosenTotal <- nPixChosenTotal + nPixChosen
            }
          } else { 
           # Here I am getting the size of the disturbance if Size is bigger than a pixel
          if (disturbanceRateRelatesToBufferedArea){ # Size bigger than a pixel and buffered
            # browser() # Implement the speeding up here
            # if (expectedDistPixels > (nPixChosen * growthStepGenerating)) # Trying to make things run faster when they can
            #   nPixChosen <- nPixChosen * growthStepGenerating
               
            # 1. Create a raster with the same resolution than RTM
            nrowscols <- 100*nPixChosen
            # Make sure the raster is more or less in the same region as we are,
            # so use the ymax and xmax of rasterToMatch to guide the extent
            tmpRas <- rast(ncols = nrowscols,
                           nrows = nrowscols,
                           res = res(rasterToMatch),
                           extent = c(xmax(rasterToMatch)-(res(rasterToMatch)*nrowscols)[1], 
                                      xmax(rasterToMatch),
                                      ymax(rasterToMatch)-(res(rasterToMatch)*nrowscols)[1], 
                                      ymax(rasterToMatch)), # xmin, xmax, ymin, ymax
                           crs = terra::crs(rasterToMatch),
                           vals = NA)
            # 2. Get the pixToChoose in the middle of the raster
            middlePoint <- getRasterMiddlePoint(tmpRas)
            if (nPixChosen < 10){
             # This means that we have to randomly choose the nearest 8 pixel neighbors
              adjCels <- c(middlePoint,
                           sample(x = terra::adjacent(tmpRas, cells = middlePoint, 
                                         directions = 8, include = FALSE), 
                                  size = (nPixChosen - 1),
                                  replace = FALSE))
            } else {
              chosenPix <- middlePoint
              while (length(chosenPix) < nPixChosen){
                adjCels <- terra::adjacent(tmpRas, cells = chosenPix,
                                           directions = 8, 
                                           include = FALSE)
                chosenPix <- c(chosenPix, adjCels)
              }
              adjCels <- chosenPix[1:nPixChosen]
            }
            # Place the chosen pixels in the raster
            tmpRas[adjCels] <- 1
            # 3. Buffer it using terra buffer by 500m (background = NA)
            bufferedTmpRas <- terra::buffer(tmpRas,
                                            width = 500, background = -1)
            bufferedTmpRas[bufferedTmpRas == -1] <- NA
            # 4. Get how many pixels became 1's, this is my 
            nPixChosenBuff <- sum(bufferedTmpRas[], na.rm = TRUE)
            #TODO Implement growthStepGenerating when buffering...
            # if (expectedDistPixels > (nPixChosenBuff * growthStepGenerating)) # Trying to make things run faster when they can
            #   nPixChosen <- nPixChosen * growthStepGenerating  
            nPixChosenTotal <- nPixChosenTotal + nPixChosen + nPixChosenBuff
            } else {# Size bigger than a pixel and NOT buffered
            if (expectedDistPixels > (nPixChosen * growthStepGenerating)) # Trying to make things run faster when they can
              nPixChosen <- nPixChosen * growthStepGenerating
            nPixChosenTotal <- nPixChosenTotal + nPixChosen
          }
          }
          chosenDistrib <- c(chosenDistrib, nPixChosen)
          IT <- IT + 1
        }
        # TODO We could here implement code to choose which iteration had less 
        # difference (i.e., past or current)
        totPixToChoose <- chosenDistrib[chosenDistrib > 0] # Need to clean up as sometimes we
        # get zero pixels if the size of each disturbance is not too big 
        # [UPDATE: This shouldn't happen anymore with the changes in the code]
      }
      # Upper-bound requested pixels by available potential
      availablePix <- sum(!is.na(terra::values(potLayF)))
      if (is.finite(availablePix)) {
        expectedDistPixels <- min(expectedDistPixels, availablePix)
        totPixToChoose <- pmin(totPixToChoose, availablePix)
      }
      # 4. Select which pixels are the most suitable and then get the number of neighbors based on 
      # the needed sizes
      # Can't forget to take pixels that already have disturbance out of the potential availability!
      # Otherwise it might choose pixels to disturb in existing disturbed ones
      # Second, bring the raster with masked disturbances into a table to choose the pixels
      if (!is.null(allLaysRas)) {
        potLayF[allLaysRas[] == 1] <- NA
      } else if (isTRUE(getOption("run_scenario.debug", FALSE))) {
        message("[generateDisturbances] allLaysRas is NULL; skipping mask of existing disturbances")
      }
      potentialDT <- na.omit(data.table(pixelID = 1:ncell(potLayF),
                                vals = terra::values(potLayF)))
      # Third, subset the best pixels
      potentialDT <- potentialDT[is.finite(vals)]
      if (!nrow(potentialDT)) {
        message(crayon::red(paste0(
          "[generateDisturbances] No finite potential values for ", ORIGIN, " (", Sector,
          "); skipping disturbance generation."
        )))
        return(NULL)
      }
      potentialVals <- sort(unique(potentialDT[["vals"]]), decreasing = TRUE)
      bestPotential <- potentialDT[vals == potentialVals[1], pixelID]
      # Then, randomly, choose the pixels. Here need to pay attention to choose one per "development"
      # and only then use the adjacent to get to the correct sizes.
      # length(totPixToChoose) --> is the number of pixels minus the ones that don't need adjacent. 
      # The totPixToChoose is then the sizes each of these need to have with adjacent cells
      # and sum(totPixToChoose) is the total number of pixels where disturbance will be generated, which
      # should match very closely to expectedDistPixels
      allSelected <- FALSE
      if (sum(totPixToChoose) > length(bestPotential)){
        # This means that the best potential is not enough. So we need to increase to the 2 best
        # potentials (and keep expanding until either we have enough or exhaust all potential).
        counter <- 1
        while (sum(totPixToChoose) > length(bestPotential) &&
               counter < length(potentialVals)) {
          classesToInclude <- 1:min(counter + 1, length(potentialVals))
          # subset the best pixels
          bestPotential <- potentialDT[vals %in% potentialVals[classesToInclude],
                                       pixelID]
          counter <- counter + 1
        }
        if (sum(totPixToChoose) > length(bestPotential)) {
          message(crayon::red(paste0("All pixels available for ", ORIGIN, " have been selected. ",
                                     "This distubance will remain static until the end of the ",
                                     "simulation.")))
          pixNewDist <- bestPotential
          allSelected <- TRUE
        } else if (length(totPixToChoose) > length(bestPotential)) {
          message(crayon::red(paste0("Number of requested disturbances for ", ORIGIN,
                                     " exceeds available potential locations (", length(totPixToChoose),
                                     " > ", length(bestPotential), "). Selecting all available cells.")))
          pixNewDist <- bestPotential
          allSelected <- TRUE
        } else {
          pixNewDist <- sample(x = bestPotential,
                               size = length(totPixToChoose),
                               replace = FALSE)
        }
      } else {
        if (length(totPixToChoose) > length(bestPotential)) {
          message(crayon::red(paste0("Number of requested disturbances for ", ORIGIN,
                                     " exceeds available potential locations (", length(totPixToChoose),
                                     " > ", length(bestPotential), "). Selecting all available cells.")))
          pixNewDist <- bestPotential
          allSelected <- TRUE
        } else {
          pixNewDist <- sample(x = bestPotential,
                               size = length(totPixToChoose),
                               replace = FALSE)
        }
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
        print("Doing focal weight...")
        nbMatrices <- lapply(whichAdj, function(X){
          m <- focalWeight(x = rasterToMatchR, d = X * res(rasterToMatch)[1], 
                      type = c('circle'))
          m[m > 0] <- 1
          return(m)
        })
        names(nbMatrices) <- paste0("NB_", whichAdj)
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
          Adj <- as.data.table(raster::adjacent(x = rasterToMatchR, cells = pix,
                                  directions = nb, pairs = TRUE,
                                  include = FALSE, id = TRUE))
          # Now we sample the number of adjacents needed from the whole table, by ID (each pix that 
          # needs those adjacents). If adjN exceeds available neighbours, fall back to sampling
          # with replacement to keep execution stable.
          chosen <- Adj[, {
            n_to_sample <- min(.N, adjN)
            .SD[sample(.N, n_to_sample, replace = n_to_sample > .N)]
          }, by = "id"]
          return(c(unique(chosen[["from"]]), unique(chosen[["to"]])))
        })
        allPixAdj <- unique(as.numeric(unlist(allPixAdj)))
        whichPixelsChosen <- c(allPixAdj, pixDontNeedAdj)
      } else {
        whichPixelsChosen <- pixNewDist
      }
      # If we ended up with more pixels than requested, trim back to expectedDistPixels
      if (length(whichPixelsChosen) > expectedDistPixels) {
        whichPixelsChosen <- sample(whichPixelsChosen, expectedDistPixels)
      }
      # Now, we make a new disturbance layer. We need to keep track of the 
      # old one because we need to connect the new stuff. We can use RTM as a template.
      newDistLay <- copy(rasterToMatch)
      names(newDistLay) <- ORIGIN
      newDistLay[newDistLay == 1] <- 0
      newDistLay[whichPixelsChosen] <- 1
      if (disturbanceRateRelatesToBufferedArea){ 
        # Because of the projection problem described below,
        # this whole section is a moot point. Using instead what was previously calculated.
        # Because of where this is (NT), the buffering below results in errors of distance.
        # TODO To fix this, I would need to implement a final buffering (i.e., when we want to 
        # buffer all disturbances for caribou) system using moving windows somehow. This would remove 
        # the effects of earth projection distortions.
        
        # If disturbanceRateRelatesToBufferedArea, we need to recalculate the total pixels!
        # newDistLayBuff <- newDistLay
        # newDistLayBuff[newDistLayBuff == 0] <- NA
        # newDistLayBuff2 <- terra::buffer(newDistLayBuff, width = 500, background = -1)
        # newDistLayBuff2[newDistLayBuff2 == -1] <- NA
        # howManyPixelsChosenBuff <- sum(newDistLayBuff2[], na.rm = TRUE)
        # totPixBuff <- sum(length(whichPixelsChosen), howManyPixelsChosenBuff)
        totPixBuff <- nPixChosenTotal
        calcPerc <- (totPixBuff - expectedDistPixels)/expectedDistPixels
      } else {
        totPixBuff <- length(whichPixelsChosen)
        calcPerc <- (length(whichPixelsChosen) - expectedDistPixels)/expectedDistPixels
      }
      
      message(paste0("Percentage of disturbed future area after buffer: ",
                     round(100*((totPixBuff*calculatedPixelSizem2)/totalstudyAreaVAreaSqm), 3), "%."))
       
      cat(crayon::yellow(paste0("Difference between expected and achieved change for ",
                                crayon::red(Sector), " -- ", crayon::red(ORIGIN), ": ",
                                crayon::red(format(100*round(calcPerc, 4), 
                                                   scientific = FALSE), " % (ideal value = 0)."),
                                "\nDisturbance achieved: ", round((totPixBuff*calculatedPixelSizem2)/1000000,3),
                                " km2 -- Disturbance expected: ", round((expectedDistPixels*calculatedPixelSizem2)/1000000,3), " km2")))
      cat(paste0(Sector, 
                 " ", 
                 ORIGIN, 
                 " ",
                 currentTime,
                 " ",
                 format(100*round(calcPerc, 8), 
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
        if (!is(oriLay, "SpatRaster"))
          oriLay <- rast(oriLay) # Using terra is faster
        # 4. Merge the layers and then connect the features
        # 4.1. Convert the raster to poly
        # Need to convert all that is zero into NA otherwise 2 polygons are created (one for 
        # 0's and one for 1's):
        
        oriLay[oriLay != 1] <- NA
        oriLayVect <- terra::as.polygons(oriLay)
        if (NROW(oriLayVect) == 0){
          # If the disturbance doesn't exist, we need to return NULL
          return(NULL)
        }
        # Here I should go over each individual disturbance and connect it iteratively.
        # This might eliminate the problem with so many lines at the same time very close to 
        # each other. I probably need to:
        # 1. Make a for loop over the number of lines, 
        # 2. select just one feature at a time
        # oriLayVectSingle <- st_as_sf(oriLayVect) # For sf and st_connect. 
        # 3. use the st_connect or terra::nearest on it # Changed to terra's native nearest. Much faster!
        classEndLay <- na.omit(unique(endLay[["Class"]])) 
        blockSize <- suppressWarnings(as.numeric(connectingBlockSize))
        if (length(blockSize) != 1) blockSize <- NA_real_
        blockIsActive <- isTRUE(is.finite(blockSize) && blockSize > 0)
        if (blockIsActive && NROW(oriLayVect) > blockSize){
          message(paste0("connectingBlockSize is ", blockSize, " and connecting layer has ", 
                         NROW(oriLayVect), " rows. Applying blocking technique to speed up",
                         "disturbance generation type connecting. ",
                         " If too many lines are connecting from the same place, decrease the",
                         "connectingBlockSize or ",
                         "set it to NULL, which will improve final result, but needs considerable",
                         "more time to run."))
          blockFullList <- 1:NROW(oriLayVect)
          amountBlocks <- ceiling(NROW(oriLayVect)/blockSize)
          blockList <- list()
          for (i in 1:amountBlocks) {
            if (length(blockFullList) <= blockSize){
              blockList[[i]] <- blockFullList
            } else {
              blockList[[i]] <- sample(x = blockFullList, 
                                       size = blockSize, replace = FALSE)
              # Update the list
              blockFullList <- setdiff(blockFullList, blockList[[i]])
            }
          }
          names(blockList) <- paste0("block_", 1:amountBlocks)
          for (i in names(blockList)){
            whichVecs <- blockList[[i]]
            message(paste0("Connecting ", unique(names(oriLayVect))," for year ", currentTime,
                           ": ", i ," of ", length(blockList),
                           " (", 100*round(length(whichVecs)/NROW(oriLayVect),3),"%)"))
            connectedOne <- terra::nearest(oriLayVect[whichVecs, ], endLay, 
                                           pairs = FALSE, 
                                           centroids = TRUE, 
                                           lines = TRUE)
            # 4. Update the endLayer with the new connection
            if (!is(connectedOne, "SpatVector"))
              connectedOne <- terra::vect(connectedOne)
        connectWidth <- tryCatch(res(rasterToMatchR)[1], error = function(...) NA_real_)
        if (exists("connected", inherits = FALSE) && !is.null(connected)){
          connected <- safeMergeSpatVectors(list(connected, connectedOne),
                                            context = "connecting block merge",
                                            bufferWidth = connectWidth)
            } else {
              connected <- connectedOne
            }
            endLay <- safeMergeSpatVectors(list(endLay, connectedOne),
                                           context = "connecting end lay",
                                           bufferWidth = connectWidth)
          }          
        } else {
          # message(paste0("connectingBlockSize is ", 
          #                ifelse(is.null(connectingBlockSize), "NULL", connectingBlockSize),
          #                " and connecting layer has ", NROW(oriLayVect),
          #                " rows. NOT applying blocking technique for disturbance generation type ",
          #                "connecting for ", Sector," -- ", disturbanceOrigin, ". If taking too long, ",
          #                "please provide connectingBlockSize as a higher number than ", 
          #                ifelse(is.null(connectingBlockSize), "1", connectingBlockSize)))
          for (i in 1:NROW(oriLayVect)){
            if(i%%100==0)
              message(paste0("Connecting ", unique(names(oriLayVect))," for year ", currentTime,
                           ": ", i ," of ", NROW(oriLayVect),
                           " (", 100*round(i/NROW(oriLayVect),4),"%)"))
            connectedOne <- terra::nearest(oriLayVect[i, ], endLay, 
                                           pairs = FALSE, 
                                           centroids = TRUE, 
                                           lines = TRUE)
            # 4. Update the endLayer with the new connection
            if (!is(connectedOne, "SpatVector"))
              connectedOne <- terra::vect(connectedOne)
            connectWidth <- tryCatch(res(rasterToMatchR)[1], error = function(...) NA_real_)
            if (exists("connected", inherits = FALSE) && !is.null(connected)){
              connected <- safeMergeSpatVectors(list(connected, connectedOne),
                                                context = "connecting block merge",
                                                bufferWidth = connectWidth)
            } else {
              connected <- connectedOne
            }
            endLay <- safeMergeSpatVectors(list(endLay, connectedOne),
                                           context = "connecting end lay",
                                           bufferWidth = connectWidth)
          }
        }
        connected[["Class"]] <- classEndLay
        Lay <- list(connected)
        names(Lay) <- disturbanceEnd
        if (disturbanceRateRelatesToBufferedArea){
          if (length(disturbanceEnd) > 1){
            warning("[generateDisturbances] Multiple disturbanceEnd targets found (",
                    paste(disturbanceEnd, collapse = ","), "); using the first entry.",
                    call. = FALSE)
            disturbanceEnd <- disturbanceEnd[1]
            names(Lay) <- disturbanceEnd
          }
          LayBuff <- terra::buffer(Lay[[disturbanceEnd]], width = 500) # Need to aggregate to avoid double counting!
          LayBuff <- terra::aggregate(x = LayBuff, dissolve = TRUE)
          currArea <- sumExpanse(LayBuff, unit = "m", transform = FALSE) # NEEDS TO BE METERS. USED BELOW
          message(paste0("Buffered (500m) area for ", Sector, 
                         " -- ", disturbanceEnd, ": ", round(currArea/1000000, 2), " Km2"))
          message(paste0("Percentage of current area: ", round(100*((currArea/1000000)/totalstudyAreaVAreaSqKm), 3), "%."))
        }
        return(Lay)
  })
  names(Connected) <- dPar[["dataName"]]
  # Mash together what has the same name
  # Which Sector has layers that are the same?

  Connected <- lapply(Connected, function(layerList) {
    if (!is.list(layerList)) return(layerList)
    layerList <- Filter(Negate(is.null), layerList)
    if (!length(layerList)) return(layerList)
    layerList <- Filter(function(x) {
      inherits(x, "SpatVector") && tryCatch(nrow(x) > 0L, error = function(...) FALSE)
    }, layerList)
    if (!length(layerList)) return(layerList)
    if (length(layerList) == 1L) return(layerList)
    mergeWidth <- tryCatch(res(rasterToMatchR)[1], error = function(...) NA_real_)
    mergedLayer <- safeMergeSpatVectors(layerList,
                                        context = "connecting layers",
                                        bufferWidth = mergeWidth)
    if (is.null(mergedLayer)) return(layerList)
    mergedName <- names(layerList)
    mergedName <- mergedName[nzchar(mergedName)]
    out <- list(mergedLayer)
    if (length(mergedName)) names(out) <- mergedName[1]
    out
  })
  
  # LAST: Updating and returning
  #############################################

  # Put all updated layers in a list to return
  individualLayers <- c(Enlarged, Generated, Connected)
  
  # Now merge all disturbances (raster format) to avoid creating new disturbances where it has 
  # already been disturbed

  currentDisturbance <- copy(individualLayers)
  # DEPRECATING BELOW: Instead of converting all to raster, convert ras to polys and keep all as vects
  # curDistRas <- lapply(1:length(currentDisturbance), function(index1){
  #   SECTOR <- names(currentDisturbance)[index1]
  #   curDistRas <- lapply(1:length(currentDisturbance[[index1]]), function(index2){
  #     Class <- names(currentDisturbance[[index1]])[index2]
  #     ras <- currentDisturbance[[index1]][[index2]]
  #     if (is.null(ras)){ 
  #       message(paste0(Class, " of ", SECTOR, " is NULL. Returning NULL"))
  #       return(NULL)
  #     }
  #     if (!class(ras) %in% c("RasterLayer", "SpatRaster")){
  #       message(paste0("Converting ", Class, " of ", SECTOR, " from vector to raster..."))
  #       # Now we fasterize so this is used to exclude pixels that have already been disturbed
  #       # As some of these are lines, we will need to buffer them so they become polygons at
  #       # the original resolution they had in the dataset
  #       if (!is(ras, "SpatVector"))
  #         ras <- vect(ras)
  #       if  (terra::geomtype(ras) %in% c("points", "lines")) {
  #         if ("resolutionVector" %in% dPar){
  #           RES <- unique(dPar[disturbanceEnd == Class, resolutionVector])/2
  #         } else RES <- NULL
  #         if (any(is.na(RES),
  #                 is.null(RES))){
  #           message(crayon::red(paste0("resolutionVector in disturbanceParameters",
  #                                      " table was not supplied for a lines file vector. Will ",
  #                                      "default to 15m total (7.5m in each direction).",
  #                                      " If this is not correct, please provide the ",
  #                                      "resolution used for developing the layer ", Class,
  #                                      " (", SECTOR,")")))
  #           RES <- 7.5
  #         }
  #         currDistT <- terra::buffer(x = ras, width = RES)
  #       }  else {
  #         currDistT <- ras
  #       }
  #       currDistSF <- sf::st_as_sf(x = currDistT)
  #       currDistSF$disturbance <- 1
  #       fld <- "disturbance"
  #       currentDisturbanceLay <- suppressWarnings(fasterize::fasterize(sf = st_collection_extract(currDistSF, "POLYGON"),
  #                                                     raster = rasterToMatchR, field = fld))
  #       # NOTE: 
  #       # Seems that the fasterize is not picking up roads from currDistSF, which makes sense considering they are 
  #       # much smaller than the resolution... However, this layer is just so we exclude the pixels that 
  #       # have been disturbed already. This means that lines should be fine, as we have space for more 
  #       # stuff to come in in a 250mx250m pixel that has 7.5 or 15m disturbed. For boo's models, we will 
  #       # use the full layers and buffer before we fasterize (PopGrowth) and do the calculations of decay 
  #       # in the polygon or at a higher resolution. So for now, this solution is fine.
  #       return(rast(currentDisturbanceLay))
  #     } else {
  #       if (is(ras, "RasterLayer")){
  #         message(paste0(Class, " of ", SECTOR, " is a RasterLayer, converting to SpatRaster..."))
  #         ras <- rast(ras)
  #       } else {
  #         message(paste0(Class, " of ", SECTOR, " is either already a SpatRaster or NULL, skipping conversion..."))
  #       }
  #       return(ras)
  #     }
  #   })
  #   if (length(curDistRas) == length(names(currentDisturbance[[index1]]))){
  #     names(curDistRas) <- names(currentDisturbance[[index1]])
  #   } else {
  #      stop("The amount of rasters doesn't match the amount of names. Please debug.")
  #    }
  #   return(curDistRas)
  # })
  # if (length(curDistRas) == length(names(currentDisturbance))){
  #   names(curDistRas) <- names(currentDisturbance)
  # } else {
  #   stop("The amount of rasters doesn't match the amount of names. Please debug.")
  # }
  curDistRas <- lapply(seq_along(currentDisturbance), function(index1){
    SECTOR <- names(currentDisturbance)[index1]
    sectorLayers <- currentDisturbance[[index1]]
    classNames <- names(sectorLayers)
    converted <- lapply(classNames, function(Class){
      ras <- sectorLayers[[Class]]
      if (is.null(ras)) {
        message(paste0(Class, " of ", SECTOR, " is NULL. Returning NULL"))
        return(NULL)
      }
      isEmpty <- FALSE
      if (inherits(ras, "SpatVector")) {
        isEmpty <- length(ras) == 0
      } else {
        vals <- tryCatch(ras[], error = function(...) NULL)
        if (is.null(vals)) {
          isEmpty <- TRUE
        } else {
          vals <- vals[!is.na(vals)]
          isEmpty <- !length(vals) || max(vals, na.rm = TRUE) == 0
        }
      }
      if (isEmpty) {
        message(paste0(Class, " of ", SECTOR, " is empty. Returning NULL"))
        return(NULL)
      }
      if (inherits(ras, c("RasterLayer", "SpatRaster"))){
        message(paste0("Converting ", Class, " of ", SECTOR, " from raster to vector..."))
        if (!inherits(ras, "SpatRaster"))
          ras <- rast(ras)
        ras[ras != 1] <- NA
        currentDisturbanceLay <- terra::as.polygons(ras)
        return(currentDisturbanceLay)
      }
      if (!inherits(ras, "SpatVector")){
        message(paste0(Class, " of ", SECTOR, " is not SpatVector, converting to it..."))
        ras <- terra::vect(ras)
      } else {
        message(paste0(Class, " of ", SECTOR, " is either already a SpatVector or NULL, skipping conversion..."))
      }
      geomType <- tryCatch(terra::geomtype(ras), error = function(...) NA_character_)
      if (!is.na(geomType) && !identical(geomType, "polygons")) {
        baseWidth <- tryCatch(res(rasterToMatch)[1], error = function(...) NA_real_)
        if (is.na(baseWidth) || baseWidth <= 0) baseWidth <- 15
        width <- if (identical(geomType, "points")) baseWidth/2 else baseWidth
        if (isTRUE(getOption("run_scenario.debug", FALSE))) {
          message("[generateDisturbances] buffering ", Class, " of ", SECTOR,
                  " (geom=", geomType, ") to width=", signif(width, 4))
        }
        ras <- terra::buffer(ras, width = width)
        ras <- terra::aggregate(ras, dissolve = TRUE)
      }
      ras
    })
    names(converted) <- classNames
    converted <- Filter(Negate(is.null), converted)
    if (!length(converted)) return(NULL)
    converted
  })
  names(curDistRas) <- names(currentDisturbance)
  curDistRas <- Filter(Negate(is.null), curDistRas)
  # BELOW IS DEPRECATED AS WE ARE NOT DEALING WITH RASTERS ANYMORE
  # Make sure to return which ones are NULL and clean them up.
  # curDistList <- unlist(curDistRas)
  # stk <- rast(curDistList)
  # stkDT <- data.table(pixelIndex = 1:ncell(stk),
  #                     values(stk))
  # stkDT[, currDisturbance := rowSums(.SD, na.rm = TRUE), .SDcols = names(stk)]
  # 
  # newDisturbanceLayers <- terra::setValues(rasterToMatch, stkDT[["currDisturbance"]])
  # 
  # # Cleanup the zeros where is actually NA
  # newDisturbanceLayers[is.na(rasterToMatch)] <- NA
  # 
  # # Some disturbances happened in the same pixel (0.001% of the area), 
  # # but this is a small problem and helps achieve the correct % in the end.  
  # newDisturbanceLayers[newDisturbanceLayers[] > 1] <- 1
  # 
  # if (!is.null(allLaysRas)) {
  #   newDisturbanceLayers[ == 1] <- 1
  # } else stop("Something went wrong. allLaysRas (with all current disturbances) is NULL. Please debug.")
  

  # [ UPDATED ]
  # curDistRas # List of newDisturbanceLayers
  
  # futureDistTb <- table(newDisturbanceLayers[])
  # rtmTable <- table(rasterToMatch[])
  # 
  # message(paste0("Unbuffered (raster) updated disturbance percent of the area (",currentTime,"): ", 
  #                round(100*(futureDistTb["1"]/rtmTable["1"]), 2), "%."))
  # 
  # newDisturbanceLayers[newDisturbanceLayers == 0] <- NA
  # newDisturbanceLayersPol <- terra::as.polygons(newDisturbanceLayers)
  # newDisturbanceLayersPol$totAreaKm2 <- terra::expanse(newDisturbanceLayersPol, 
  #                                                      unit = "km", transform = FALSE) # Unbuffered to 500 but new layers + old
  
  ########################### FINAL LAYERS ###########################

  # individualLayers: mixed rasters and vectors, coming directly from the disturbances generated
   
  # curDistRas: all vectors, unbuffered new disturbances, no old ones here except for seismic lines and settlements
  #             While seismic lines are not lines anymore (they got buffered to increase the disturbance), we still have lines
  #             belonging to roads and pipelines.     
  
  ########################### FINAL LAYERS ###########################
  
  # Now I want to test it with 500m buffer --> For knowldge. The layers below are not 
  # returned!!
  
  if (disturbanceRateRelatesToBufferedArea){
    curDistVcs <- unlist(curDistRas)
    
    curDistVcsAll <- lapply(names(curDistVcs), function(eachVectNm){
      eachVect <- curDistVcs[[eachVectNm]]
      if (length(eachVect) == 0) return(NULL)
      message(paste0("Buffering and/or merging polygons for ", eachVectNm))
      buffVect <- terra::buffer(eachVect, width = 500)
      buffVect <- terra::aggregate(buffVect, dissolve = TRUE)
      return(buffVect)
    })
    curDistVcsAll <- Filter(Negate(is.null), curDistVcsAll)
    mergeWidth <- tryCatch(res(rasterToMatchR)[1], error = function(...) NA_real_)
    newDisturbanceLayers <- safeMergeSpatVectors(curDistVcsAll,
                                                 context = "buffered disturbance vectors",
                                                 bufferWidth = mergeWidth,
                                                 targetCRS = terra::crs(rasterToMatchR))
    if (!is.null(newDisturbanceLayers)) {
      newDisturbanceLayers <- terra::aggregate(newDisturbanceLayers, dissolve = TRUE)
    }
    
    # Now I get the previous layers and do the same
    oldDisturbanceLayers <- createBufferedDisturbances(disturbanceList = disturbanceList, 
                                                                   bufferSize = 500,
                                                                   rasterToMatch = rasterToMatch,
                                                                   studyArea = studyArea,
                                                                   currentTime = currentTime,
                                                                   convertToRaster = FALSE)
    
    newDisturbanceLayers$totAreaKm2 <- expanse_vec(newDisturbanceLayers, 
                                                   unit = "km", transform = FALSE) 
    oldDisturbanceLayers$totAreaKm2 <- expanse_vec(oldDisturbanceLayers, 
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
  
  list(individualLayers = individualLayers, 
       currentDisturbanceLayer = curDistRas)
}
