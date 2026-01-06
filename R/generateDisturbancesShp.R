generateDisturbancesShp <- function(disturbanceParameters,
                                    disturbanceList,
                                    rasterToMatch,
                                    studyArea,
                                    fires,
                                    currentTime,
                                    firstTime,
                                    growthStepGenerating, #not used yet
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

  if (maskWaterAndMountainsFromLines)(Require("spaths"))
  if (useRoadsPackage)(Require("geodata"))
  
  # Extracting layers from previous ones
  # Total study area
  rasterToMatchR <- raster::raster(rasterToMatch)
  probabilityDisturbance <- if (is.null(probabilityDisturbance)) list() else probabilityDisturbance
  studyAreaHash <- digest(studyArea)
  resolveSeismicLengthStats <- function(clusterVec,
                                        fallbackVec,
                                        disturbanceRow,
                                        defaultMean = 1000,
                                        defaultSd   = 300) {
    extractLengths <- function(vec) {
      if (!inherits(vec, "SpatVector")) return(numeric())
      if (tryCatch(nrow(vec) == 0L, error = function(...) TRUE)) return(numeric())
      geomTypes <- tryCatch(unique(terra::geomtype(vec)), error = function(...) character(0))
      if (length(geomTypes) && !any(geomTypes %in% c("lines", "polygons"))) return(numeric())
      len <- numeric()
      if ("calculatedLength" %in% names(vec)) {
        vals <- vec[["calculatedLength"]]
        if (is.list(vals)) vals <- unlist(vals, use.names = FALSE)
        len <- suppressWarnings(as.numeric(vals))
      }
      if (!length(len) || !any(is.finite(len))) {
        len <- tryCatch(as.numeric(terra::perim(vec)), error = function(...) numeric())
      }
      len <- len[is.finite(len) & len > 0]
      len
    }
    statsFromVec <- function(len, sourceTag) {
      if (!length(len)) return(NULL)
      meanLen <- mean(len)
      sdLen <- stats::sd(len)
      if (!is.finite(sdLen) || sdLen <= 0) {
        sdLen <- max(1, meanLen * 0.1)
      }
      list(mean  = as.numeric(meanLen),
           sd    = as.numeric(sdLen),
           lower = max(0, min(len, na.rm = TRUE)),
           upper = max(len, na.rm = TRUE),
           source = sourceTag)
    }
    clusterStats <- statsFromVec(extractLengths(clusterVec), "cluster")
    if (!is.null(clusterStats)) return(clusterStats)
    currentStats <- statsFromVec(extractLengths(fallbackVec), "current")
    if (!is.null(currentStats)) return(currentStats)
    if (!is.null(disturbanceRow) && nrow(disturbanceRow) > 0) {
      meanCol <- "lineLengthMean"
      sdCol   <- "lineLengthSd"
      mu <- if (meanCol %in% names(disturbanceRow))
        suppressWarnings(as.numeric(disturbanceRow[[meanCol]])) else NA_real_
      sigma <- if (sdCol %in% names(disturbanceRow))
        suppressWarnings(as.numeric(disturbanceRow[[sdCol]])) else NA_real_
      if (is.finite(mu) && mu > 0 && is.finite(sigma) && sigma > 0) {
        return(list(mean = mu, sd = sigma, lower = 0, upper = Inf, source = "parameters"))
      }
    }
    list(mean = defaultMean,
         sd   = defaultSd,
         lower = 0,
         upper = defaultMean + 3 * defaultSd,
         source = "default")
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
  clipToStudyArea <- function(vec, area = studyArea) {
    if (!inherits(vec, "SpatVector") || tryCatch(nrow(vec) == 0, error = function(...) TRUE)) return(vec)
    clipped <- tryCatch(terra::intersect(vec, area), error = function(...) NULL)
    if (inherits(clipped, "SpatVector") && tryCatch(nrow(clipped) > 0, error = function(...) FALSE)) {
      return(clipped)
    }
    # fallback to crop if intersect fails or empties unexpectedly
    tryCatch(terra::crop(vec, area), error = function(...) vec)
  }
  bindSpatVectors <- function(lhs, rhs, referenceRaster = rasterToMatchR) {
    lhsValid <- inherits(lhs, "SpatVector") && tryCatch(nrow(lhs) > 0, error = function(...) FALSE)
    rhsValid <- inherits(rhs, "SpatVector") && tryCatch(nrow(rhs) > 0, error = function(...) FALSE)
    if (!lhsValid && !rhsValid) return(NULL)
    if (!lhsValid) return(rhs)
    if (!rhsValid) return(lhs)
    lhsGeom <- tryCatch(unique(terra::geomtype(lhs)), error = function(...) character(0))
    rhsGeom <- tryCatch(unique(terra::geomtype(rhs)), error = function(...) character(0))
    if (!setequal(lhsGeom, rhsGeom)) {
      toPolygons <- function(vec) {
        geom <- tryCatch(unique(terra::geomtype(vec)), error = function(...) character(0))
        if (length(geom) && all(geom == "polygons")) return(vec)
        width <- tryCatch(res(referenceRaster)[1], error = function(...) NA_real_)
        if (!is.finite(width) || width <= 0) width <- 1
        buff <- terra::buffer(vec, width = width/2)
        terra::aggregate(buff, dissolve = TRUE)
      }
      lhs <- toPolygons(lhs)
      rhs <- toPolygons(rhs)
      lhsGeom <- tryCatch(unique(terra::geomtype(lhs)), error = function(...) character(0))
      rhsGeom <- tryCatch(unique(terra::geomtype(rhs)), error = function(...) character(0))
    }
    tryCatch(rbind(lhs, rhs), error = function(e) {
      warning(paste0("bindSpatVectors: failed to merge geometries (lhs=", paste(lhsGeom, collapse = "/"),
                     ", rhs=", paste(rhsGeom, collapse = "/"), "): ", conditionMessage(e)),
              immediate. = TRUE)
      lhs
    })
  }
  estimateSeismicGridCount <- function(targetAreaM2,
                                       lengthStats,
                                       spacingRange = c(50, 100),
                                       lineWidth = 6,
                                       maxGrids = 10000L) {
    if (!is.finite(targetAreaM2) || targetAreaM2 <= 0) return(1L)
    spacingVals <- spacingRange[is.finite(spacingRange) & spacingRange > 0]
    meanSpacing <- if (length(spacingVals)) mean(spacingVals) else 75
    meanLength <- lengthStats$mean
    if (!is.finite(meanLength) || meanLength <= 0) meanLength <- 1000
    nCols <- max(1, ceiling(meanLength/meanSpacing))
    nRows <- max(1, ceiling(meanLength/meanSpacing))
    meanCross <- 1 # expected value from sample(0:2, 1)
    estLength <- meanLength * (nCols + nRows + meanCross)
    estAreaPerGrid <- max(estLength * lineWidth, 1)
    estGrids <- ceiling(targetAreaM2/estAreaPerGrid)
    estGrids <- max(1L, as.integer(estGrids))
    if (estGrids > maxGrids) {
      warning("Estimated seismicLineGrids (", estGrids, ") exceeds maxGrids (", maxGrids,
              "); capping to avoid runaway grid generation.", immediate. = TRUE)
      estGrids <- as.integer(maxGrids)
    }
    estGrids
  }
  resolveSeismicClusterStep <- function(step, clusterAreas, targetArea,
                                        totalAreaSqm = NA_real_,
                                        targetIterations = 25L) {
    userSupplied <- !is.null(step) && length(step) == 1L && !is.na(step)
    if (userSupplied) return(max(1L, as.integer(round(step))))
    meanArea <- if (length(clusterAreas)) mean(clusterAreas, na.rm = TRUE) else NA_real_
    if (!is.finite(meanArea) || meanArea <= 0) return(1L)
    estNeeded <- targetArea/meanArea
    auto <- max(1L, as.integer(round(estNeeded/targetIterations)))
    if (!is.finite(auto) || auto <= 0) auto <- 1L
    if (is.finite(totalAreaSqm) && totalAreaSqm > 0) {
      targetFrac <- targetArea/totalAreaSqm
      targetFrac <- min(1, targetFrac)
      scaleFactor <- min(1, max(0.01, 1/(1 + targetFrac * 25)))
    } else {
      scaleFactor <- 1/20
    }
    scaled <- max(1L, as.integer(round(auto * scaleFactor)))
    attr(scaled, "scaleFactor") <- scaleFactor
    scaled
  }
  
  if (!is(studyArea, "SpatVector"))
    studyArea <- terra::vect(studyArea)
  
  studyArea <- terra::project(x = studyArea, y = terra::crs(rasterToMatch))
  uniStudyArea <- terra::aggregate(studyArea)
  totalstudyAreaVAreaSqKm <- sumExpanse(uniStudyArea, unit = "km", transform = FALSE)
  totalstudyAreaVAreaSqm <- sumExpanse(uniStudyArea, unit = "m", transform = FALSE)
  totNPix <- sum(rasterToMatch[], na.rm = TRUE)
  calculatedPixelSizem2 <- totalstudyAreaVAreaSqm/totNPix # in m2
  growthStepEnlargingLinesRaw <- growthStepEnlargingLines
  growthStepEnlargingLines <- if (is.null(growthStepEnlargingLines) || is.na(growthStepEnlargingLines)) 1 else growthStepEnlargingLines
  
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
          message(
            paste0(
              "Layer of dataName = ", Sector,
              " and disturbanceOrigin = ", ORIGIN,
              " is missing; skipping Enlarging for this class."
            )
          )
          return(NULL)
        } # gracefully skip if no current layer available for enlarging
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
            if ("resolutionVector" %in% names(dPar)){
              if (originalForm == "lines")
                RES <- dParOri[["resolutionVector"]]
              if (isTRUE(originalForm == "points")) RES <- RES/2
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
          LayBuff <- clipToStudyArea(LayBuff)
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
  }
  #############################################
  
  
  # SECOND: Generating
  #############################################
  message(crayon::white("Generating disturbances..."))
  # First, mask the current disturbances on the potential layer so we don't choose them again
  # Get the resolution for lines
  RES_list <- disturbanceParameters[["resolutionVector"]]
  if (is.list(RES_list) && any(lengths(RES_list) > 1)) {
    stop(paste0("Different resolutions have not yet been implemented.",
                "Please modify the code to allow for it"))
  }
  RES <- suppressWarnings(as.numeric(unlist(RES_list, use.names = FALSE)))
  if (length(RES) == 0L) {
    RES <- tryCatch({
      rres <- terra::res(rasterToMatch)
      as.numeric(mean(rres, na.rm = TRUE))
    }, error = function(...) NA_real_)
    if (terra::is.lonlat(rasterToMatch)) {
      warning("rasterToMatch has a lon/lat CRS; buffer width will be interpreted in degrees. ",
              "Consider projecting inputs to a metric CRS (e.g., EPSG:3005).")
    }
  }
  RES <- RES[1]
  if (!is.finite(RES) || RES <= 0) RES <- 7.5
  toSpat <- function(x) {
    if (inherits(x, "SpatVector")) return(x)
    if (inherits(x, "sf")) return(terra::vect(x))
    x
  }
  hasFeatures <- function(x) inherits(x, "SpatVector") &&
    tryCatch(terra::nrow(x) > 0L, error = function(...) FALSE)

  if (is.null(currentDisturbanceLayer)){
    message(paste0("currentDisturbanceLayer is NULL. This is likely the first year of the",
                   " simulation. Creating current disturbed polygons..."))
    # Current layer is null, which indicates this might be the first year of the simulation. I
    # will need to create it then, based on the existing disturbances.
    # Get all current disturbances
    # only pull out the spatial layers, drop anything else (e.g. numeric 1)
    
    # flatten the disturbanceList
    unDL <- unlist(disturbanceList, recursive = FALSE)
    # keep only the “current” layers (drop anything with "potential" in its name)
    validNames <- names(unDL)[ ! grepl("potential", names(unDL)) ]
    unDL <- unDL[ validNames ]
    # from those, keep only spatial objects (drop numerics, etc.)
    currDist <- Filter(function(x) inherits(x, c("SpatVector", "sf", "Spatial", "RasterLayer", "SpatRaster")), unDL)
    # Combine all layers, both polygons and lines, separately
    linesLays <- currDist[sapply(currDist, geomtype) == "lines"]
    polyLays  <- currDist[sapply(currDist, geomtype) == "polygons"]
    # Drop zero-row layers to avoid empty rbind() calls
    if (length(linesLays)) {
      linesLays <- Filter(function(x) {
        inherits(x, "SpatVector") && tryCatch(nrow(x) > 0L, error = function(...) FALSE)
      }, linesLays)
    }
    if (length(polyLays)) {
      polyLays <- Filter(function(x) {
        inherits(x, "SpatVector") && tryCatch(nrow(x) > 0L, error = function(...) FALSE)
      }, polyLays)
    }
    if (length(linesLays)) {
      linesLays <- lapply(linesLays, toSpat)            # keep raw line geometry
    }
    polyLays  <- lapply(polyLays, toSpat)
    hasFeatures <- function(x) inherits(x, "SpatVector") &&
      tryCatch(terra::nrow(x) > 0L, error = function(...) FALSE)
    linesLays <- Filter(hasFeatures, linesLays)
    polyLays  <- Filter(hasFeatures, polyLays)
    # buffered copies of lines for masking only
    maskLinesLays <- lapply(linesLays, function(x) {
      if (inherits(x, "SpatVector")) terra::buffer(x, width = RES) else x
    })
    maskLinesLays <- Filter(hasFeatures, maskLinesLays)
    if (!length(polyLays) && !length(linesLays)) {
      # No current disturbances — proceed with generation using full potential
      allLays <- NULL
    } else {
      linesAndPolys <- Filter(Negate(is.null), c(polyLays, maskLinesLays))
      if (length(linesAndPolys)) {
        msgClasses <- paste(vapply(linesAndPolys, function(x) paste(class(x), collapse = "/"), character(1)), collapse = ", ")
        message("Combining current disturbance layers (", length(linesAndPolys), ") classes: ", msgClasses)
      }
      if (length(linesAndPolys) == 0L) {
        allLays <- NULL
      } else {
        sfList <- lapply(linesAndPolys, function(x) {
          if (inherits(x, "SpatVector")) x <- sf::st_as_sf(x)
          if (inherits(x, "sf")) x[, 0, drop = FALSE] else NULL
        })
        sfList <- Filter(function(x) inherits(x, "sf") && nrow(x) > 0L, sfList)
        if (!length(sfList)) {
          allLays <- NULL
        } else {
          combinedSF <- do.call(rbind, sfList)
          allLays <- terra::vect(combinedSF)
        }
      }
    }
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
            curDistVcsAll <- lapply(curDistVcsAll, toSpat)
            curDistVcsAll <- Filter(hasFeatures, curDistVcsAll)
            if (!length(curDistVcsAll)) {
              allLays <- NULL
            } else {
              sfList <- lapply(curDistVcsAll, function(x) {
                if (inherits(x, "SpatVector")) x <- sf::st_as_sf(x)
                if (inherits(x, "sf")) x[, 0, drop = FALSE] else NULL
              })
              sfList <- Filter(function(x) inherits(x, "sf") && nrow(x) > 0L, sfList)
              if (!length(sfList)) {
                allLays <- NULL
              } else {
                message("Combining existing disturbance vectors (", length(sfList), ")")
                allLays <- terra::vect(do.call(rbind, sfList))
              }
            }
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
            
            # guard against empty layers
            if (is.null(potLay) || (inherits(potLay, "SpatVector") && terra::nrow(potLay) == 0)) {
              message("Forestry: potential layer is empty or missing; skipping forestry generation for this year.")
              return(NULL)
            }
            
            # First: Select only productive forests
            # Previous pixel strategy for Generating Disturbances
            potField <- dParOri[["potentialField"]]
            if (any(is.na(potField), 
                    potField == "",
                    is.null(potField))){ 
              # If NA, it doesn't matter, but need to create a 
              # field so not the whole thing becomes one big 1 map
              potField <- "Potential"
              potLay[["Potential"]] <- 1
            }
            if (!potField %in% names(potLay)) {
              potLay[[potField]] <- 1
            }
            validIdx <- which(!is.na(potLay$ORIGIN) & potLay$ORIGIN < (currentTime - 50))
            if (!length(validIdx)) {
              message("Forestry: no features within the ORIGIN window; skipping for this year.")
              return(NULL)
            }
            potLay <- potLay[validIdx]
            areas <- tryCatch(expanse_vec(potLay, unit = "m", transform = FALSE), error = function(e) NA_real_)
            if (any(is.na(areas) | areas <= 0)) {
              dropIdx <- which(is.na(areas) | areas <= 0)
              if (length(dropIdx))
                potLay <- potLay[-dropIdx]
            }
            if (terra::nrow(potLay) == 0) {
              message("Forestry: no polygons remain after removing zero-area features; skipping for this year.")
              return(NULL)
            }
            potVals <- tryCatch(terra::values(potLay, mat = FALSE)[[potField]], error = function(e) NULL)
            if (is.null(potVals)) potVals <- rep(1, terra::nrow(potLay))
            if (is.list(potVals)) {
              potVals <- vapply(potVals, function(x) {
                if (length(x)) {
                  suppressWarnings(as.numeric(x[[1]]))
                } else {
                  NA_real_
                }
              }, numeric(1))
            }
            potValsNum <- suppressWarnings(as.numeric(potVals))
            if (all(is.na(potValsNum))) {
              potValsNum <- rep(1, terra::nrow(potLay))
            } else if (any(is.na(potValsNum))) {
              potValsNum[is.na(potValsNum)] <- 0
            }
            potLay[[potField]] <- potValsNum
            
            potLayF <- tryCatch(
              terra::rasterize(potLay, terra::rast(rasterToMatchR), field = potField,
                               background = NA, touches = TRUE),
              error = function(e) {
                message(paste("Forestry: rasterize failed ->", conditionMessage(e),
                              "Skipping forestry generation this iteration."))
                return(NULL)
              }
            )
            if (is.null(potLayF)) return(NULL)
            # Second: remove fires
            if (!is.null(fires)) {
              fireVals <- tryCatch(fires[], error = function(e) NULL)
              if (!is.null(fireVals)) {
                potVals <- potLayF[]
                burnMask <- fireVals == 1
                if (any(burnMask, na.rm = TRUE)) {
                  potVals[burnMask & !is.na(burnMask)] <- NA
                  potLayF[] <- potVals
                }
              }
            }
            # Convert the fire raster to polygons
            potLay <- terra::as.polygons(potLayF, values = TRUE, na.rm = TRUE)
            
            #safe exit if nothing left ----
            if (terra::nrow(potLay) == 0) {
              message("Forestry: no potential after fire masking; skipping for this year.")
              return(NULL)
            }
            
            # ensure the attribute is named "Potential"
            if (!"Potential" %in% names(potLay)) {
              if (length(names(potLay)) == 1) names(potLay) <- "Potential" else potLay[["Potential"]] <- 1
            }
            
            # drop non-positive/NA potential just in case
            potNm <- names(potLay)[1L]  # "Potential" after the block above
            potLay <- terra::subset(potLay, subset = !is.na(potLay[[potNm]]) & potLay[[potNm]] > 0)
            
            # then keep the existing aggregate(dissolve) logic
            if (NROW(potLay) > 1)
              potLay <- terra::aggregate(potLay, by = "Potential", dissolve = TRUE, count = FALSE)
            
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
          if (!"Potential" %in% names(potLay)) { # If we don't have Potential, we need to create it
            potLay[["Potential"]] <- 1
          }
          if (NROW(potLay) > 1)
            potLay <- aggregate(potLay, by = aggby, dissolve = TRUE, count = FALSE)
          if (length(potLay) == 0){# In case the cropped area doens't have anything
            message(paste0("The potential area for ", Sector, " class ", ORIGIN, " is NULL.",
                           " Likely cropped out from studyArea. Returning NULL."))
            browser()
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
          
          useBufferedArea <- isTRUE(disturbanceRateRelatesToBufferedArea) || identical(ORIGIN, "seismicLines")
          # Before I do the iterations I wanna know what is the currently disturbed area for this specific disturbance type
          if (all(!is.null(Lay),
                  useBufferedArea)){
            if (is(Lay, "RasterLayer"))
              Lay <- rast(Lay)
            if (is(Lay, "SpatRaster")){ # Need to convert to polygon for area
              Lay[Lay == 0] <- NA # Otherwise buffers weird places!
              Lay <- terra::as.polygons(Lay)
            }
            LayBuff <- terra::buffer(Lay, width = 500) # Need to aggregate to avoid double counting!
            LayBuff <- clipToStudyArea(LayBuff)
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
          # --- before the while: clamp request & set iteration guards ---
          nAvail <- nrow(potLay)
          
          # Total potential area that could ever be used
          availArea <- suppressWarnings(
            sumExpanse(terra::aggregate(potLay, dissolve = TRUE),
                       transform = FALSE, unit = "m")
          )
          if (is.finite(availArea) && !is.na(availArea) && availArea < expectedNewDisturbAreaSqM) {
            warning(sprintf(
              "Requested area (%.2f m²) exceeds available potential (%.2f m²); clamping.",
              as.numeric(expectedNewDisturbAreaSqM), as.numeric(availArea)
            ))
            expectedNewDisturbAreaSqM <- availArea
          }
          # Hard stop so we never loop forever even if nothing changes
          maxITR   <- max(5L, nAvail + 5L)
          prevTotal <- -Inf
          tol       <- 1e-6
          
          # Make sure we select the areas to disturb in a way we have enough space to disturb the necessary sizes
          ITR <- 1
          totalAreaAvailable <- 0
          
          # 1.2. We need to choose the best places:
          while (totalAreaAvailable < expectedNewDisturbAreaSqM) {
            # Get the possible places for the disturbances. Best potential are the places with max values
            potVals <- terra::values(potLay, mat = FALSE)[, "Potential", drop = TRUE]
            if (is.list(potVals))   potVals <- unlist(potVals, use.names = FALSE)
            if (is.factor(potVals)) potVals <- as.character(potVals)
            potVals <- suppressWarnings(as.numeric(potVals))
            
            # guard
            if (!length(potVals) || all(is.na(potVals))) {
              warning("No usable numeric values in 'Potential' after extraction. Stopping.")
              break
            }
            
            # keep sort if you use it elsewhere; row-windowing is based on feature count
            valuesAvailable <- sort(potVals)
            
            if (ORIGIN %in% siteSelectionAsDistributing) {
              rowsToChoose <- seq_len(nAvail)
            } else {
              rowsToChoose <- max(1, nAvail - (ITR - 1)) : nAvail
            }
            if (isTRUE(rowsToChoose == 0L) && target_m2 > 0) {
              # ensure at least one feature if any potential exists
              rowsToChoose <- 1L
            }
            potLayTop <- potLay[rowsToChoose]
            
            # always start with full potential
            potLayTopValid <- potLayTop
            
            # Now exclude where disturbance already exists
            # Use intersect to see if they overlap.
            if (ORIGIN != "seismicLines") # For seismic lines, intersection happens below
              message(paste0("Intersecting existing disturbances with potential for development for ",
                             Sector, " -- ", ORIGIN))
            if (Sector == "forestry") {
              # only try to erase if there _are_ existing disturbances
              if (inherits(allLays, "SpatVector") && nrow(allLays) > 0) {
                # remove any existing disturbances
                croppedAllLays <- reproducible::cropInputs(allLays, potLayTop)
                potLayTopValid <- terra::erase(potLayTop, croppedAllLays)
              } else {
                # no existing disturbances → full potential, then exit loop
                potLayTopValid <- potLayTop
                totalAreaAvailable <- sumExpanse(
                  terra::aggregate(potLayTopValid, dissolve = TRUE),
                  unit = "m", transform = FALSE
                )
                break
              }
            } else if (ORIGIN != "seismicLines") { # If "seismicLines", we don't need to erase as they can overlap
              if (is.null(allLays) || !inherits(allLays, "SpatVector") || nrow(allLays) == 0) {
                potLayTopValid <- potLayTop
              } else {
                doIntersect <- terra::intersect(allLays, potLayTop)
                if (inherits(doIntersect, "SpatVector") && nrow(doIntersect) > 0) {
                  croppedAllLays <- reproducible::cropInputs(allLays, potLayTop)
                  potLayTopValid <- terra::erase(potLayTop, croppedAllLays)
                } else {
                  potLayTopValid <- potLayTop
                }
              }
            } else {
              potLayTopValid <- potLayTop
            }
            # If everything got erased, fall back to the un-erased potential; if that
            # is also empty, stop.
            if (!inherits(potLayTopValid, "SpatVector") || nrow(potLayTopValid) == 0) {
              warning(paste0(
                "No valid potential features remain after erasing/intersection for ",
                Sector, " -- ", ORIGIN,
                ". Falling back to full potential for this iteration."
              ), immediate. = TRUE)
              potLayTopValid <- potLayTop
              if (!inherits(potLayTopValid, "SpatVector") || nrow(potLayTopValid) == 0) {
                warning("Fallback potential is also empty. Stopping.")
                return(NULL)
              }
            }

            # Now check if there is enough space for the expected disturbance!
            totalAreaAvailable <- tryCatch(
              sumExpanse(terra::aggregate(potLayTopValid, dissolve = TRUE),
                         transform = FALSE, unit = "m"),
              error = function(e) {
                warning("Failed to compute expanse; stopping. ", conditionMessage(e))
                NA_real_
              }
            )
            
            # --- termination guard: no progress or too many iterations ---
            if (!is.finite(totalAreaAvailable) || is.na(totalAreaAvailable)) {
              warning("Total available area could not be computed (NA/Inf). Stopping.")
              break
            }
            if (ITR > maxITR || abs(totalAreaAvailable - prevTotal) < tol) {
              warning(sprintf(
                "Stopping after %d iterations: no more usable potential (%.2f m² of %.2f m² requested).",
                ITR - 1, as.numeric(totalAreaAvailable), as.numeric(expectedNewDisturbAreaSqM)
              ))
              break
            }
            prevTotal <- totalAreaAvailable
            
            ITR <- ITR + 1
          }
          
          # If probabilityDisturbance is NOT provided, calculate
          if (all(ORIGIN %in% siteSelectionAsDistributing,
                  is.null(probabilityDisturbance[[ORIGIN]]))){ # NOTE: Potentiall yery time consuming!
            message(paste0("probabilityDisturbance for ", ORIGIN, " is NULL. Calculating from data..."))
            # 1. Extract the total area of each polygon type
            potLayTopValid$areaPerPoly <- expanse_vec(potLayTopValid, unit = "m", transform = FALSE)
            # 2. Sum all to calculate the total area of all polys
            totAreaPolys <- sum(potLayTopValid$areaPerPoly)
            # 3. Get the total disturbance area for each polygon type       
            areaDistPerPoly <- terra::intersect(Lay, potLayTopValid) # This is super time demanding. 
            # Allow to pass a proportion. This would also avoid completely selecting all polygons if passed! 
            # 4. Calculate the total percentage of the disturbance per polygon type
            buffSeisL <- terra::buffer(areaDistPerPoly, width = 3) # Need to aggregate to avoid double counting!
            areaDistPerPoly$area <- expanse_vec(buffSeisL, unit = "m", transform = FALSE)
            
            # --- Ensure polygons and required attributes for area accounting ----
            if (inherits(areaDistPerPoly, "SpatVector")) {
              gt <- terra::geomtype(areaDistPerPoly)[1]
              if (gt == "lines") {
                # use precomputed buffered layer if available (faster and consistent)
                if (exists("buffSeisL", inherits = TRUE) &&
                    inherits(buffSeisL, "SpatVector") &&
                    terra::geomtype(buffSeisL)[1] == "polygons") {
                  areaDistPerPoly <- buffSeisL
                } else {
                  # fallback: buffer now; width in METRES from params
                  w <- if (!is.null(getOption("seismic_erase_width"))) getOption("seismic_erase_width") else 50
                  if (exists("seismic_erase_width", inherits = TRUE)) w <- seismic_erase_width
                  areaDistPerPoly <- terra::buffer(areaDistPerPoly, width = w)
                }
              }
              # keep Potential and compute area by class
              if (!"Potential" %in% names(areaDistPerPoly)) areaDistPerPoly$Potential <- 1L
              areaDistPerPoly <- terra::aggregate(areaDistPerPoly, by = "Potential", fun = sum, dissolve = TRUE)
              areaDistPerPoly$area <- expanse_vec(areaDistPerPoly)
            }
            
            # Build the table explicitly (terra versions differ on `geom` arg)
            df <- sf::st_drop_geometry(sf::st_as_sf(areaDistPerPoly))
            areaDT <- data.table::as.data.table(df[, c("Potential","area")])
            
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
          if (ORIGIN == "seismicLines") {
            if (firstTime) {
              message("First time generating seismicLines, proceeding with clustering...")
              out_createCropLay <- Cache(
                createCropLayFinalYear1,
                Lay                    = Lay,
                potLayTopValid         = potLayTopValid,
                runClusteringInParallel= runClusteringInParallel,
                clusterDistance        = clusterDistance,
                studyAreaHash          = studyAreaHash
              )
              cropLayFinal   <- out_createCropLay$lines
              potLayTopValid <- out_createCropLay$availableArea
              if (!inherits(potLayTopValid, "SpatVector") || nrow(potLayTopValid) == 0L) {
                potLayTopValid <- potLay
              }
              finalPotLay <- potLayTopValid
              
              if (inherits(cropLayFinal, "SpatVector") && nrow(cropLayFinal) > 0) {
                suppressWarnings(terra::writeVector(
                  cropLayFinal,
                  file.path(outputsFolder, paste0("seismicLinesYear", currentTime, "_", studyAreaHash, ".shp")),
                  overwrite = TRUE
                ))
              } else {
                warning("First-year seismic clustering produced no lines; skipping shapefile write.")
              }
            } else {
              # reuse last year's clustered lines (with Pot_Clus)
              if (
                !is.null(currentDisturbanceLayer) &&
                !is.null(currentDisturbanceLayer[[Sector]]) &&
                !is.null(currentDisturbanceLayer[[Sector]][[ORIGIN]])
              ) {
                cropLayFinal <- currentDisturbanceLayer[[Sector]][[ORIGIN]]
              } else {
                cropLayFinal <- Lay
                warning("No prior seismicLines in currentDisturbanceLayer; falling back to raw `Lay` (no Pot_Clus).")
              }
            }
            if (!exists("finalPotLay", inherits = FALSE) ||
                !inherits(finalPotLay, "SpatVector") ||
                nrow(finalPotLay) == 0L) {
              finalPotLay <- potLayTopValid
            }
            lengthStats <- resolveSeismicLengthStats(
              clusterVec    = cropLayFinal,
              fallbackVec   = Lay,
              disturbanceRow = dParOri
            )
            Mean <- lengthStats$mean
            Sd   <- lengthStats$sd
            lengthLower <- lengthStats$lower
            lengthUpper <- lengthStats$upper
            if (!is.finite(lengthLower) || lengthLower < 0) lengthLower <- 0
            if (!is.finite(lengthUpper) || lengthUpper <= lengthLower) {
              lengthUpper <- Inf
            }
            if (isTRUE(getOption("run_scenario.debug", FALSE))) {
              message(
                "[generateDisturbancesShp] seismic line length stats source=",
                lengthStats$source, ", mean=",
                signif(Mean, 3), ", sd=", signif(Sd, 3),
                ", lower=", signif(lengthLower, 3),
                ", upper=", if (is.finite(lengthUpper)) signif(lengthUpper, 3) else "Inf"
              )
            }
            if (!useClusterMethod) {
              if (is.null(seismicLineGrids) ||
                  !is.finite(seismicLineGrids) ||
                  is.na(seismicLineGrids) ||
                  seismicLineGrids <= 0) {
                seismicLineGrids <- estimateSeismicGridCount(
                  targetAreaM2 = expectedNewDisturbAreaSqM,
                  lengthStats  = lengthStats
                )
                rateMsg <- if (is.finite(disturbRate)) paste0(signif(disturbRate, 3), "%") else "unknown"
                message(crayon::green(paste0(
                  "Auto-derived seismicLineGrids = ", seismicLineGrids,
                  " using study-area seismic lines (mean length ~", round(Mean, 1),
                  " m), rate ", rateMsg,
                  " over ", dParOri[["disturbanceInterval"]], " year interval."
                )))
              } else {
                seismicLineGrids <- max(1L, as.integer(round(seismicLineGrids)))
                message(crayon::blue(paste0(
                  "Using provided seismicLineGrids = ", seismicLineGrids,
                  " (no auto-estimation applied)."
                )))
              }
            } else {
              message(crayon::blue(paste0(
                "useClusterMethod is TRUE; seismicLineGrids = ",
                ifelse(is.null(seismicLineGrids), "NULL", seismicLineGrids),
                " (grid count not used in clustering workflow)."
              )))
            }
          }
          
          # 1. Make the iteration, and while the area is not achieved, continue 
          areaChosenTotal <- 0
          IT <- 1
          newDisturbs <- NULL
          alreadyReduced <- FALSE
          lineClusterStepAutoMsg <- TRUE
          
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
            # normal size, we only disturb what is expected. Ensure we don't depend on bare rtnorm.
            expr <- dParOri[["disturbanceSize"]]
            # Ensure we always call the fully-qualified msm::rtnorm
            expr <- gsub("(?<!msm::)rtnorm\\(", "msm::rtnorm(", expr, perl = TRUE)
            expr <- gsub("msm::msm::", "msm::", expr, fixed = TRUE)
            Size <- round(eval(parse(text = expr)), 0)
            remaining <- max(0, expectedNewDisturbAreaSqM - areaChosenTotal)
            if (remaining <= 0) break
            if (Size > remaining) Size <- remaining
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
                while (is.list(cropLayFinal) && length(cropLayFinal) == 1) {
                  cropLayFinal <- cropLayFinal[[1]]
                }
                
                if (!inherits(cropLayFinal, "SpatVector") || nrow(cropLayFinal) == 0) {
                  warning("Seismic duplication received no clustered lines; returning empty layer.", immediate. = TRUE)
                  newDisturbs <- cropLayFinal
                  areaChosenTotal <- expectedNewDisturbAreaSqM
                  break
                }
                
                if (geomtype(cropLayFinal) != "lines"){
                  debug_dir <- file.path(getwd(), "scratch", "debug")
                  dir.create(debug_dir, recursive = TRUE, showWarnings = FALSE)
                  base_name <- paste0("cropLayFinal_not_lines_", currentTime, "_", Sector, "_", ORIGIN)
                  rds_path  <- file.path(debug_dir, paste0(base_name, ".rds"))
                  gpkg_path <- file.path(debug_dir, paste0(base_name, ".gpkg"))
                  summary_obj <- list(
                    sector = Sector,
                    origin = ORIGIN,
                    currentTime = currentTime,
                    firstTime = firstTime,
                    class = class(cropLayFinal),
                    geomtype = tryCatch(unique(terra::geomtype(cropLayFinal)), error = function(...) NA_character_),
                    nrow = tryCatch(nrow(cropLayFinal), error = function(...) NA_integer_),
                    names = tryCatch(names(cropLayFinal), error = function(...) NULL)
                  )
                  try(saveRDS(list(meta = summary_obj, object = cropLayFinal), rds_path), silent = TRUE)
                  try(terra::writeVector(cropLayFinal, gpkg_path, filetype = "GPKG", overwrite = TRUE), silent = TRUE)
                  warning(paste0("Seismic duplication expected lines but got geomtype ",
                                 paste(summary_obj$geomtype, collapse = "/"),
                                 ". Dumped to ", rds_path, " and ", gpkg_path, " for inspection; returning empty layer."),
                          immediate. = TRUE)
                  newDisturbs <- cropLayFinal
                  areaChosenTotal <- expectedNewDisturbAreaSqM
                  break
                }
                
                # 0. Calculate the area of each line buffered at the scale used for targets
                # Use 500 m when rates relate to buffered area; otherwise keep the 3 m footprint.
                buffWidth <- if (isTRUE(disturbanceRateRelatesToBufferedArea)) 500 else 3
                areaCol <- "buffAreaM2"
                cropLayFinal[[areaCol]] <- expanse_vec(buffer(x = cropLayFinal, width = buffWidth))
                
                # 1. Add a column in how much each cluster represents in the total in m2 -- not by individual line!
                cropLayFinalDT <- as.data.table(as.data.frame(cropLayFinal))
                
                # --- Standardize types / names to avoid year-to-year drift ---
                # Potential must be integer (not factor/character)
                if (is.factor(cropLayFinalDT$Potential)) {
                  cropLayFinalDT[, Potential := as.integer(as.character(Potential))]
                } else if (is.character(cropLayFinalDT$Potential)) {
                  suppressWarnings(
                    cropLayFinalDT[, Potential := as.integer(Potential)]
                  )
                }
                
                #    Ensure we have a single 'Class' column.
                if (!"Class" %in% names(cropLayFinalDT)) {
                  if ("Class_1" %in% names(cropLayFinalDT)) data.table::setnames(cropLayFinalDT, "Class_1", "Class")
                  if ("Class_2" %in% names(cropLayFinalDT) && !"Class" %in% names(cropLayFinalDT))
                    data.table::setnames(cropLayFinalDT, "Class_2", "Class")
                }
                
                if (!"Pot_Clus" %in% names(cropLayFinalDT)){
                  message("Pot_Clus not found in cropLayFinalDT. Debug")
                  browser()
                }
                cropLayFinalDT[, sumBuffAreaM2 := sum(get(areaCol)), by = "Pot_Clus"] 
                # Need to do by potential as for each potential, cluster numbers are repeated
                totalBuffArea <- sum(cropLayFinalDT[[areaCol]], na.rm = TRUE)
                clusterAreas <- unique(cropLayFinalDT[, .(Pot_Clus, sumBuffAreaM2)])$sumBuffAreaM2
                
                # Guard against division by zero and compute a single % per cluster
                if (totalBuffArea <= 0) {
                  cropLayFinalDT[, PercBuffAreaOfTotalM2 := 0]
                } else {
                  # Compute % once per cluster, then propagate
                  clusPct <- cropLayFinalDT[, .(PercBuffAreaOfTotalM2 = 100 * sumBuffAreaM2[1] / totalBuffArea),
                                            by = "Pot_Clus"]
                  cropLayFinalDT <- clusPct[cropLayFinalDT, on = "Pot_Clus"]
                }
                
                # Robust guard: NA-safe and tolerant of minor floating-point drift
                pct_sum <- sum(
                  unique(cropLayFinalDT[, c("Pot_Clus", "PercBuffAreaOfTotalM2")]$PercBuffAreaOfTotalM2),
                  na.rm = TRUE
                )
                if (!is.na(pct_sum) && pct_sum > 100.001) {
                  message(paste0(
                    "Total contribution of clusters in total area is higher than 100% (",
                    round(pct_sum, 3),
                    "%). Something may be wrong."
                  ))
                  if (isTRUE(getOption("anthroDisturbance.debugClusters", FALSE))) {
                    message("Entering debug mode because anthroDisturbance.debugClusters = TRUE")
                    browser()
                  }
                } # TODO test
                # 2. Choose randomly clusters (with different probabilities) that sum to the total expected new. 
                # sumBuffAreaM2 --> by cluster
                # PercBuffAreaOfTotalM2 --> representation if each cluster over the total
                # Potential --> Represents the highest potential for being chosen.
                # --> Draw for all potentials, the probability a cluster within these will be chosen (higher potential, higher chances)
                # Normalize probabilities to sum to 1
                lineClusterStep <- resolveSeismicClusterStep(
                  step = growthStepEnlargingLinesRaw,
                  clusterAreas = clusterAreas,
                  targetArea = expectedNewDisturbAreaSqM,
                  totalAreaSqm = totalstudyAreaVAreaSqm
                )
                if (is.null(growthStepEnlargingLinesRaw) || is.na(growthStepEnlargingLinesRaw)) {
                  if (isTRUE(lineClusterStepAutoMsg)) {
                    scaleFactor <- attr(lineClusterStep, "scaleFactor")
                    if (is.null(scaleFactor)) scaleFactor <- NA_real_
                    message(crayon::blue(paste0(
                      "Auto-derived growthStepEnlargingLines for seismic clustering = ",
                      lineClusterStep, " (scaleFactor=", signif(scaleFactor, 3),
                      ", aiming to reduce iterations)."
                    )))
                    lineClusterStepAutoMsg <- FALSE
                  }
                }

                probabilities <- unique(cropLayFinalDT$Potential) / sum(unique(cropLayFinalDT$Potential))
                if (length(unique(cropLayFinalDT$Potential)) == 1){
                  sampledClusters <- rep(unique(cropLayFinalDT$Potential), times = lineClusterStep)
                } else {
                  sampledClusters <- sample(unique(cropLayFinalDT$Potential), 
                                            size = lineClusterStep, 
                                            replace = TRUE, 
                                            prob = probabilities)
                }
                selectedClusters <- NULL
                for (uniqueSampClus in unique(sampledClusters)){
                  toChoseFrom <- unique(cropLayFinalDT[Potential == uniqueSampClus, Pot_Clus])
                  howManyINeed <- sum(sampledClusters == uniqueSampClus)
                  if (length(toChoseFrom) == 1){
                    sampledOnes <- rep(toChoseFrom, times = howManyINeed)
                  } else if (length(toChoseFrom) > 1) {
                    sampledOnes <- sample(toChoseFrom,
                                          replace = TRUE,
                                          size = howManyINeed)
                  } else {
                    next
                  }
                  selectedClusters <- c(selectedClusters, sampledOnes)
                }
                # Guard against degenerate selection: no clusters to duplicate
                candLines <- NULL
                if (!is.null(selectedClusters) && length(selectedClusters)) {
                  candLines <- cropLayFinal[cropLayFinal$Pot_Clus %in% selectedClusters, ]
                }
                if (is.null(candLines) ||
                    !inherits(candLines, "SpatVector") ||
                    nrow(candLines) == 0L) {
                  warning(paste0(
                    "Unable to identify candidate seismic line clusters for duplication in Sector '",
                    Sector, "' (ORIGIN = '", ORIGIN, "', time = ", currentTime,
                    "). Skipping new seismic line generation for this disturbance."
                  ), immediate. = TRUE)
                  break
                }
                #TODO Here I can parallelize using future!
                # Clean potentially invalid geometries to avoid low-level C++ aborts in terra
                candLines <- tryCatch({
                  bad <- tryCatch(!terra::is.valid(candLines), error = function(...) rep(FALSE, nrow(candLines)))
                  if (any(bad, na.rm = TRUE)) candLines <- candLines[!bad, ]
                  sf_obj <- tryCatch(sf::st_as_sf(candLines), error = function(e) NULL)
                  if (!is.null(sf_obj)) {
                    sf_obj <- tryCatch(sf::st_make_valid(sf_obj), error = function(...) sf_obj)
                    candLines <- tryCatch(terra::vect(sf_obj), error = function(e) candLines)
                  }
                  candLines
                }, error = function(e) candLines)

                # Wrap simulateLines so refinedStructure can be requested safely; on failure, fall back.
                newLines <- tryCatch(
                  simulateLines(Lines = candLines,
                                distThreshold = clusterDistance,
                                distNewLinesFact = distanceNewLinesFactor,
                                refinedStructure = refinedStructure),
                  error = function(e) {
                    warning(
                      paste0(
                        "simulateLines failed (refinedStructure=", refinedStructure,
                        "); falling back to refinedStructure=FALSE. Error: ",
                        conditionMessage(e)
                      ),
                      immediate. = TRUE
                    )
                    simulateLines(Lines = candLines,
                                  distThreshold = clusterDistance,
                                  distNewLinesFact = distanceNewLinesFactor,
                                  refinedStructure = FALSE)
                  }
                )
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
                  valsToExclude <- potLayTopValid$Potential
                  if (is.list(valsToExclude)) {
                    valsToExclude <- unlist(valsToExclude, recursive = TRUE, use.names = FALSE)
                  }
                  valsToExclude <- suppressWarnings(as.numeric(valsToExclude))
                  valsToExclude <- valsToExclude[is.finite(valsToExclude)]
                  candidates <- setdiff(valuesAvailable, valsToExclude)
                  if (!length(candidates)) {
                    warning(
                      "Unable to identify an additional potential class for seismic line placement; ",
                      "falling back to the remaining potential pool."
                    )
                    candidates <- valuesAvailable
                  }
                  nextBestValue <- max(candidates)
                  rowsToChoose <- which(valuesAvailable == nextBestValue)
                  potLayTopValid <- potLay[rowsToChoose]
                  # 1. UPDATING THE LAYER: 
                  cropLay <- postProcessTo(Lay, potLayTopValid)
                  cropLayBuf <- buffer(cropLay, width = 50)
                  cropLayAg <- aggregate(cropLayBuf, dissolve = TRUE)
                  finalPotLay <- erase(potLayTopValid, cropLayAg)
                  centerPointToAdd <- terra::spatSample(finalPotLay, size = howManyMissing, method = "random")
                  centerPoint <- if (is.null(centerPoint)) centerPointToAdd else rbind(centerPoint, centerPointToAdd)
                }
                # 5. Draw the new grid based on the total length expected
                # 5.1. Draw a square based on the centerPoint, where the distance from point to the lines 
                # is the diagonal of a square of the lineLenght you want.
                
                # print("Currently not using, but should test! Something is likely not working")
                # browser() # HERE is where Mean and SD is used from cropLay
                lineLength <- msm::rtnorm(
                  length(centerPoint),
                  Mean,
                  Sd,
                  lower = lengthLower,
                  upper = lengthUpper
                )
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
                  rowsOfPoints <- if (length(rowsOfPointsL)) do.call(rbind, rowsOfPointsL) else NULL
                  # Cols:
                  colsOfPointsL <- lapply(0:(dim(polRas)[1]-1), function(colIndex){
                    thePair <- c(1, dim(polRas)[2])+(colIndex*dim(polRas)[2])
                    # create the line based on the points by extracting the points based on thePair 
                    theLine <- as.lines(gridPoints[thePair, ])
                    return(theLine)
                  })
                  colsOfPoints <- if (length(colsOfPointsL)) do.call(rbind, colsOfPointsL) else NULL
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
                    crossedLines <- if (length(crossLines)) do.call(rbind, crossLines) else NULL
                  } else {
                    crossedLines <- howManyCross <- NULL
                  }
                  parts <- Filter(Negate(is.null), list(rowsOfPoints, colsOfPoints, crossedLines))
                  if (length(parts)) {
                    gridReady <- do.call(rbind, parts) # Newly created grid
                    return(gridReady)
                  } else {
                    return(NULL)
                  }
                })
                gridReadyB <- if (length(gridReady)) do.call(rbind, gridReady) else NULL
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
              
              cn <- names(centerPoint)
              if (length(cn) == ncol(newDist)) names(newDist) <- cn
              # 1.4. If we relate to buffered area, we need to make an inner buffer, as it is included in 
              # areaChosenTotal.
              if (useBufferedArea){
                # Then, if disturbanceRateRelatesToBufferedArea, we buffer this point to identify what 
                # if the minimum area we can have. Then, we also calculate the area of the new disturbance 
                # and we check if the area of the new disturbance is smaller than the minimum area of the 
                # buffered point. 
                minDistBuff <- terra::buffer(centerPoint, width = 500)
                areaMinDB <- sumExpanse(minDistBuff)
                areaNewDist <- sumExpanse(newDist)
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
            
            # Keep new disturbance inside the study area/potential footprint to avoid oversized polygons
            #newDist <- reproducible::postProcess(newDist, studyArea)
            #if (inherits(potLayTopValid, "SpatVector") && tryCatch(nrow(potLayTopValid) > 0, error = function(...) FALSE)) {
            #  clipped <- tryCatch(terra::intersect(newDist, potLayTopValid), error = function(...) NULL)
            #  if (inherits(clipped, "SpatVector") && tryCatch(nrow(clipped) > 0, error = function(...) FALSE)) {
            #    newDist <- clipped
            #  }
            #}
            
            if (useBufferedArea){
              newDistBuff <- terra::buffer(newDist, width = 500) 
              newDistBuff <- clipToStudyArea(newDistBuff)
              if (nrow(newDistBuff) > 1){
                newDistBuff <- terra::aggregate(newDistBuff)
              }
              areaChosenTotal <- areaChosenTotal + sumExpanse(newDistBuff, unit = "m", transform = FALSE)
            } else { 
              areaChosenTotal <- areaChosenTotal + sumExpanse(newDist, unit = "m", transform = FALSE)
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
              # Crop again because features may extend outside the study area
              toCrop <- if (isTRUE(useClusterMethod)) newLines else newDist
              newDist <- reproducible::postProcess(toCrop, studyArea)
            }
            
            if (IT == 1){
              newDisturbs <- newDist
            } else {
              newDisturbs <- if (is.null(newDisturbs)) newDist else rbind(newDisturbs, newDist)
            }
            IT <- IT + 1
          }
          
          # If we overshot the buffered-area target for polygonal Generating classes,
          # shrink the last generated polygon so that the final 500 m–buffered area
          # matches the expected target more closely. This is applied only for
          # forestry/mining/oilGas polygon disturbances (not seismic lines).
          if (isTRUE(disturbanceRateRelatesToBufferedArea) &&
              !identical(ORIGIN, "seismicLines") &&
              inherits(newDisturbs, "SpatVector")) {
            gt_all <- tryCatch(unique(terra::geomtype(newDisturbs)), error = function(...) character(0))
            if (length(gt_all) && "polygons" %in% gt_all) {
              buff_all <- tryCatch(terra::buffer(newDisturbs, width = 500), error = function(...) NULL)
              buff_all <- clipToStudyArea(buff_all)
              if (inherits(buff_all, "SpatVector") && nrow(buff_all) > 0) {
                featAreas <- sumExpanse(buff_all, unit = "m", transform = FALSE)
                totArea <- suppressWarnings(as.numeric(sum(featAreas, na.rm = TRUE)))
                overshoot <- suppressWarnings(as.numeric(totArea - expectedNewDisturbAreaSqM))
                if (is.finite(overshoot) && overshoot > 0) {
                  lastIdx <- length(featAreas)
                  lastArea <- featAreas[lastIdx]
                  prevArea <- totArea - lastArea
                  targetLast <- max(0, expectedNewDisturbAreaSqM - prevArea)
                  
                  if (targetLast <= 0) {
                    message(paste0("Overshoot for ", Sector, " -- ", ORIGIN,
                                   " is fully attributable to last polygon; dropping it to meet target."))
                    if (nrow(newDisturbs) > 1L) {
                      newDisturbs <- newDisturbs[seq_len(nrow(newDisturbs) - 1L), ]
                      buff_all <- tryCatch(terra::buffer(newDisturbs, width = 500), error = function(...) NULL)
                      areaChosenTotal <- if (inherits(buff_all, "SpatVector") && nrow(buff_all) > 0) {
                        sumExpanse(buff_all, unit = "m", transform = FALSE)
                      } else {
                        0
                      }
                    } else {
                      newDisturbs <- NULL
                      areaChosenTotal <- 0
                    }
                  } else {
                    geomLast <- newDisturbs[nrow(newDisturbs), ]
                    gt_last <- tryCatch(terra::geomtype(geomLast)[1], error = function(...) NA_character_)
                    if (isTRUE(gt_last == "polygons")) {
                      bestGeom <- geomLast
                      bestArea <- lastArea
                      # Rough upper bound for shrink distance based on equivalent-circle radius
                      upper <- suppressWarnings(sqrt(lastArea / pi))
                      if (!is.finite(upper) || upper <= 0) upper <- 1000
                      lower <- 0
                      for (it_shrink in seq_len(12L)) {
                        mid <- (lower + upper) / 2
                        if (!is.finite(mid) || mid <= 0) break
                        shrunk <- tryCatch(terra::buffer(geomLast, width = -mid), error = function(...) NULL)
                        if (is.null(shrunk) || nrow(shrunk) == 0L ||
                            any(is.na(as.vector(terra::ext(shrunk))))) {
                          upper <- mid
                          next
                        }
                        buff_mid <- tryCatch(terra::buffer(shrunk, width = 500), error = function(...) NULL)
                        if (is.null(buff_mid) || nrow(buff_mid) == 0L) {
                          upper <- mid
                          next
                        }
                        areaMid <- suppressWarnings(as.numeric(sumExpanse(buff_mid, unit = "m", transform = FALSE)))
                        if (!is.finite(areaMid) || areaMid <= 0) {
                          upper <- mid
                          next
                        }
                        if (abs(areaMid - targetLast) < abs(bestArea - targetLast)) {
                          bestGeom <- shrunk
                          bestArea <- areaMid
                        }
                        if (areaMid > targetLast) {
                          # still too big → shrink more
                          lower <- mid
                        } else {
                          # too small → shrink less
                          upper <- mid
                        }
                      }
                      # Replace the last feature with the shrunken geometry
                      if (!is.null(bestGeom) && nrow(bestGeom) > 0L) {
                        if (nrow(newDisturbs) > 1L) {
                          keep <- newDisturbs[seq_len(nrow(newDisturbs) - 1L), ]
                          newDisturbs <- rbind(keep, bestGeom)
                        } else {
                          newDisturbs <- bestGeom
                        }
                        buff_all <- tryCatch(terra::buffer(newDisturbs, width = 500), error = function(...) NULL)
                        buff_all <- clipToStudyArea(buff_all)
                        areaChosenTotal <- if (inherits(buff_all, "SpatVector") && nrow(buff_all) > 0) {
                          sumExpanse(buff_all, unit = "m", transform = FALSE)
                        } else {
                          0
                        }
                      }
                    }
                  }
                }
              }
            }
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
        # If nothing was generated (NULL/empty) but an input layer exists, keep it
        # to avoid wiping layers when buffered masking eliminates all potential.
        if ((is.null(newDisturbs) ||
             (inherits(newDisturbs, "SpatVector") && tryCatch(nrow(newDisturbs) == 0, error = function(...) FALSE))) &&
            !is.null(Lay) &&
            inherits(Lay, "SpatVector") &&
            tryCatch(nrow(Lay) > 0, error = function(...) FALSE)) {
          message(paste0("Generated disturbances for ", Sector, " -- ", ORIGIN,
                         " are empty; retaining existing layer."))
          newDisturbs <- Lay
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
        connected <- NULL
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
            connectedOne <- terra::vect(connectedOne)
          }
          connected <- bindSpatVectors(connected, connectedOne)
          endLay <- bindSpatVectors(endLay, connectedOne)
        }          
      } else {
        # Hash target to avoid reusing mismatched cached rasters across different rasterToMatch/studyArea pairs
        targetHash <- reproducible::.robustDigest(list(studyArea = studyArea, rasterToMatch = rasterToMatch))
        targetHash <- paste0(unlist(targetHash), collapse = "_")
        if (isTRUE(useRoadsPackage) || (isTRUE(maskWaterAndMountainsFromLines) && is.null(featuresToAvoid))) {
          DEMraw <- geodata::elevation_30s(country = "CAN", path = Paths[["inputPath"]])
          DEM <- postProcessTo(from = DEMraw, to = rasterToMatch, cropTo = studyArea, 
                               projectTo = rasterToMatch, maskTo = studyArea, 
                               writeTo = file.path(Paths[["inputPath"]], 
                                                   paste0("DEM_", targetHash, ".tif")))
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
              waterraw <- geodata::landcover("water", path = Paths[["inputPath"]])
              water <- postProcessTo(from = waterraw, to = rasterToMatch, cropTo = studyArea, 
                                     projectTo = rasterToMatch, maskTo = studyArea, 
                                     writeTo = file.path(Paths[["inputPath"]], 
                                                         paste0("water_", 
                                                                targetHash, ".tif")))
              wetlandsraw <- geodata::landcover("wetland", path = Paths[["inputPath"]])
              wet <- postProcessTo(from = wetlandsraw, to = rasterToMatch, cropTo = studyArea, 
                                   projectTo = rasterToMatch, maskTo = studyArea, 
                                   writeTo = file.path(Paths[["inputPath"]], 
                                                       paste0("wet_", 
                                                               targetHash, ".tif")))
              # Now we put all layers together and exclude the pixels for building the roads (but without modifying the final maps!)
              DEM[DEM[] < altitudeCut] <- 0
              DEM[DEM[] >= altitudeCut] <- 1
              DEM[wet[] > 0.8] <- 1
              DEM[water[] > 0.8] <- 1
              featuresToAvoid <- copy(DEM)
              if (inherits(featuresToAvoid, "SpatRaster")) {
                featuresToAvoid <- tryCatch(
                  terra::extend(featuresToAvoid, rasterToMatch),
                  error = function(e) {
                    warning("Failed to align featuresToAvoid to rasterToMatch extent: ", conditionMessage(e),
                            immediate. = TRUE, call. = FALSE)
                    featuresToAvoid
                  }
                )
              }
              featuresToAvoid[featuresToAvoid == 0] <- NA
              featuresToAvoid[is.na(rasterToMatch)] <- 1
            }
          }
          # Here comes what is below...
          connected <- NULL
          for (i in 1:NROW(oriLayVect)){
            if(i%%100==0)
              message(paste0("Connecting ", Sector," for year ", currentTime,
                             ": ", i ," of ", NROW(oriLayVect),
                             " (", 100*round(i/NROW(oriLayVect),4),"%)"))
            # avoid error on invalid geometries from cropping
            oi <- oriLayVect[i, ]
            # If terra has makeValid(), use it; otherwise go via sf
            valid_oi <- try(terra::is.valid(oi), silent = TRUE)
            if (inherits(valid_oi, "try-error") || isFALSE(valid_oi)) {
              # fallback through sf to repair
              oi_sf <- sf::st_as_sf(oi)
              if (any(!sf::st_is_valid(oi_sf))) {
                oi_sf <- sf::st_make_valid(oi_sf)
              }
              oi <- terra::vect(oi_sf)
            }
            endLayCropped <- terra::crop(endLay, oi)
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
              # crop the avoidance mask around the current origin
              cost <- terra::rast(rasterToMatch)
              cost[] <- 1
              
              if (inherits(featuresToAvoid, "SpatRaster")) {
                cost[is.na(featuresToAvoid[])] <- NA
              } else if (inherits(featuresToAvoid, "SpatVector")) {
                # rasterize obstacles (polygons) to NA
                forbid <- terra::rasterize(featuresToAvoid, cost, field = 1)
                cost[!is.na(forbid[])] <- NA
              }
              
              # Destination: nearest road/target geometry as a line-to-point mapping
              cEndLay <- terra::nearest(oriLayVect[i, ], endLay, pairs = FALSE, centroids = TRUE, lines = FALSE)
              
              connectedOne <- tryCatch(
                spaths::shortest_paths(
                  rst         = cost,
                  origins     = oriLayVect[i, ],
                  destinations= cEndLay,
                  output      = "lines",
                  show_progress = FALSE
                ),
                error = function(e) NULL
              )
              
              # If nothing came back, skip this origin
              if (!inherits(connectedOne, "SpatVector") || nrow(connectedOne) == 0) next
              
              # Some spaths builds add a 'layer' field; if present, keep only > 0
              if ("layer" %in% names(connectedOne)) {
                connectedOne <- connectedOne[connectedOne$layer > 0, ]
                if (nrow(connectedOne) == 0) next
              }
              
              # Ensure SpatVector and proceed (no ext() checks)
              if (!inherits(connectedOne, "SpatVector")) connectedOne <- terra::vect(connectedOne)
              
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
              connectedOne <- terra::vect(connectedOne)
            }
            connected <- bindSpatVectors(connected, connectedOne)
            endLay <- bindSpatVectors(endLay, connectedOne)
          }
          message(paste0("For loop for ", Sector," -- ", disturbanceEnd, " finished."))
        }
        
      }
      if (is.null(connected)){ # It's likely because the one or few created disturbances were 
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
        LayBuff <- clipToStudyArea(LayBuff)
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
      merged <- list(if (length(toMerge)) do.call(rbind, toMerge) else NULL)
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
  individualLayers <- c(Enlarged, Generated, Connected)
  
  if (is.null(individualLayers)){
    stop("The study area has no potential for disturbances. Please enlarge or change the area.")
  }
  # Now merge all disturbances (raster format) to avoid creating new disturbances where it has 
  # already been disturbed
  currentDisturbance <- copy(individualLayers)
  
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
        isMAX0ras <- if (is(ras, "SpatVector")) {
          length(ras) == 0
        } else {
          vals <- tryCatch(ras[], error = function(e) NULL)
          if (is.null(vals) || !any(is.finite(vals))) TRUE else max(vals, na.rm = TRUE) == 0
        }
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
  if (firstTime) {
    f <- file.path(outputsFolder, paste0("seismicLinesYear", currentTime, "_", studyAreaHash, ".shp"))
    seismicLinesFirstYear <- if (file.exists(f)) terra::vect(f) else NULL
  } else {
    seismicLinesFirstYear <- NULL
  }
  
  individualLayers <- cleanupList(individualLayers, 
                                  outter = TRUE, 
                                  inner = TRUE, 
                                  cleanEmpty = TRUE, 
                                  nullEmpty = FALSE)
  
  return(list(individualLayers = individualLayers, 
              currentDisturbanceLayer = curDistRas,
              seismicLinesFirstYear = seismicLinesFirstYear))
}
