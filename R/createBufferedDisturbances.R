#' Buffer & combine multiple disturbance layers (vector + raster) into one output
#'
#' @param disturbanceList   A named list of disturbance‐type sub‐lists. Each sub‐list
#'                          has elements that are either SpatVector/SpatRaster or RasterLayer.
#'                          Any element whose name contains “potential” will be excluded.
#' @param bufferSize        Numeric: buffer distance (in the same CRS units as `studyArea`).
#' @param rasterToMatch     A terra::SpatRaster (template) defining extent/resolution/CRS.
#' @param studyArea         A terra::SpatVector polygon of the entire study area (for area %).
#' @param currentTime       (NOT USED) Placeholder to match module signature.
#' @param convertToRaster   Logical: if TRUE, return a raster; if FALSE, return a combined polygon.
#'
#' @return  If convertToRaster=TRUE: a SpatRaster with 1 = buffered disturbance, 0 = no disturbance, NA outside extent.
#'          If convertToRaster=FALSE: a SpatVector polygon of all buffered disturbances (dissolved).
#' @export
createBufferedDisturbances <- function(disturbanceList,
                                       bufferSize,
                                       rasterToMatch,
                                       studyArea,
                                       currentTime,
                                       convertToRaster) {

  if (is(rasterToMatch, "SpatRaster"))
    rasterToMatchR <- raster::raster(rasterToMatch) else
      rasterToMatchR <- rasterToMatch
    
    studyAreaTotalArea <- terra::expanse(studyArea, transform = FALSE, unit = "km")
    
    # Unlist and filter out potential layers
    unlistedDL <- unlist(disturbanceList)
    nonPotentialLayers <- names(unlistedDL)[!grepl("potential", names(unlistedDL), fixed = FALSE)]
    
    # Early return if only potential layers and vector mode
    if (length(nonPotentialLayers) == 0) {
      if (convertToRaster) {
        # Return a zero-filled raster that matches the template
        zeroRaster <- rasterToMatch
        zeroRaster[] <- 0
        message("Returning all-zero raster due to no non-potential disturbance layers")
        return(zeroRaster)
      } else {
        stop("No non-potential disturbance layers")
      }
    }
    
    allBuffered <- lapply(nonPotentialLayers, function(INDEX){
      lay <- unlistedDL[[INDEX]]
      skipBuffer <- FALSE  # Initialize flag
      
      # Skip non-spatial layers
      if (!inherits(lay, c("SpatVector", "SpatRaster", "RasterLayer"))) {
        warning("skipping non-spatial layer: ", INDEX)
        return(NULL)
      }
      
      # Handle empty vector
      if (inherits(lay, "SpatVector") && nrow(lay) == 0) {
        message(paste0("Layer ", INDEX, " is empty. Returning NULL."))
        return(NULL)
      }
      
      # Convert non-polygon vectors to polygons
      if (inherits(lay, "SpatVector")) {
        geom_type <- unique(terra::geomtype(lay))
        if (!all(geom_type %in% c("polygons", "multipolygons"))) {
          if (bufferSize == 0) {
            lay <- tryCatch(
              terra::buffer(lay, width = 1e-10),
              error = function(e) {
                warning("Geometry conversion failed for ", INDEX, ": ", e$message)
                return(NULL)
              }
            )
          } else {
            lay <- tryCatch(
              terra::buffer(lay, width = bufferSize),
              error = function(e) {
                warning("Buffering failed for ", INDEX, ": ", e$message)
                return(NULL)
              }
            )
          }
          skipBuffer <- TRUE  # Skip later buffering
        }
      }
      
      # Rest of processing remains the same...
      whichClass <- class(lay)
      
      if (whichClass %in% c("SpatRaster", "RasterLayer")) {
        if (is(lay, "RasterLayer"))
          lay <- rast(lay)
        if (sum(lay[], na.rm = TRUE) == 0) {
          message(paste0("Layer ", INDEX, " seems empty. Returning NULL."))
          return(NULL)
        }
        lay[lay != 1] <- NA
        lay <- as.polygons(lay, values=TRUE, na.rm=TRUE, digits=5)
        if (length(lay) > 0) {
          lay <- terra::subset(x = lay, subset = lay[[names(lay)]] == 1)
        }
      }
      
      # Modified buffer handling
      if (exists("skipBuffer") && skipBuffer) {
        bLay <- lay
      } else if (bufferSize == 0) {
        bLay <- lay
      } else {
        if (length(lay) > 0) {
          bLay <- tryCatch({
            terra::buffer(lay, width = bufferSize)
          }, warning = function(w) {
            warning("Buffer operation warning for ", INDEX, ": ", w$message)
            lay
          }, error = function(e) {
            warning("Buffer operation error for ", INDEX, ": ", e$message)
            lay
          })
        } else {
          bLay <- lay
        }
      }
      
      if (length(bLay) > 0) {
        bLay <- terra::aggregate(bLay, dissolve = TRUE)
        totalAreaKm2 <- terra::expanse(bLay, transform = FALSE, unit = "km")
        message(crayon::green(paste0("Total area disturbed for ", INDEX, " as vector: ",
                                     round(totalAreaKm2, digits = 3), " km2")))
        message(paste0("This represents ", round(100*(totalAreaKm2/studyAreaTotalArea), 4),
                       "% of total area."))
        
        bLaySF <- sf::st_as_sf(bLay)
        if (nrow(bLaySF) > 0) {
          bLayRas <- tryCatch({
            fasterize::fasterize(sf = bLaySF, raster = rasterToMatchR)
          }, error = function(e) {
            warning("fasterize failed: ", e$message)
            NULL
          })
          if (!is.null(bLayRas)) names(bLayRas) <- INDEX
          distTb <- table(bLayRas[])
          totDistKm <- ifelse("1" %in% names(distTb), distTb["1"] * 0.0625, 0)
          message(crayon::green(paste0("Total area disturbed for ", INDEX,
                                       " as raster: ", totDistKm, " km2")))
        } else {
          bLayRas <- NULL
        }
      } else {
        bLayRas <- NULL
      }
      
      if (convertToRaster) {
        return(bLayRas)
      } else {
        if (exists("bLay") && length(bLay) > 0) {
          return(bLay)
        } else {
          return(NULL)
        }
      }
    })
    
    # Filter out NULL results from allBuffered
    allBuffered <- allBuffered[!sapply(allBuffered, is.null)]
    
    if (convertToRaster) {
      if (length(allBuffered) == 0) {
        # Return a zero-filled raster
        zeroRaster <- rasterToMatch
        zeroRaster[] <- 0
        return(zeroRaster)
      }
      
      bufferedAnthropogenicDisturbance500m <- rasterToMatch
      bufferedAnthropogenicDisturbance500m[!is.na(bufferedAnthropogenicDisturbance500m[])] <- 0
      for (i in seq_along(allBuffered)) {
        if (!is.null(allBuffered[[i]])) {
          bufferedAnthropogenicDisturbance500m[which(allBuffered[[i]][] == 1)] <- 1
        }
      }
      if (inherits(bufferedAnthropogenicDisturbance500m, "RasterLayer")) {
        bufferedAnthropogenicDisturbance500m <- terra::rast(bufferedAnthropogenicDisturbance500m)
      }
      distTable <- table(bufferedAnthropogenicDisturbance500m[], useNA = "no")
      totDistKm <- ifelse("1" %in% names(distTable), distTable["1"]*0.0625, 0)
      percDist <- ifelse("0" %in% names(distTable) && "1" %in% names(distTable),
                         round(100*(distTable["1"]/(distTable["0"]+distTable["1"])), 2), 0)
      totalArea <- ifelse("0" %in% names(distTable) && "1" %in% names(distTable),
                          0.0625*(distTable["0"]+distTable["1"]), NA)
      message(paste0("Buffered disturbances for all disturbances. Total area disturbed: ",
                     totDistKm, "km2 -- ", percDist, "% of the total area (",
                     totalArea,"km2)"))
      return(bufferedAnthropogenicDisturbance500m)
    } else {
      if (length(allBuffered) > 0) {
        disturbanceLayerAll <- do.call(rbind, allBuffered)
        disturbanceLayerAll <- terra::aggregate(disturbanceLayerAll, dissolve = TRUE)
        totalDisturbedArea <- terra::expanse(disturbanceLayerAll, transform = FALSE, unit = "km")
        message(paste0("Buffered disturbances for all disturbances. Total area disturbed: ",
                       round(totalDisturbedArea, 2), " km2 -- ",
                       round(100*(totalDisturbedArea/studyAreaTotalArea), 2), "% of the total area (",
                       round(studyAreaTotalArea, 2)," km2)"))
        return(disturbanceLayerAll)
      } else {
        stop("No valid layers to combine in vector mode")
      }
    }
}
