replaceListFast <- function(disturbanceList, 
                        updatedLayers, 
                        currentTime,
                        disturbanceParameters){
  
  correctMatching <- match(names(updatedLayers), names(disturbanceList))
  
  externalLays <- lapply(unique(sort(correctMatching)), function(INDEX) {
    Sector <- names(disturbanceList)[INDEX]
    updatedLayerIndex <- which(correctMatching==INDEX)
    # Each sector might have more than one disturbance type. We don't need to 
    # aggregate them, just rbind. We can put a time marker on the past (if they don't have) 
    # and new layers: createdInSimulationTime
    # Get all potential names. Will be used internally and externally in the lapply
    potentialLayName <- names(disturbanceList[[Sector]])[grepl(x = names(disturbanceList[[Sector]]), 
                                                               pattern = "potential")]
    # Layers available for updates
    upToDate <- updatedLayers[updatedLayerIndex] 
    # As I have repeated names (i.e., for 
    # energy I have 2 layers named Energy, one for windTurbines and one for powerLines), I need a 
    # more specific approach to select the layers. Just referencing to the Sector leaves the second
    # one out.
    if (is.null(unlist(upToDate, use.names = FALSE))){
      # layName <- names(upToDate)
      # # Currently disturbed NULL
      # unified <- disturbanceList[[Sector]][[layName]]
      return(NULL) #TODO NEEDS TESTS!
    } else {
      names(upToDate) <- NULL # We need to remove the SECTOR names and unlist in case these two lists repeat the names.
      upToDate <- unlist(upToDate)
      # Merge the updated lists with the current disturbance ones:
      if (is.null(names(upToDate))){
        message(paste0("names(upToDate) for SECTOR = ",Sector," is NULL. Entering browser mode..."))
        browser()
      }
      internalLays <- lapply(names(upToDate), function(layName){
        # Updated disturbances
        currDist <- upToDate[[layName]]
        if (all(!is.null(currDist),
                !"createdInSimulationTime" %in% names(currDist),
                is(currDist, "SpatVector")))
          currDist[["createdInSimulationTime"]] <- currentTime        
        # The only problem is then the ones that are enlarging. We need to pass the table here too to detect
        # which ones need to just be replaced. In the current version of the module, this is the case of
        # settlements. So here we would look for any "enlarged and just replace if it is" 
        # If an enlarging, replace and move on 
        shouldReplace <- disturbanceParameters[dataName == Sector &
                                                 disturbanceOrigin == layName &
                                                 disturbanceEnd == "",
                                               disturbanceType] == "Enlarging"
        onlyReplace <- isTRUE(shouldReplace)
        if (onlyReplace){
          message(paste0("Disturbance layer ", layName," for ", Sector, " is of type Enlarging. ",
                         "This means it contains previous disturbances and can just replace the",
                         " previous layer."))
          return(currDist)
        }
        # Currently disturbed 
        pastDist <- disturbanceList[[INDEX]][[layName]]
        if (all(!is.null(pastDist),
                !"createdInSimulationTime" %in% names(pastDist),
                is(pastDist, "SpatVector")))
          pastDist[["createdInSimulationTime"]] <- currentTime - 10
        classPast <- class(pastDist)
        classCurr <- class(currDist)
        if (any(is.null(pastDist),
                is.null(currDist))) { # If one is missing, just return the other
          unified <- c(pastDist, currDist)[[1]]
        } else { # Need to merge them!
          if (any(classPast != classCurr,
                  classPast %in% c("RasterLayer", "SpatRaster"),
                  classCurr %in% c("RasterLayer", "SpatRaster"))){ # Objects are of different classes
            # This is likely for generated stuff (raster) from polygons of potential
            print("Objects are of different classes. Harmonizing...")
            if (any(classPast %in% c("RasterLayer", "SpatRaster"),
                    classCurr %in% c("RasterLayer", "SpatRaster"))) {
              # [UPDATE] Althought it is considerably quicker to convert polys to rasters than the other way around
              # I don't expect this to happen with the provided dataset. So here is only for new datasets used.
              # Convert the spatRaster into spatVector
              if (classPast == "SpatRaster"){
                pastDist <- as.polygons(pastDist)
                pastDist[["createdInSimulationTime"]] <- currentTime - 10
              } else {
                currDist <- as.polygons(currDist)
                currDist[["createdInSimulationTime"]] <- currentTime
              }
              unified <- rbind(pastDist, currDist)
            } else {
              stop(paste0("New case of combination of past and present layers. ",
                          "Past layer of class ", class(pastDist), 
                          " while new layer of class ", class(currDist),
                          ". Code development is needed!"))
            }
          } else { 
            # Objects share the same class and vectors
              if (geomtype(pastDist) != geomtype(currDist)) { # Works only if both are vectors!
                # Make sure that all are buffered, even if just a little 
                wid <- 0.0000000003
                if (geomtype(pastDist) != "polygons"){
                  pastDistBuf <- terra::buffer(pastDist, width = wid)
                  while (any(is.na(as.vector(ext(pastDistBuf))))){
                    message(paste0("Minimum buffering size current disturbances (", wid,
                                   ") failed. Increasing buffering size..."))
                    wid <- wid + wid
                    pastDistBuf <- terra::buffer(pastDist, width = wid)
                  }
                  pastDist <- pastDistBuf
                }
                if (geomtype(currDist) != "polygons"){
                  currDistBuf <- terra::buffer(currDist, width = wid)
                  while (any(is.na(as.vector(ext(currDistBuf))))){
                    message(paste0("Minimum buffering size for past disturbances (", wid,
                                   ") failed. Increasing buffering size..."))
                    wid <- wid + wid
                    currDistBuf <- terra::buffer(currDist, width = wid)
                  }
                  currDist <- currDistBuf
                }
              }
                  tictoc::tic(paste0("Elapsed time for merging ", layName, " (", Sector, ")"))
                  unified <- rbind(pastDist, currDist) # Doubles the features when we need to overlay/unify these
                  # But that is not a problem per se. In fact, we might want to see all geometries!
                  if (!"Class" %in% names(currDist)){
                    if ("Class" %in% names(pastDist)){
                      unified[["Class"]] <- as.character(unique(pastDist[["Class"]]))
                    } else {
                      unified[["Class"]] <- layName
                    }
                  }
                  tictoc::toc()
          }
        }
        return(unified)
      })
    }
    # Need to unlist whatever is listed here
    if (length(potentialLayName) > 0) {
      message(paste0("Appending potantial layer back to '", Sector, "' sector."))
      # If we have potential Layers, put them in the list with the internal layers
      internalLays <- append(disturbanceList[[Sector]][potentialLayName], internalLays)
    }
    internalLays <- unlist(internalLays, use.names = FALSE)
    names(internalLays) <- c(potentialLayName, names(upToDate))
    return(internalLays)
  })
  names(externalLays) <- names(disturbanceList)
  return(externalLays)
}