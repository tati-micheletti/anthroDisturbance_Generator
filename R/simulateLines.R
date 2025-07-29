simulateLines <- function(Lines, distThreshold = 5000,
                         distNewLinesFact = 2,
                         refinedStructure = FALSE){

  buffDist <- distNewLinesFact*distThreshold # Distance buffer for new lines from a center point
  # THIS IS WHERE WE START SIMULATING THEM!
  for (potclus in unique(Lines$Pot_Clus)){
    
    # Getting the subset
    exSet <- Lines[Lines$Pot_Clus == potclus,]
    
    # Number of lines to simulate
    clCentr <- centroids(exSet)

    # Need to make a buffer around the centroid to get a random point for the center of the lines
    centBuff <- terra::buffer(x = clCentr, width = buffDist)
    centB <- terra::spatSample(centBuff, size = NROW(clCentr), 
                               method = "random") # Here we need to get a random point within this buffered area!
    cent <- terra::crds(centB)
    nLines <- NROW(exSet)

    if (length(nLines) > 1){
      print(paste0("length(nLines) is > 1. Maybe there is a problem? Investigate.",
                   "This should be one value per Cluster_Potential as total lines for that",
                   "cluster... Entering debug mode."))
      browser()
    }
    xlim <- c(min(cent[,"x"]), max(cent[,"x"]))
    ylim <- c(min(cent[,"y"]), max(cent[,"y"]))
    
    angles <- numeric(NROW(exSet))
    for (ii in 1:NROW(exSet)) {
      angles[ii] <- calculateLineAngle(exSet[ii, ])
    }
    
    if (nLines > 1){
      # Specify the number of parallel and perpendicular pairs
      # Identifying lines' configuration
      # Define tolerance for angles to consider lines parallel or perpendicular
      tolerance <- 5  # Degrees
      # Initialize counters
      parallel_count <- 0
      perpendicular_count <- 0

      if (length(angles)>1){
        # Determine the number of pairs
        num_pairs <- floor(length(angles) / 2)
        
        
        # Compare each pair of lines
        for (ii in seq(1, 2 * num_pairs, by = 2)) {
        # for (ii in 1:(length(angles) - 1)) {
          # for (j in (ii + 1):length(angles)) {
            # angle_diff <- abs(angles[ii] - angles[j])
              angle_diff <- abs(angles[ii] - angles[ii + 1])
            if (angle_diff <= tolerance || 
                abs(angle_diff - 180) <= tolerance ||
                abs(angle_diff - 365) <= tolerance) {
              parallel_count <- parallel_count + 1
            } else {
              perpendicular_count <- perpendicular_count + 1
            }
          # }
        }
      }

      SD <- sd(perim(exSet))
      # Catch error
      if (any(is.na(SD))){
        print(paste0("Standard deviation of lines for Pot_Clus ", potclus," is NA.",
                     "This should not happen. ",
                     "Entering debug mode."))
        browser()
      }

      n_parallel_pairs <- parallel_count
      n_perpendicular_pairs <- perpendicular_count
      lineLengths <- truncnorm::rtruncnorm(nLines,
                                           a = min(perim(exSet)),
                                           b = max(perim(exSet)),
                                           mean = mean(perim(exSet)),
                                           sd = if (!is.na(SD)) SD else 0)
    } else {
      lineLengths <- perim(exSet)
      if (length(lineLengths) > 1) stop("Something went wrong. Please debug.")
    }

    # Create an empty list to store the simulated lines
    simulatedLines <- vector("list", nLines)

    if (refinedStructure){
      # Initialize a vector to store the spacing
      spacing <- as.numeric(nrow(clCentr))
      
      # If more than one, then add the information about minimum distance, etc. 
      # If only one, no need.
      if (spacing > 1){
        for (ii in seq_len(nrow(clCentr)-1)) {
          # Exclude the current centroid
          nextCentroid <- clCentr[(ii+1), ]
          
          # Calculate the distance to the nearest other centroid
          dists <- distance(clCentr[ii, ], nextCentroid)
          spacing[ii] <- min(dists)
        }
      }
      if (nLines > 1){
        # Generate parallel pairs
        if (n_parallel_pairs > 0){
          for (k in 1:(n_parallel_pairs / 2)) {
            angle <- runif(1, 0, 180)
            ind1 <- k * 2 - 1
            ind2 <- k * 2
            simulatedLines[[ind1]] <- generateLine(angle = angle, length = lineLengths[ind1], 
                                                   xlim = xlim, ylim = ylim, mCrs = exSet)
            # Get a buffer around the first line
            buffLine1 <- terra::buffer(simulatedLines[[ind1]], width = spacing[k])
            # Get the boundary of the buffer (this will be exactly 'spacing' distance from the line)
            boundaryL1 <- as.points(buffLine1)
            # Get a random point on this boundary for the start point
            n_points <- nrow(boundaryL1)
            random_index <- sample(1:n_points, 1)
            newLineStart <- boundaryL1[random_index]
            # Update xlim/ylim for the generateLine function based on this point
            new_coords <- terra::crds(newLineStart)
            simulatedLines[[ind2]] <- generateLine(angle = angle + runif(1, -tolerance, tolerance), 
                                                   length = lineLengths[ind2], 
                                                   xlim = c(new_coords[1], new_coords[1]), 
                                                   ylim = c(new_coords[2], new_coords[2]), 
                                                   mCrs = exSet)
          }
        }
        # Generate perpendicular pairs
        if (n_perpendicular_pairs > 0){
          for (k in 1:(n_perpendicular_pairs / 2)) {
            angle <- runif(n = 1, min = 0, max = 180)
            ind1 <- n_parallel_pairs + k * 2 - 1
            ind2 <- n_parallel_pairs + k * 2
            simulatedLines[[ind1]] <- generateLine(angle, length = lineLengths[ind1],
                                                   xlim = xlim, ylim = ylim, mCrs = exSet)
            # Get a buffer around the first line
            buffLine1 <- terra::buffer(simulatedLines[[ind1]], width = spacing[k])
            # Get the boundary of the buffer (this will be exactly 'spacing' distance from the line)
            boundaryL1 <- as.points(buffLine1)
            
            # Get a random point on this boundary for the start point
            n_points <- nrow(boundaryL1)
            random_index <- sample(1:n_points, 1)
            newLineStart <- boundaryL1[random_index]
            # Update xlim/ylim for the generateLine function based on this point
            new_coords <- terra::crds(newLineStart)
            simulatedLines[[ind2]] <- generateLine(angle = angle + 90 + runif(1, -tolerance, tolerance), 
                                                   length = lineLengths[ind2], 
                                                   xlim = c(new_coords[1], new_coords[1]), 
                                                   ylim = c(new_coords[2], new_coords[2]), 
                                                   mCrs = exSet)
          }
        }
        # Generate remaining random lines
        if (2*(n_parallel_pairs + n_perpendicular_pairs) < nLines){ # Means we still need to simulate Lines
          for (k in sum(sapply(simulatedLines, is.null)):nLines) {
            angle <- runif(1, 0, 180)
            simulatedLines[[k]] <- generateLine(angle, length = lineLengths[k],
                                                xlim = xlim, ylim = ylim, mCrs = exSet)
          }
        } 
      } else {
        # Just one line
        angle <- runif(1, 0, 180)
        simulatedLines[[1]] <- generateLine(angle, length = lineLengths[1],
                                            xlim = xlim, ylim = ylim, mCrs = exSet)
      }
    } else {
      angle <- runif(nLines, 0, 180)
      simulatedLines <- lapply(seq_along(angles), function(indexAngle){
        if (is.na(lineLengths[indexAngle])) {
          print(paste0("lineLenghts index ",indexAngle," is NA (simulateLines.R). Debug."))
          browser()
        }
      generatedLines <-  generateLine(angle[indexAngle], length = lineLengths[indexAngle], 
                                      xlim = xlim, 
                                      ylim = ylim, 
                                      mCrs = exSet)
      return(generatedLines)
      })
    }
    
    # Combine the simulated lines into a single `SpatVector`
    sLines <- do.call(rbind, simulatedLines)
    sLines$Pot_Clus <- unique(exSet$Pot_Clus)
    sLines$lineLength <- perim(sLines)
    # if (exists("createdLines")){ # How can they exist?!?!
      # createdLines <- rbind(createdLines, sLines)
    # } else {
      createdLines <- sLines 
    # }
    # createdLines[createdLines$Pot_Clus == potclus, "Pot_Clus"] <- paste0(unique(exSet$Potential), "_",potclus)
    # createdLines[createdLines$Pot_Clus == potclus, "calculatedLength"] <- perim(createdLines[createdLines$Pot_Clus == potclus, ])
    # createdLines[createdLines$Pot_Clus == potclus, "angles"] <- angles
      createdLines[, "Pot_Clus"] <- paste0(unique(exSet$Potential), "_",potclus)
      createdLines[, "calculatedLength"] <- perim(createdLines[createdLines$Pot_Clus == potclus, ])
      createdLines[, "angles"] <- angles
      createdLines[, "Potential"] <- sub("_.*", "", createdLines$Pot_Clus)
  }
  return(createdLines)
}

