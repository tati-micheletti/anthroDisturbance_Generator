simulateLines <- function(Lines, distThreshold = 5000,
                              distNewLinesFact = 2,
                              refinedStructure = FALSE) {
  # Input validation at the start
  if (nrow(Lines) == 0) stop("simulateLines: no input features")
  if (!"Pot_Clus" %in% names(Lines)) stop("simulateLines: missing Pot_Clus field")
  
  # Initialize createdLines properly
  createdLines <- Lines[0, ]  # Empty SpatVector with same structure as input
  
  buffDist <- distNewLinesFact * distThreshold # Distance buffer for new lines from a center point
  
  # Debug output
  #cat("Starting simulation with", nrow(Lines), "input lines\n")
  
  for (potclus in unique(Lines$Pot_Clus)) {
    #cat("Processing cluster:", potclus, "\n")
    
    # Getting the subset
    exSet <- Lines[Lines$Pot_Clus == potclus,]
    
    # Number of lines to simulate
    clCentr <- centroids(exSet)
    # Need to make a buffer around the centroid to get a random point for the center of the lines
    centBuff <- terra::buffer(x = clCentr, width = buffDist)
    centB <- terra::spatSample(centBuff, size = NROW(clCentr),
                               method = "random") # Here we need to get a random point within this buffered area!
    cent <- crds(centB)
    nLines <- NROW(exSet)
    xlim <- c(min(cent[,"x"]), max(cent[,"x"]))
    ylim <- c(min(cent[,"y"]), max(cent[,"y"]))
    
    angles <- numeric(NROW(exSet))
    for (ii in 1:NROW(exSet)) {
      angles[ii] <- calculateLineAngle(exSet[ii, ])
    }
    
    if (nLines > 1) {
      # Specify the number of parallel and perpendicular pairs
      # Identifying lines' configuration
      # Define tolerance for angles to consider lines parallel or perpendicular
      tolerance <- 5  # Degrees
      # Initialize counters
      parallel_count <- 0
      perpendicular_count <- 0
      if (length(angles) > 1) {
        # Determine the number of pairs
        num_pairs <- floor(length(angles) / 2)
        
        # Compare each pair of lines
        for (ii in seq(1, 2 * num_pairs, by = 2)) {
          angle_diff <- abs(angles[ii] - angles[ii + 1])
          if (angle_diff <= tolerance || abs(angle_diff - 180) <= tolerance) {
            parallel_count <- parallel_count + 1
          } else {
            perpendicular_count <- perpendicular_count + 1
          }
        }
      }
      SD <- sd(lengths(exSet))
      n_parallel_pairs <- parallel_count
      n_perpendicular_pairs <- perpendicular_count
      lineLengths <- truncnorm::rtruncnorm(nLines,
                                           a = min(lengths(exSet)),
                                           b = max(lengths(exSet)),
                                           mean = mean(lengths(exSet)),
                                           sd = if (!is.na(SD)) SD else 0)
    } else {
      lineLengths <- lengths(exSet)
      if (length(lineLengths) > 1) stop("Something went wrong. Please debug.")
    }
    # Initialize with random lines for all
    simulatedLines <- lapply(1:nLines, function(i) {
      generateLine(angle = runif(1, 0, 180),
                   length = lineLengths[i],
                   xlim = xlim,
                   ylim = ylim,
                   mCrs = exSet)
    })
    
    if (refinedStructure) {
      # Initialize a vector to store the spacing
      spacing <- rep(NA_real_, nrow(clCentr))
      
      # If more than one, then add the information about minimum distance, etc.
      # If only one, no need.
      if (nrow(clCentr) > 1) {
        for (ii in seq_len(nrow(clCentr)-1)) {
          # Exclude the current centroid
          nextCentroid <- clCentr[(ii+1), ]
          
          # Calculate the distance to the nearest other centroid
          dists <- distance(clCentr[ii, ], nextCentroid)
          spacing[ii] <- min(dists)
        }
        # For the last centroid, use median spacing or some reasonable value
        spacing[nrow(clCentr)] <- median(spacing[1:(nrow(clCentr)-1)], na.rm = TRUE)
      } else {
        # For single centroid, use default spacing
        spacing <- buffDist / 2
      }
      
      spacing <- pmax(1, spacing)
      
      if (any(is.na(spacing))) spacing[is.na(spacing)] <- buffDist / 2
      
      # Generate parallel pairs with bounds checking
      if (n_parallel_pairs > 0) {
        for (k in 1:n_parallel_pairs) {
          if (2*k > nLines) break
          angle <- runif(1, 0, 180)
          simulatedLines[[2*k-1]] <- generateLine(angle, lineLengths[2*k-1], xlim, ylim, exSet)
          
          # Generate nearby parallel line
          buff_width <- if (k <= length(spacing)) max(1, abs(spacing[k])) else buffDist/2
          buffLine <- terra::buffer(simulatedLines[[2*k-1]], width = buff_width)
          boundary_pts <- as.points(buffLine)
          if (nrow(boundary_pts) > 0) {
            pt <- boundary_pts[sample(nrow(boundary_pts)), 1]
            pt_coords <- crds(pt)
            simulatedLines[[2*k]] <- generateLine(
              angle + runif(1, -tolerance, tolerance),
              lineLengths[2*k],
              xlim = c(pt_coords[1], pt_coords[1]),
              ylim = c(pt_coords[2], pt_coords[2]),
              mCrs = exSet
            )
          }
        }
      }
      
      # Generate perpendicular pairs with proper indexing
      if (n_perpendicular_pairs > 0) {
        base_index <- 2 * n_parallel_pairs
        for (k in 1:n_perpendicular_pairs) {
          idx <- base_index + 2*k - 1
          if (idx + 1 > nLines) break
          
          angle <- runif(1, 0, 180)
          simulatedLines[[idx]] <- generateLine(angle, lineLengths[idx], xlim, ylim, exSet)
          
          # Generate nearby perpendicular line
          buff_width <- if (k <= length(spacing)) max(1, abs(spacing[k])) else buffDist/2
          buffLine <- terra::buffer(simulatedLines[[idx]], width = buff_width)
          boundary_pts <- as.points(buffLine)
          if (nrow(boundary_pts) > 0) {
            pt <- boundary_pts[sample(nrow(boundary_pts)), 1]
            pt_coords <- crds(pt)
            simulatedLines[[idx+1]] <- generateLine(
              angle + 90 + runif(1, -tolerance, tolerance),
              lineLengths[idx+1],
              xlim = c(pt_coords[1], pt_coords[1]),
              ylim = c(pt_coords[2], pt_coords[2]),
              mCrs = exSet
            )
          }
        }
      }
    }
    
    # Safely combine lines and handle degenerates
    validLines <- simulatedLines[!sapply(simulatedLines, is.null)]
    sLines <- do.call(rbind, simulatedLines[!sapply(simulatedLines, is.null)])
    if (nrow(sLines) > 0) {
      sLines$Pot_Clus <- potclus
      sLines$lineLength <- lengths(sLines)
      sLines$angles <- sapply(1:nrow(sLines), function(i) calculateLineAngle(sLines[i, ]))
      createdLines <- rbind(createdLines, sLines)
    }
  }
  
  createdLines$calculatedLength <- lengths(createdLines)
  createdLines$Potential      <- sub("_.*", "", createdLines$Pot_Clus)
  return(createdLines)
}
