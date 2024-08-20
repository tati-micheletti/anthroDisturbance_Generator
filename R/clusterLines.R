# First, per polygon (I have the area divided into polygons of oil potential), 
# I cluster the Lines to identify “individual grids” (or as close to it as possible). 
# Once I know the clusters from data, I can identify the average 
# 1) number of Lines, as well 
# 2) line distance, 
# 3) length, 
# 4) angles --> how many parallel and how many perpendicular Lines 
# respective SDs per cluster. 
# Then I can just randomly choose which type of cluster the generated next will be and draw the grid based on the data.
 
clusterLines <- function(Lines, distThreshold = 5000, plotClusterDiagnostic = FALSE){
  # Calculate the centroids of each line
  centroids <- centroids(Lines)
  
  xy <- geom(centroids)
  
  D <- dist(data.frame(rownames = centroids$individualID, 
                       x = xy[,"x"],
                       y = xy[,"y"]))
  chc <- hclust(D, method="complete")
  
  # Distance with a 5km threshold as default  
  chc.d <- cutree(chc, h=distThreshold)
  if (plotClusterDiagnostic){
    tb <- as.data.table(table(chc.d))
    plot(tb$chc.d, tb$N, xlab = "Number of Clusters", ylab = "Number of Lines")
  }
  Lines$cluster <- chc.d
  ## Plotting:
  # LinesSF <- st_as_sf(Lines)
  # LinesSFc <- LinesSF[, "cluster"]
  # plot(LinesSFc, col = rainbow(length(unique(LinesSFc$cluster))))
  for (i in unique(Lines$cluster)){
    Lines[Lines$cluster == i, "Potential_Cluster"] <- paste0(unique(Lines[Lines$cluster == i,]$Potential),
                                                         "_",i)
    Lines[Lines$cluster == i, "lengthMean"] <- mean(perim(Lines[Lines$cluster == i, ]))
    SD <- sd(perim(Lines[Lines$cluster == i, ]))
    Lines[Lines$cluster == i, "lengthSd"] <- if (!is.na(SD)) SD else 0
    Lines[Lines$cluster == i, "lengthMin"] <- min(perim(Lines[Lines$cluster == i, ]))
    Lines[Lines$cluster == i, "lengthMax"] <- max(perim(Lines[Lines$cluster == i, ]))
    Lines[Lines$cluster == i, "totalLines"] <- nrow(Lines[Lines$cluster == i, ])
    # AND SPACING!
    clCentr <- centroids(Lines[Lines$cluster == i, ])
    # Initialize a vector to store the spacing
    spacing <- numeric(nrow(clCentr))
    # Loop through each centroid to find the minimum distance to other centroids
    for (i in seq_len(nrow(clCentr))) {
      # Exclude the current centroid
      other_centroids <- clCentr[-i, ]
      
      # Calculate the distance to the nearest other centroid
      dists <- distance(clCentr[i, ], other_centroids)
      spacing[i] <- min(dists)
    }
    # Add spacing as a new column to the lines data
    Lines[Lines$cluster == i, "meanSpacing"] <- mean(spacing)
    SDsp <- sd(spacing)
    Lines[Lines$cluster == i, "sdSpacing"] <- if (!is.na(SDsp)) SDsp else 0
    Lines[Lines$cluster == i, "minSpacing"] <- min(spacing)
    Lines[Lines$cluster == i, "maxSpacing"] <- max(spacing)
    angles <- numeric(NROW(Lines[Lines$cluster == i, ]))
    for (ii in 1:NROW(Lines[Lines$cluster == i, ])){
      angles[ii] <- calculateLineAngle(Lines[Lines$cluster == i, ][ii,])
    }
    Lines[Lines$cluster == i, "angles"] <- angles
    # Identifying lines' configuration
    # Define tolerance for angles to consider lines parallel or perpendicular
    tolerance <- 5  # Degrees
    # Initialize counters
    parallel_count <- 0
    perpendicular_count <- 0
    
    # Compare each pair of lines
    for (i in 1:(length(angles) - 1)) {
      for (j in (i + 1):length(angles)) {
        angle_diff <- abs(angles[i] - angles[j])
        if (angle_diff <= tolerance || abs(angle_diff - 180) <= tolerance) {
          parallel_count <- parallel_count + 1
        } else if (abs(angle_diff - 90) <= tolerance) {
          perpendicular_count <- perpendicular_count + 1
        }
      }
    }
    Lines[Lines$cluster == i, "noParallelLines"] <- parallel_count
    Lines[Lines$cluster == i, "noPerpendicularLines"] <- perpendicular_count

    # TEST APPROACH!!
    # Number of lines to simulate
    exSet <- Lines[Lines$cluster == i]
    cent <- crds(clCentr)

    nLines <- as.numeric(unique(as.data.frame(exSet[,"totalLines"])))
    # Specify the number of parallel and perpendicular pairs
    n_parallel_pairs <- parallel_count
    n_perpendicular_pairs <- perpendicular_count
    # Define the area where lines will be placed
    xlim <- c(min(cent[,"x"]), max(cent[,"x"]))
    ylim <- c(min(cent[,"y"]), max(cent[,"y"]))
    
    # HERE!!!!!! From here down, test!!!
    
    # THIS BELOW IS THE CORRECT WAY OF GETTING LINE SIZES AND SPACING!
    lineLengths <- truncnorm::rtruncnorm(nLines,
                          a = unique(as.data.frame(exSet[, "lengthMin"])),
                          b = unique(as.data.frame(exSet[, "lengthMax"])),
                          mean = unique(as.data.frame(exSet[, "lengthMean"])),
                          sd = unique(as.data.frame(exSet[, "lengthSd"])))
    

    browse() #<~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ I AM HERE
    
    # HERE!!! I still need to control for the distance (spacing) between the lines
    # It is something I probably need to do at the Line 145 using the info meanSpacing...
    
    generateLine <- function(angle, length, xlim, ylim, mCrs) {
      # Generate a random starting point within the limits
      start_x <- runif(1, xlim[1], xlim[2])
      start_y <- runif(1, ylim[1], ylim[2])
      
      # Calculate the endpoint based on the angle
      end_x <- start_x + length * cos(angle * pi / 180)
      end_y <- start_y + length * sin(angle * pi / 180)
      
      # Return the line as a `SpatVector` line feature
      point1 <- as.points(terra::ext(cbind(c(start_x, start_x), c(start_y, start_y))),
                          crs = terra::crs(mCrs))
      point2 <- as.points(terra::ext(cbind(c(end_x, end_x), c(end_y, end_y))),
                          crs = terra::crs(mCrs))
      bothPoints <- rbind(point1, point2)
      lin <- terra::as.lines(bothPoints)
      lin$Class  <- "Seismic"
      return(lin)
    }
    
    # Create an empty list to store the simulated lines
    simulatedLines <- vector("list", nLines)
    
    # Generate parallel pairs
    for (i in 1:(n_parallel_pairs / 2)) {
      angle <- runif(1, 0, 180)
      ind1 <- i * 2 - 1
      ind2 <- i * 2
      simulatedLines[[ind1]] <- generateLine(angle = angle, length = lineLengths[ind1], 
                                                   xlim = xlim, ylim = ylim, mCrs = exSet)
      simulatedLines[[ind2]] <- generateLine(angle = angle + runif(1, -tolerance, tolerance), 
                                               length = lineLengths[ind2], 
                                               xlim = xlim, ylim = ylim, mCrs = exSet)
    }
    
    # Generate perpendicular pairs
    for (i in 1:(n_perpendicular_pairs / 2)) {
      angle <- runif(1, 0, 180)
      ind1 <- n_parallel_pairs + i * 2 - 1
      ind2 <- n_parallel_pairs + i * 2
      simulatedLines[[ind1]] <- generateLine(angle, length = lineLengths[ind1],
                                              xlim = xlim, ylim = ylim, mCrs = exSet)
      simulatedLines[[ind2]] <- generateLine(angle + 90, length = lineLengths[ind2],
                                          xlim = xlim, ylim = ylim, mCrs = exSet)
    }
    
    # Generate remaining random lines
    for (i in (n_parallel_pairs + n_perpendicular_pairs + 1):nLines) {
      angle <- runif(1, 0, 180)
      simulatedLines[[i]] <- generateLine(angle, length = lineLengths[i],
                                          xlim = xlim, ylim = ylim, mCrs = exSet)
    }
    
    # Combine the simulated lines into a single `SpatVector`
    sLines <- do.call(rbind, simulatedLines)
    
    # NOW Remember to add to this new cluster the same information of the others (class at least? 
    # I think the whole cluster info will be redone in the next step again...)
    


  }

  # Add line width per cluster
  

  browser()
  # merge(p, dfr, all.x=TRUE, by.x=c('NAME_1', 'NAME_2'), by.y=c('District', 'Canton'))
  return(Lines)
}
