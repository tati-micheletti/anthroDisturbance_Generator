seismicGrids <- function(Line,
                         centerPoint,
                         originalLayer,
                         linesLength,
                         randomLinesN = runif(n = 1, min = 0, max = 4),
                         fraction = 0.002, # Percentage of the distanceBetweenLines for the jittering
                         distanceBetweenLines = rtnorm(n = 1,
                                                       mean = mean(c(50, 100)),
                                                       sd = sd(c(50, 100)),
                                                       lower = 50),
                         howMany,
                         existingLine) {
  
  # Unfortunately, we can't get groups of lines from data so I 
  # can know how many lines I should have. Therefore, we are 
  # gonna use this howMany instead.
  
  # Function based on the code available at 
  # https://rdrr.io/github/RodrigoAgronomia/PAR/src/R/trial_geom.R 
  
  calc_angle <- function(x, y) {
    dst_diff <- as.numeric(x - y)
    return(atan2(dst_diff[1], dst_diff[2]) + pi)
  }
  
  if (!any(is(Line, "SpatVector"),
           is(Line, "sf"))) stop("Line musst be sf or SpatVector") 
  if (!is(centerPoint, "SpatVector")) stop("centerPoint musst be a SpatVector") 
  
  if (length(distanceBetweenLines) == 1)
    distanceBetweenLines <- rep(distanceBetweenLines, times = 2)
  
  if (is(Line, "SpatVector")){
    LineSF <- sf::st_as_sf(Line)
  } else {
    LineSF <- Line
  }

  cl <- sf::st_coordinates(LineSF)[, 1:2]
  angle <- calc_angle(cl[1, ], cl[nrow(cl), ]) + pi / 2
  dist2 <- distance(Line, centerPoint)
  no <- ceiling(dist2 / distanceBetweenLines[1])
  startOffset <- length(1:(no-howMany[1]-1))
  if (existingLine){
    offset_path <- list(prll = c(rep(0, times = startOffset), 
                                 runif(n = howMany[1]+1, # The plus 1 is due to the original line being included 
                                       min = -fraction*distanceBetweenLines[1], 
                                       max = fraction*distanceBetweenLines[1])))
  } else {
    offset_path <- list(prll = c(runif(n = howMany[1]+1, # The plus 1 is due to the original line being included 
                                       min = -fraction*distanceBetweenLines[1], 
                                       max = fraction*distanceBetweenLines[1]),
                                 rep(0, times = startOffset)))
  }
  # angleJitter <- data.frame(runif(no, 0, 3), 
  # # Didn't work as I expected. I thought it would change lines angles, but just displaces them 
  # over their own axis
  nl <- sapply(seq(0, no), simplify = FALSE, function(w) {
    off <- if (w == 0) 0 else offset_path[["prll"]][w]
    # jitSin <- if (w == 0) 0 else angleJitter[w, 1]
    # jitCos <- if (w == 0) 0 else angleJitter[w, 2]
    lineTempl <- LineSF
    # I can use linesLength to fiddle with the tail/head of the lines. 
    # First I check the line size, 
    LineSFSize <- perim(Line)
    # then I get a new size. 
    LineSFNewSize <- round(eval(parse(text = linesLength)), 0)
    # Reduce from the original 
    LineSFSizeDiff <- LineSFNewSize-LineSFSize
    # randomly split the leftover between head and tail. TADA!
    spl <- sample(x = seq(0, LineSFSizeDiff, by = 0.01), size =  1)
    leadLineSplit <- c(spl, LineSFSizeDiff-spl)
    
    sf::st_geometry(lineTempl) <- st_extend_line(sf::st_geometry(lineTempl), 
                                                 distance = c(leadLineSplit[1], leadLineSplit[2]),
                                                 end = "BOTH")
    sf::st_geometry(lineTempl) + (off + w) * distanceBetweenLines[1] * c(sin(angle), cos(angle))
  })
  path_lines <- sf::st_sfc(do.call(rbind, nl), crs = sf::st_crs(Line))
  if (existingLine){
    path_lines <- vect(path_lines[(NROW(path_lines)-howMany[1]):NROW(path_lines)])
  } else {
    path_lines <- vect(path_lines[1:howMany[1]])
  }
  # terra::plot(path_lines, col = "forestgreen")
  
  # Now the crossing lines
  
  # Find the center of the first line
  x1 <- (xmin(path_lines[1, ])+xmax(path_lines[1, ]))/2
  y1 <- (ymin(path_lines[1, ])+ymax(path_lines[1, ]))/2
  # Do the same for the last line
  x2 <- (xmin(path_lines[nrow(path_lines), ])+xmax(path_lines[nrow(path_lines), ]))/2
  y2 <- (ymin(path_lines[nrow(path_lines), ])+ymax(path_lines[nrow(path_lines), ]))/2
  # Create the points
  centPoint1 <- vect(matrix(c(x1, y1), ncol = 2), type="points", atts=NULL, crs=crs(path_lines))
  centPoint2 <- vect(matrix(c(x2, y2), ncol = 2), type="points", atts=NULL, crs=crs(path_lines))
  # terra::plot(rbind(centPoint1, centPoint2), add = TRUE, col = "red")
  leadLine <- as.lines(x = rbind(centPoint1, centPoint2))
  cl <- sf::st_coordinates(sf::st_as_sf(leadLine))[, 1:2]
  angle <- calc_angle(cl[1, ], cl[nrow(cl), ]) + pi / 2
  no2 <- ceiling(howMany[2]/2)
  offset_path[["orth"]] <- runif(n = length(seq(-no2, no2)),
                                 min = -fraction*distanceBetweenLines[2],
                                 max = fraction*distanceBetweenLines[2])
  nl2 <- sapply(seq(-no2, no2), simplify = FALSE, function(w) {
    off <- offset_path[["orth"]][w+no2+1]
    lineTempl <- st_as_sf(leadLine)
    # I can use linesLength to fiddle with the tail/head of the lines. 
    # First I check the line size, 
    leadLineSize <- perim(leadLine)
    # then I get a new size. 
    leadLineNewSize <- round(eval(parse(text = linesLength)), 0)
    # Reduce from the original 
    leadLineDiff <- leadLineNewSize-leadLineSize
    # randomly split the leftover between head and tail. TADA!
    spl <- sample(x = seq(0, leadLineDiff, by = 0.01), size =  1)
    leadLineSplit <- c(spl, leadLineDiff-spl)
    sf::st_geometry(lineTempl) <- st_extend_line(sf::st_geometry(lineTempl), 
                                                 distance = c(leadLineSplit[1], leadLineSplit[2]), 
                                                 end = "BOTH")
    sf::st_geometry(lineTempl) + (off + w) * distanceBetweenLines[2] * c(sin(angle), cos(angle))
  })
  path_lines2 <- vect(sf::st_sfc(do.call(rbind, nl2), crs = sf::st_crs(leadLine)))
  terra::plot(path_lines2, add = TRUE, col = "blue")
  
  # Now throw X random lines in there
  browser()
  randomLinesN
  path_rand
  # round(eval(parse(text = dParOri[["disturbanceSize"]])), 0)
  
  gridLines <-  rbind(path_lines, path_lines2, path_rand)
  return()
}

# Functions below copied/modified from 
# from https://stackoverflow.com/questions/65928126/increase-polyline-length-r-sf
st_ends_heading <- function(line){
  M <- sf::st_coordinates(line)
  i <- c(2, nrow(M) - 1)
  j <- c(1, -1)
  
  headings <- mapply(i, j, FUN = function(i, j) {
    Ax <- M[i-j,1]
    Ay <- M[i-j,2]
    Bx <- M[i,1]
    By <- M[i,2]
    unname(atan2(Ay-By, Ax-Bx))
  })
  
  return(headings)
}
st_extend_line <- function(line, distance, end = "BOTH"){
  if (!(end %in% c("BOTH", "HEAD", "TAIL")) | length(end) != 1) stop("'end' must be 'BOTH', 'HEAD' or 'TAIL'")
  
  M <- sf::st_coordinates(line)[,1:2]
  keep <- !(end == c("TAIL", "HEAD"))
  
  ends <- c(1, nrow(M))[keep]
  headings <- st_ends_heading(line)[keep]
  distances <- if (length(distance) == 1) rep(distance, 2) else rev(distance[1:2])
  
  M[ends,] <- M[ends,] + distances[keep] * c(cos(headings), sin(headings))
  newline <- sf::st_linestring(M)
  
  # If input is sfc_LINESTRING and not sfg_LINESTRING
  if (is.list(line)) newline <- sf::st_sfc(newline, crs = sf::st_crs(line))
  
  return(newline)
}
