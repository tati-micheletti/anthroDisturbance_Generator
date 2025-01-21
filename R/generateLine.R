generateLine <- function(angle, length, xlim, ylim, mCrs) {
  
  if (any(is.na(angle))) stop("angle is NA. Please debug")
  if (any(is.na(xlim))) stop("xlim is NA. Please debug")
  if (any(is.na(ylim))) stop("ylim is NA. Please debug")
  if (any(is.na(mCrs))) stop("mCrs is NA. Please debug")
  if (any(is.na(length))) stop("length is NA. Please debug")
  
  # Generate a random starting point within the limits
  start_x <- runif(1, xlim[1], xlim[2])
  start_y <- runif(1, ylim[1], ylim[2])
  
  # Calculate the endpoint based on the angle
  end_x <- as.numeric(start_x + length * cos(angle * pi / 180))
  end_y <- as.numeric(start_y + length * sin(angle * pi / 180))
  
  
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
