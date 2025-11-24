generateLine <- function(angle, length, xlim, ylim, mCrs) {
  if (any(is.na(angle))) stop("angle is NA. Please debug")
  if (any(is.na(xlim))) stop("xlim is NA. Please debug")
  if (any(is.na(ylim))) stop("ylim is NA. Please debug")
  if (any(is.na(length))) stop("length is NA. Please debug")

  crs_val <- tryCatch(terra::crs(mCrs, proj = TRUE), error = function(e) "")
  if (is.null(crs_val) || !nzchar(crs_val)) stop("mCrs is NA or missing. Please debug")
  
  # additional safeguards
  if (!is.numeric(angle) || length(angle) != 1) {
    stop("angle must be a single numeric value")
  }
  if (!is.numeric(length) || length(length) != 1 || length < 0) {
    stop("length must be a single non-negative numeric")
  }
  if (!is.numeric(xlim) || length(xlim) != 2) {
    stop("xlim must be a numeric vector of length 2")
  }
  if (!is.numeric(ylim) || length(ylim) != 2) {
    stop("ylim must be a numeric vector of length 2")
  }
  
  # Generate a random starting point within the limits
  start_x <- runif(1, xlim[1], xlim[2])
  start_y <- runif(1, ylim[1], ylim[2])
  
  # Calculate the endpoint based on the angle
  end_x <- as.numeric(start_x + length * cos(angle * pi / 180))
  end_y <- as.numeric(start_y + length * sin(angle * pi / 180))
  
  
  # Return the line as a `SpatVector` line feature
  coords_mat <- matrix(
    c(start_x, start_y,
      end_x,   end_y),
    ncol = 2,
    byrow = TRUE
  )
  pts <- terra::vect(coords_mat, type = "points", crs = crs_val)
  lin <- terra::as.lines(pts)
  lin$Class  <- "Seismic"
  
  return(lin)
}
