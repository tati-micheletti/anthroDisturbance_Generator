calculateLineAngle <- function(Line) {
  coords <- terra::crds(Line)  # Extract coordinates of the line
  delta_y <- coords[2, 2] - coords[1, 2]
  delta_x <- coords[2, 1] - coords[1, 1]
  angle <- terra::atan2(delta_y, delta_x) * 180 / pi  # Convert to degrees
  return(angle)
}
