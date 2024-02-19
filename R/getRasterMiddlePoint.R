getRasterMiddlePoint <- function(ras){
  row <- trunc(terra::nrow(ras)/2)
  col <- trunc(terra::ncol(ras)/2)
  cell <- terra::cellFromRowCol(ras, row, col)
  return(cell) 
}
