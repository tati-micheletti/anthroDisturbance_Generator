checkOverlap <- function(line1, line2) {
  # # Calculate the intersection
  # buff1 <- terra::buffer(line1, 10)
  # buff2 <- terra::buffer(line2, 10)
  # doIntersec <- terra::intersect(buff1, buff2)
  # if (length(doIntersec) == 0) return(FALSE) 
  # return(TRUE)
  overlap_matrix <- terra::relate(x = line1, y = line2, relation = "T********")
  overlap_matrix[1,2]
  
  }
