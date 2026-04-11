cleanupList <- function(aList, 
                        outter = FALSE, 
                        inner = TRUE, 
                        cleanEmpty = FALSE, 
                        nullEmpty = TRUE) {
  # Step 1: Deep copy input once (protect spatial data integrity)
  tempList <- copy(aList)
  
  # Step 2: Inner-level NULL removal (preserve NULL elements)
  if (inner) {
    tempList <- lapply(tempList, function(innerList) {
      if (is.null(innerList)) return(NULL)
      if (!is.list(innerList)) return(innerList)
      Filter(Negate(is.null), innerList)
    })
  }
  
  # Step 3: Top-level NULL removal
  if (outter) {
    tempList <- Filter(Negate(is.null), tempList)
  }
  
  # Step 4: Convert empty lists to NULL
  if (nullEmpty) {
    tempList <- lapply(tempList, function(x) {
      if (is.list(x) && length(x) == 0) {
        NULL
      } else {
        x
      }
    })
  }
  
  # Step 5: Remove empty-list and NULL elements if requested
  if (cleanEmpty) {
    tempList <- Filter(function(x) {
      !((is.list(x) && length(x) == 0) || is.null(x))
    }, tempList)
  }
  
  # Final: Deep copy before returning (protect modified list)
  result <- copy(tempList)
  return(result)
}
