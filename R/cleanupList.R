cleanupList <- function(aList, 
                        outter = FALSE, 
                        inner = TRUE, 
                        cleanEmpty = FALSE, 
                        nullEmpty = TRUE){
  if (inner){
    cleanedList1 <- lapply(aList, function(innerList) {
      filtered_list <- Filter(Negate(is.null), innerList)
      if (length(filtered_list) == 0) {
        NULL # Return NULL if the list is empty after filtering
      } else {
        filtered_list # Otherwise, return the filtered list
      }
    })
  } else {
    cleanedList1 <- copy(aList)
  }
  if (outter){
    cleanedList2 <- Filter(Negate(is.null), cleanedList1)
  } else {
    cleanedList2 <- copy(cleanedList1)
  }
  
  if (nullEmpty){
    cleanedListFinal1 <- lapply(cleanedList2, function(x) {
      if (is.list(x) && length(x) == 0) {
        # If it is, return NULL
        return(NULL)
        } else {
        # Otherwise, return the element unchanged
        return(x)
      }
    })  } else {
    cleanedListFinal1 <- copy(cleanedList2)
  }
  
  if (cleanEmpty){
    cleanedListFinal <- Filter(function(x) {
      # Keep the element if it's NOT (a list AND has length 0)
      ! (is.list(x) && length(x) == 0)
    }, cleanedList2)
  } else {
    cleanedListFinal <- copy(cleanedList2)
  }
  return(cleanedListFinal)
}
