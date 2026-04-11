#' Extract Non-Potential Layers
#'
#' @param disturbanceList A named list of sectors, each containing a named list of data classes.
#' @return A data.table with columns Sector and dataClass for layers not starting with 'potential'.

extractNonPotentialLayers <- function(disturbanceList) {
  # Validate input
  if (!is.list(disturbanceList) || length(disturbanceList) == 0) {
    return(data.table(Sector=character(), dataClass=character()))
  }
  
  # Build list of tables per sector
  dt_list <- lapply(seq_along(disturbanceList), function(i) {
    sector <- names(disturbanceList)[i]
    # Skip unnamed or non-list entries
    sublist <- disturbanceList[[i]]
    if (is.null(sector) || is.na(sector) || !is.list(sublist)) return(NULL)
    
    classes <- names(sublist)
    # Filter out NA names and those starting with 'potential'
    valid <- !is.na(classes) & !startsWith(classes, "potential")
    classes <- classes[valid]
    if (length(classes) == 0) return(NULL)
    
    data.table(Sector = rep(sector, length(classes)),
               dataClass = classes)
  })
  
  # Combine and dedupe
  result <- unique(rbindlist(dt_list, use.names = TRUE, fill = TRUE))
  return(result)
}
