extractNonPotentialLayers <- function(disturbanceList) {
  # Initialize vectors to store the results
  all_sectors <- character()
  all_data_classes <- character()
  
  # Iterate through the top-level list names (Sectors)
  for (sector_name in names(disturbanceList)) {
    category_list <- disturbanceList[[sector_name]]
    
    # Iterate through the second-level list names (dataClass names)
    for (data_class_name in names(category_list)) {
      
      # Condition: Check if the layer name does NOT start with "potential"
      if (!startsWith(data_class_name, "potential")) {
        # If the condition is met, add the names to our vectors
        all_sectors <- c(all_sectors, sector_name)
        all_data_classes <- c(all_data_classes, data_class_name)
      }
    } # End inner loop
  } # End outer loop
  
  # Create the final data.table
  result_dt <- data.table(
    Sector = all_sectors,
    dataClass = all_data_classes
  )
  
  # Remove duplicate rows (handles cases like the duplicate 'mining/mining')
  result_dt <- unique(result_dt)
  
  return(result_dt)
}
