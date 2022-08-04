replaceList <- function(disturbanceList, 
                        updatedLayers){
  
  print("Develop replaceList function from anthroDisturbance_Generator. Make it generic and replace 
        the same function in potentialResourcesNT_DataPrep")
  browser()  
  # Oil and Gas
  disturbanceList$oilGas[names(disturbanceList$oilGas) == "potentialOilGas"] <- NULL
  disturbanceList$oilGas[["potentialOilGas"]] <- potentialOil
  
  # Mining
  disturbanceList$mining[names(disturbanceList$mining) == "potentialMining"] <- NULL
  disturbanceList$mining[["potentialMining"]] <- potentialMining
  
  # Cutblocks/Forestry
  disturbanceList$forestry[names(disturbanceList$forestry) == "potentialCutblocks"] <- NULL
  disturbanceList$forestry[["potentialCutblocks"]] <- potentialCutblocks
  
  return(disturbanceList)
}