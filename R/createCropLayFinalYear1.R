createCropLayFinalYear1 <- function(Lay, potLayTopValid, runClusteringInParallel, clusterDistance){
  # 1. Crop and mask the seismicLines layer to the actual most probable area (potLayTopValid)
  cropLay <- Cache(postProcessTo, Lay, potLayTopValid)
  # 2. Convert all lines to polygons to exclude them from selecting the random point
  # DOUBLE CHECK I CAN'T MAKE POLYGONS OUT OF THE LINES USING extent=TRUE!
  cropLayBuf <- Cache(buffer, cropLay, width = 50)
  cropLayAg <- Cache(aggregate, cropLayBuf, dissolve = TRUE)
  finalPotLay <- Cache(erase, potLayTopValid, cropLayAg)
  # 4. Calculate the total length of lines in the "chosen area", their mean and sd 
  # NOTE: Although we have for the whole area, it is better to improve it with data from specific
  # the specific area we are working on.
  message("Getting the line's length, average and sd...")
  # C1. We have lots of duplicated objectids, which makes it harder later to identify lines 
  # individually (see C2). So we need to give a new individual value for each before the 
  # intersection
  cropLay$individualID <- 1:NROW(cropLay)
  tic("Time elapsed for intersecting potential and seismic lines layer: ")
  cropLay <- reproducible::Cache(terra::intersect(potLayTopValid, cropLay)) # Takes about 20min? 
  toc()
  # # C2. Some lines cross potential regions, generating 2 rows (i.e., 2 different Potential values)  
  # # per OBJECTID. I should deal with this by assigning these specific lines to one of them, randomly.
  ###### TO REMOVE OVERLAPPING
  tic("Time elapsed to identify overlapping features: ")
  overlap_matrix <- terra::relate(x = cropLay, relation = "T********", pairs = TRUE)
  overlapMatrix <- as.data.table(overlap_matrix)
  overlapMatrix[, self := fifelse(id.1 == id.2, TRUE, FALSE)] # Exclude "self-overlapping"
  overlapMatrix <- overlapMatrix[self == FALSE,]
  overlapMatrix[, self := NULL]
  
  # Ensure id.1 is always the smaller and id.2 is the larger
  overlapMatrix <- overlapMatrix[, .(id.1 = pmin(id.1, id.2), id.2 = pmax(id.1, id.2))] 
  overlapMatrix <- unique(overlapMatrix) # Remove duplicates
  
  to_remove <- numeric()
  for (i in 1:nrow(overlapMatrix)) {
    if (i %in% c(1, 1000)) message(paste0("Feature ", i, " of ", nrow(overlapMatrix)))
    if (i %% 10000 == 0) message(paste0("Feature ", i, " of ", nrow(overlapMatrix)))
    angle1 <- calculateLineAngle(cropLay[overlapMatrix[i,"id.1"]])
    angle2 <- calculateLineAngle(cropLay[overlapMatrix[i,"id.2"]])
    if (round(angle1, 4) != round(angle2, 4)){
      next
      # If different angles, they intersect but don't overlap
    } else {
      # If they overlap, by how much? We should exclude the shortest
      line1 <- perim(cropLay[overlapMatrix[i,"id.1"]])
      line2 <- perim(cropLay[overlapMatrix[i,"id.2"]])
      smallest <- if (line1 > line2) overlapMatrix[i, id.2] else overlapMatrix[i, id.1]
      to_remove <- c(to_remove, smallest)
    }
  }
  toc()
  
  cropLay <- cropLay[-to_remove, ]
  ###### TO REMOVE OVERLAPPING
  
  # Loop through the potentials, and make clusters
  # While looping, also get the average line lengths and sds 
  cropLayFinal <- do.call(rbind, lapply(unique(cropLay$Potential), function(Pot){
    cropLayInt <- Cache(clusterLines, cropLay[cropLay$Potential == Pot,], 
                        distThreshold = clusterDistance, 
                        currPotential = Pot,
                        runInParallel = runClusteringInParallel,
                        totPotential = length(unique(cropLay$Potential)))
    return(cropLayInt)
  }))
  
  return(cropLayFinal)
}