createCropLayFinalYear1_new <- function(Lay, 
                                        potLayTopValid, 
                                        runClusteringInParallel, 
                                        clusterDistance, 
                                        studyAreaHash){
  
  # 0. Input validation
  if (!is.numeric(clusterDistance) || length(clusterDistance) != 1) {
    stop("`clusterDistance` must be a single numeric value")
  }
  if (clusterDistance < 0) {
    stop("`clusterDistance` must be non‐negative")
  }
  
  # 1. Crop and mask the seismicLines layer to the actual most probable area (potLayTopValid)
  cropLay <- Cache(postProcessTo, Lay, potLayTopValid)
  
  # if no lines fall inside the potential mask, return an empty SpatVector with a warning
  if (nrow(cropLay) == 0L) {
    warning("createCropLayFinalYear1: no input lines fall within potLayTopValid — returning empty SpatVector")
    return(Lay[0, ])
  }
  
  # 2. Convert all lines to polygons to exclude them from selecting the random point
  cropLayBuf <- Cache(buffer, cropLay, width = 50)
  cropLayAg <- Cache(aggregate, cropLayBuf, dissolve = TRUE)
  finalPotLay <- Cache(erase, potLayTopValid, cropLayAg)
  potLayTopValid <- finalPotLay
  
  # 3. Assign unique ID
  cropLay$individualID <- seq_len(nrow(cropLay))
  
  # 4. Calculate the total length of lines in the "chosen area", their mean and sd
  # Remove overlapping
  tic("Time elapsed to identify overlapping features: ")
  overlap_matrix <- terra::relate(x = cropLay, relation = "T********", pairs = TRUE)
  overlapMatrix <- as.data.table(overlap_matrix)
  overlapMatrix[, self := fifelse(id.1 == id.2, TRUE, FALSE)] # Exclude "self-overlapping"
  overlapMatrix <- overlapMatrix[self == FALSE,]
  overlapMatrix[, self := NULL]
  
  # Ensure id.1 is always the smaller and id.2 is the larger
  overlapMatrix <- overlapMatrix[, .(id.1 = pmin(id.1, id.2), id.2 = pmax(id.1, id.2))] 
  overlapMatrix <- unique(overlapMatrix) # Remove duplicates
  
  to_remove <- integer()
  for (i in seq_len(nrow(overlapMatrix))) {
    if (i %in% c(1, 1000)) message(paste0("Feature ", i, " of ", nrow(overlapMatrix)))
    if (i %% 10000 == 0) message(paste0("Feature ", i, " of ", nrow(overlapMatrix)))
    angle1 <- calculateLineAngle(cropLay[overlapMatrix[i, "id.1"]])
    angle2 <- calculateLineAngle(cropLay[overlapMatrix[i, "id.2"]])
    # Skip if angle is NA
    if (is.na(angle1) || is.na(angle2)) next
    if (round(angle1, 4) != round(angle2, 4)){
      next
    } else {
      line1 <- perim(cropLay[overlapMatrix[i,"id.1"]])
      line2 <- perim(cropLay[overlapMatrix[i,"id.2"]])
      smallest <- if (line1 > line2) overlapMatrix[i, id.2] else overlapMatrix[i, id.1]
      to_remove <- c(to_remove, smallest)
    }
  }
  toc()
  if (length(to_remove) > 0) {
    cropLay <- cropLay[-unique(to_remove), ]
  }
  # End overlapping removal
  
  # Loop through the potentials, and make clusters
  cropLayFinal <- do.call(rbind, lapply(unique(cropLay$Potential), function(Pot){
    cropLayInt <- Cache(clusterLines, cropLay[cropLay$Potential == Pot,], 
                        distThreshold = clusterDistance, 
                        currPotential = Pot,
                        runInParallel = runClusteringInParallel,
                        totPotential = length(unique(cropLay$Potential)), 
                        userTags = studyAreaHash)
    return(cropLayInt)
  }))
  
  # now returns list of cropLayFinal AND updated line masked potLayTopValid
  return(list(
    lines           = cropLayFinal,
    availableArea   = potLayTopValid
  ))
}
