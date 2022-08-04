whichDisturbancesToGenerate <- function(startTime,
                                        endTime,
                                        currentTime,
                                        disturbanceParameters){
    whichOnes <- na.omit(unlist(lapply(1:NROW(disturbanceParameters), function(INDEX){
      fullSeq <- seq(startTime, endTime, by = disturbanceParameters[INDEX, disturbanceInterval])
      if (currentTime %in% fullSeq){
        return(INDEX)
      } else {
        return(NA)
      }
    })))
  return(whichOnes)
}