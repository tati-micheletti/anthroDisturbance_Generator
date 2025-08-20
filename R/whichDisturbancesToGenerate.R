whichDisturbancesToGenerate <- function(startTime,
                                        endTime,
                                        currentTime,
                                        disturbanceParameters) {
  intervals <- suppressWarnings(as.numeric(disturbanceParameters[["disturbanceInterval"]]))
  bad <- !is.finite(intervals) | intervals <= 0
  if (any(bad)) {
    warning(sprintf(
      "Skipping %d row(s) with non-positive or non-numeric disturbanceInterval: %s",
      sum(bad), paste(which(bad), collapse = ", ")
    ), call. = FALSE)
  }
  ok_idx <- which(!bad)
  hits <- vapply(ok_idx, function(i) {
    # avoids huge sequences; works cleanly for integer year steps
    if (currentTime < startTime || currentTime > endTime) return(NA_integer_)
    if (((currentTime - startTime) %% intervals[i]) == 0) i else NA_integer_
  }, integer(1))
  as.integer(na.omit(hits))
}
