writeDiagnostics <- function(sim, scheduledIdx = NULL) {
  # Summarize gating for all rows; for Generating rows include origin/potential counts
  dpar <- sim$disturbanceParameters
  if (!NROW(dpar)) return(invisible(NULL))

  # Helper to count features/raster cells
  nonEmptyCount <- function(x) {
    if (is.null(x)) return(0L)
    if (inherits(x, c("SpatVector", "sf", "Spatial"))) {
      return(as.integer(tryCatch(terra::nrow(x), error = function(...) 0L)))
    }
    if (inherits(x, c("RasterLayer", "SpatRaster"))) {
      rr <- if (inherits(x, "RasterLayer")) terra::rast(x) else x
      return(as.integer(tryCatch(sum(!is.na(terra::values(rr))), error = function(...) 0L)))
    }
    0L
  }

  # Scheduled set at this time — allow caller override to avoid double-computing
  if (is.null(scheduledIdx)) {
    scheduledIdx <- try(
      whichDisturbancesToGenerate(
        startTime = SpaDES.core::start(sim),
        currentTime = SpaDES.core::time(sim),
        endTime = SpaDES.core::end(sim),
        disturbanceParameters = dpar
      ),
      silent = TRUE
    )
    if (inherits(scheduledIdx, "try-error")) scheduledIdx <- integer(0)
  }

  rows <- lapply(seq_len(NROW(dpar)), function(i) {
    dn <- dpar$dataName[i]
    dc <- dpar$dataClass[i]
    do <- dpar$disturbanceOrigin[i]
    originLay <- try(sim$disturbanceList[[dn]][[do]], silent = TRUE)
    potLay    <- try(sim$disturbanceList[[dn]][[dc]], silent = TRUE)
    originCnt <- if (inherits(originLay, "try-error")) 0L else nonEmptyCount(originLay)
    potCnt    <- if (inherits(potLay, "try-error")) 0L else nonEmptyCount(potLay)
    rateVal   <- suppressWarnings(as.numeric(dpar$disturbanceRate[i]))
    hasSz     <- !is.na(dpar$disturbanceSize[i])
    sched     <- i %in% scheduledIdx
    # Reason classification
    reason <- if (dpar$disturbanceType[i] != "Generating") "not_generating" else
              if (potCnt == 0L) "no_potential" else
              if (is.na(rateVal) || rateVal <= 0) "no_rate" else
              if (!hasSz) "no_size" else
              if (!sched) "not_scheduled" else "ready"
    list(
      row = i,
      dataName = dn,
      dataClass = dc,
      disturbanceType = dpar$disturbanceType[i],
      origin = do,
      originCount = originCnt,
      potentialCount = potCnt,
      hasSize = hasSz,
      rate = rateVal,
      scheduled = sched,
      reason = reason
    )
  })
  df <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  # Print to console (first 20 rows)
  msg <- capture.output(print(as.data.frame(head(df, 20))))
  message(paste(c("Generation diagnostics (first 20):", msg), collapse = "\n"))
  message(sprintf("writeDiagnostics: snapshot rows=%d scheduled=%s", nrow(df), paste0(scheduledIdx, collapse = ",")))

  message("writeDiagnostics: sim class = ", paste(class(sim), collapse = "|"))
  # Write to CSV for later analysis
  outDir <- tryCatch({
    pths <- SpaDES.core::paths(sim)
    message("writeDiagnostics: available path names = ", paste(names(pths), collapse = ","))
    pths[["outputPath"]]
  }, error = function(e) {
    message("writeDiagnostics: failed retrieving paths -> ", conditionMessage(e))
    NULL
  })
  message("writeDiagnostics: resolved outDir = ", if (is.null(outDir)) "NULL" else as.character(outDir),
          " (length=", length(outDir), ", class=", paste(class(outDir), collapse = "|"), ")")
  runName <- tryCatch(SpaDES.core::P(sim)$.runName, error = function(e) "run")
  if (is.null(runName) || !nzchar(runName)) runName <- "run"
  if (length(outDir) == 1 && nzchar(outDir)) {
    writeErr <- tryCatch({
      timeFun <- get("time", envir = asNamespace("SpaDES.core"))
      currentTime <- timeFun(sim)
      fn <- file.path(outDir, sprintf("diagnosticsFull_%s_%s.csv", currentTime, runName))
      message("writeDiagnostics: attempting to write to ", fn)
      data.table::fwrite(df, fn)
      message("writeDiagnostics: wrote diagnostics CSV")
      FALSE
    }, error = function(e) {
      message("writeDiagnostics: failed to write diagnostics CSV -> ", conditionMessage(e))
      TRUE
    })
    if (isTRUE(writeErr)) message("writeDiagnostics: skipping CSV for this snapshot due to previous error.")
  }
  invisible(df)
}
