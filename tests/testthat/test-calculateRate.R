library(testthat)
library(terra)
library(data.table)
library(digest)
library(withr)

# Helper: create a minimal SpatVector study area (one 1x1 cell)
create_study_area <- function() {
  r <- rast(nrows=1, ncols=1, xmin=0, xmax=1, ymin=0, ymax=1, crs="EPSG:4326")
  terra::as.polygons(r)
}

# Dummy disturbanceDT mapping
disturbanceDT <- data.table(
  fieldToSearch = "Class",
  classToSearch = "C1",
  dataName      = "dn1",
  dataClass     = "orig1"
)

# Dummy non-empty layer stub
dummy_layer <- list("not_null")

# Common inputs for tests
params_template <- function() {
  data.table(
    dataName          = "dn1",
    dataClass         = "orig1",
    disturbanceType   = "Generating",
    disturbanceOrigin = "orig1",
    disturbanceRate   = NA_real_
  )
}

rast_template <- function() {
  rast(nrows=1, ncols=1, xmin=0, xmax=1, ymin=0, ymax=1, crs="EPSG:4326")
}


# 1. Error if both DisturbanceRate and totalDisturbanceRate are provided
test_that("throws error when both DisturbanceRate and totalDisturbanceRate are non-null", {
  params <- params_template()
  expect_error(
    calculateRate(
      disturbanceParameters           = params,
      disturbanceDT                   = disturbanceDT,
      disturbanceList                 = list(dn1 = list(orig1 = dummy_layer)),
      whichToUpdate                   = 1L,
      RTM                             = rast_template(),
      studyArea                       = create_study_area(),
      destinationPath                 = tempdir(),
      DisturbanceRate                 = data.table(x=1),
      totalDisturbanceRate            = 10
    ),
    "Both DisturbanceRate and totalDisturbanceRate were provided"
  )
})

# 2. Using explicit DisturbanceRate table overrides calculation
test_that("applies user-supplied DisturbanceRate correctly and preserves other fields", {
  params <- params_template()
  override <- data.table(
    dataName          = "dn1",
    dataClass         = "orig1",
    disturbanceType   = "Generating",
    disturbanceOrigin = "orig1",
    disturbanceRate   = 7.5
  )
  
  result <- calculateRate(
    disturbanceParameters           = params,
    disturbanceDT                   = disturbanceDT,
    disturbanceList                 = list(dn1 = list(orig1 = dummy_layer)),
    whichToUpdate                   = 1L,
    RTM                             = rast_template(),
    studyArea                       = create_study_area(),
    destinationPath                 = tempdir(),
    DisturbanceRate                 = override,
    totalDisturbanceRate            = NULL
  )
  
  # Only one row returned, rate set correctly
  expect_equal(nrow(result), 1)
  expect_equal(result$disturbanceRate, 7.5)
  # Other columns unchanged
  expect_equal(result$dataName, params$dataName)
  expect_equal(result$dataClass, params$dataClass)
  expect_equal(result$disturbanceType, params$disturbanceType)
  expect_equal(result$disturbanceOrigin, params$disturbanceOrigin)
})

# 3. Layer absent leads to row being dropped when a dummy DisturbanceRate is supplied (bypassing ECCC fetch)
#    Also checks zero-length list is treated same as NULL
test_that("skips updating when the disturbance layer is NULL or empty list", {
  params <- params_template()
  # supply a non-null DisturbanceRate stub to skip empirical ECCC block
  stubRates <- data.table(
    dataName          = "dn1",
    dataClass         = "orig1",
    disturbanceType   = "Generating",
    disturbanceOrigin = "orig1",
    disturbanceRate   = NA_real_
  )
  for (lay in list(NULL, list())) {
    result <- calculateRate(
      disturbanceParameters             = params,
      disturbanceDT                     = disturbanceDT,
      disturbanceList                   = list(dn1 = list(orig1 = lay)),
      whichToUpdate                     = 1L,
      RTM                               = rast_template(),
      studyArea                         = create_study_area(),
      destinationPath                   = tempdir(),
      DisturbanceRate                   = stubRates,
      totalDisturbanceRate              = NULL,
      disturbanceRateRelatesToBufferedArea = FALSE,
      maskOutLinesFromPolys               = FALSE,
      aggregateSameDisturbances           = FALSE,
      archiveNEW                          = NULL,
      targetFileNEW                       = NULL,
      urlNEW                              = NULL,
      archiveOLD                          = NULL,
      targetFileOLD                       = NULL,
      urlOLD                              = NULL
    )
    expect_equal(nrow(result), 0)
  }
})

# 4. Empirical fallback: negative change is clamped to zero
test_that("negative empirical change yields zero disturbance rate without warning", {
  # stub disturbanceInfoFromECCC to return negative growth
  stub_data <- data.table(Class = "C1", yearOLD = 5, yearNEW = 0)
  original_fn <- disturbanceInfoFromECCC
  assign("disturbanceInfoFromECCC", function(...) list(AD_changed = stub_data), envir = .GlobalEnv)
  on.exit(assign("disturbanceInfoFromECCC", original_fn, envir = .GlobalEnv), add = TRUE)
  
  params <- params_template()
  # suppress messages rather than expect_silent, since calculateRate uses message()
  suppressMessages(
    result <- calculateRate(
      disturbanceParameters             = params,
      disturbanceDT                     = disturbanceDT,
      disturbanceList                   = list(dn1 = list(orig1 = dummy_layer)),
      whichToUpdate                     = 1L,
      RTM                               = rast_template(),
      studyArea                         = create_study_area(),
      destinationPath                   = tempdir(),
      DisturbanceRate                   = NULL,
      totalDisturbanceRate              = NULL,
      disturbanceRateRelatesToBufferedArea = FALSE,
      maskOutLinesFromPolys               = FALSE,
      aggregateSameDisturbances           = FALSE,
      archiveNEW                          = NULL,
      targetFileNEW                       = NULL,
      urlNEW                              = NULL,
      archiveOLD                          = NULL,
      targetFileOLD                       = NULL,
      urlOLD                              = NULL
    )
  )
  expect_equal(nrow(result), 1)
  expect_equal(result$disturbanceRate, 0)
})

# 5. Empirical fallback: positive change yields correct disturbance rate
test_that("positive empirical change yields correct disturbance rate within tolerance", {
  # stub disturbanceInfoFromECCC to return positive growth
  stub_data <- data.table(Class = "C1", yearOLD = 0, yearNEW = 5)
  original_fn <- disturbanceInfoFromECCC
  assign("disturbanceInfoFromECCC", function(...) list(AD_changed = stub_data), envir = .GlobalEnv)
  on.exit(assign("disturbanceInfoFromECCC", original_fn, envir = .GlobalEnv), add = TRUE)
  
  params <- params_template()
  # suppress messages
  suppressMessages(
    result <- calculateRate(
      disturbanceParameters             = params,
      disturbanceDT                     = disturbanceDT,
      disturbanceList                   = list(dn1 = list(orig1 = dummy_layer)),
      whichToUpdate                     = 1L,
      RTM                               = rast_template(),
      studyArea                         = create_study_area(),
      destinationPath                   = tempdir(),
      DisturbanceRate                   = NULL,
      totalDisturbanceRate              = NULL,
      disturbanceRateRelatesToBufferedArea = FALSE,
      maskOutLinesFromPolys               = FALSE,
      aggregateSameDisturbances           = FALSE,
      archiveNEW                          = NULL,
      targetFileNEW                       = NULL,
      urlNEW                              = NULL,
      archiveOLD                          = NULL,
      targetFileOLD                       = NULL,
      urlOLD                              = NULL
    )
  )
  # calculate expected: ((NEW - OLD)/5 years)/area * 100
  studyV <- terra::aggregate(create_study_area())
  totalArea <- terra::expanse(studyV, unit = "km", transform = FALSE)
  expected_rate <- ((stub_data$yearNEW - stub_data$yearOLD)/5)/totalArea*100
  expect_equal(result$disturbanceRate, expected_rate, tolerance = 1e-8)
})

# 6. Test for totalDisturbanceRate
test_that("correctly applies totalDisturbanceRate when DisturbanceRate is NULL", {
  params <- params_template()
  stub_data <- data.table(Class = "C1", yearOLD = 0, yearNEW = 5)
  proportionTable <- data.table(
    dataClass = "orig1",
    proportionOfTotalDisturbance = 0.5,
    calculatedDisturbanceProportion = NA_real_
  )
  
  original_fn <- disturbanceInfoFromECCC
  assign("disturbanceInfoFromECCC", function(...) list(AD_changed = stub_data, proportionTable = proportionTable), envir = .GlobalEnv)
  on.exit(assign("disturbanceInfoFromECCC", original_fn, envir = .GlobalEnv), add = TRUE)
  
  expect_warning(  
    result <- calculateRate(
      disturbanceParameters             = params,
      disturbanceDT                     = disturbanceDT,
      disturbanceList                   = list(dn1 = list(orig1 = dummy_layer)),
      whichToUpdate                     = 1L,
      RTM                               = rast_template(),
      studyArea                         = create_study_area(),
      destinationPath                   = tempdir(),
      DisturbanceRate                   = NULL,
      totalDisturbanceRate              = 10,  # Example total disturbance rate
      disturbanceRateRelatesToBufferedArea = FALSE,
      maskOutLinesFromPolys               = FALSE,
      aggregateSameDisturbances           = FALSE,
      archiveNEW                          = NULL,
      targetFileNEW                       = NULL,
      urlNEW                              = NULL,
      archiveOLD                          = NULL,
      targetFileOLD                       = NULL,
      urlOLD                              = NULL
    ),
    "The totalDisturbanceRate was supplied as 10"
  )
  
  # Check if the calculated rate is proportional to totalDisturbanceRate
  expected_rate <- 10 * proportionTable$proportionOfTotalDisturbance
  expect_equal(result$disturbanceRate, expected_rate, tolerance = 1e-8)
})

# 7. Test for warning when diffYears is OLD_2015
test_that("issues warning when diffYears is OLD_2015", {
  original <- get("disturbanceInfoFromECCC", envir = .GlobalEnv)
  assign("disturbanceInfoFromECCC", function(...) list(AD_changed = data.table(Class=character(), yearOLD=integer(), yearNEW=integer())), envir = .GlobalEnv)
  on.exit(assign("disturbanceInfoFromECCC", original, envir = .GlobalEnv), add = TRUE)
  
  params <- params_template()
  expect_warning(
    calculateRate(
      disturbanceParameters             = params,
      disturbanceDT                     = disturbanceDT,
      disturbanceList                   = list(dn1 = list(orig1 = dummy_layer)),
      whichToUpdate                     = 1L,
      RTM                               = rast_template(),
      studyArea                         = create_study_area(),
      destinationPath                   = tempdir(),
      diffYears                         = "OLD_2015",  # Trigger warning
      DisturbanceRate                   = NULL,
      totalDisturbanceRate              = NULL,
      disturbanceRateRelatesToBufferedArea = FALSE,
      maskOutLinesFromPolys               = FALSE,
      aggregateSameDisturbances           = FALSE,
      archiveNEW                          = NULL,
      targetFileNEW                       = NULL,
      urlNEW                              = NULL,
      archiveOLD                          = NULL,
      targetFileOLD                       = NULL,
      urlOLD                              = NULL
    ),
    "While ECCC OLD footprint layer covers the whole country"
  )
})

# 8. Test for zero empirical change
test_that("zero empirical change yields zero disturbance rate", {
  stub_data <- data.table(Class = "C1", yearOLD = 5, yearNEW = 5)
  original_fn <- disturbanceInfoFromECCC
  assign("disturbanceInfoFromECCC", function(...) list(AD_changed = stub_data), envir = .GlobalEnv)
  on.exit(assign("disturbanceInfoFromECCC", original_fn, envir = .GlobalEnv), add = TRUE)
  
  params <- params_template()
  suppressMessages(
    result <- calculateRate(
      disturbanceParameters             = params,
      disturbanceDT                     = disturbanceDT,
      disturbanceList                   = list(dn1 = list(orig1 = dummy_layer)),
      whichToUpdate                     = 1L,
      RTM                               = rast_template(),
      studyArea                         = create_study_area(),
      destinationPath                   = tempdir(),
      DisturbanceRate                   = NULL,
      totalDisturbanceRate              = NULL,
      disturbanceRateRelatesToBufferedArea = FALSE,
      maskOutLinesFromPolys               = FALSE,
      aggregateSameDisturbances           = FALSE,
      archiveNEW                          = NULL,
      targetFileNEW                       = NULL,
      urlNEW                              = NULL,
      archiveOLD                          = NULL,
      targetFileOLD                       = NULL,
      urlOLD                              = NULL
    )
  )
  expect_equal(result$disturbanceRate, 0)
})

# 9. Test for multiple disturbances
test_that("correctly handles multiple disturbances", {
  params <- rbindlist(list(
    params_template(),
    data.table(
      dataName          = "dn2",
      dataClass         = "orig2",
      disturbanceType   = "Generating",
      disturbanceOrigin = "orig2",
      disturbanceRate   = NA_real_
    )
  ))
  
  disturbanceDT_multi <- rbindlist(list(
    disturbanceDT,
    data.table(
      fieldToSearch = "Class",
      classToSearch = "C2",
      dataName      = "dn2",
      dataClass     = "orig2"
    )
  ))
  
  dummy_layer2 <- list("not_null")
  disturbanceList_multi <- list(
    dn1 = list(orig1 = dummy_layer),
    dn2 = list(orig2 = dummy_layer2)
  )
  
  stub_data <- data.table(
    Class = c("C1", "C2"),
    yearOLD = c(0, 0),
    yearNEW = c(5, 10)
  )
  
  original_fn <- disturbanceInfoFromECCC
  assign("disturbanceInfoFromECCC", function(...) list(AD_changed = stub_data), envir = .GlobalEnv)
  on.exit(assign("disturbanceInfoFromECCC", original_fn, envir = .GlobalEnv), add = TRUE)
  
  result <- calculateRate(
    disturbanceParameters             = params,
    disturbanceDT                     = disturbanceDT_multi,
    disturbanceList                   = disturbanceList_multi,
    whichToUpdate                     = 1:2,
    RTM                               = rast_template(),
    studyArea                         = create_study_area(),
    destinationPath                   = tempdir(),
    DisturbanceRate                   = NULL,
    totalDisturbanceRate              = NULL,
    disturbanceRateRelatesToBufferedArea = FALSE,
    maskOutLinesFromPolys               = FALSE,
    aggregateSameDisturbances           = FALSE,
    archiveNEW                          = NULL,
    targetFileNEW                       = NULL,
    urlNEW                              = NULL,
    archiveOLD                          = NULL,
    targetFileOLD                       = NULL,
    urlOLD                              = NULL
  )
  
  studyV <- terra::aggregate(create_study_area())
  totalArea <- terra::expanse(studyV, unit = "km", transform = FALSE)
  expected_rate1 <- ((5 - 0)/5)/totalArea*100
  expected_rate2 <- ((10 - 0)/5)/totalArea*100
  
  expect_equal(result$disturbanceRate[1], expected_rate1, tolerance = 1e-8)
  expect_equal(result$disturbanceRate[2], expected_rate2, tolerance = 1e-8)
})

# 10. Test for cache file handling
test_that("uses cached file if available", {
  originalECCC <- get("disturbanceInfoFromECCC", envir = .GlobalEnv)
  assign("disturbanceInfoFromECCC", function(...){
    # return only AD_changed so calculateRate_new reads your CSV but never calls prepInputs()
    list(AD_changed = data.table(Class="C1", yearOLD=0, yearNEW=5))
  }, envir = .GlobalEnv)
  on.exit(assign("disturbanceInfoFromECCC", originalECCC, envir = .GlobalEnv), add = TRUE)
  
  params <- params_template()
  key   <- digest::digest(create_study_area())
  fname <- file.path(tempdir(), paste0("anthropogenicDisturbance_ECCC_2010_2015_", key, ".csv"))
  write.csv(data.table(Class="C1", yearOLD=0, yearNEW=5), fname, row.names=FALSE)
  
  result <- calculateRate(
    disturbanceParameters             = params,
    disturbanceDT                     = disturbanceDT,
    disturbanceList                   = list(dn1 = list(orig1 = dummy_layer)),
    whichToUpdate                     = 1L,
    RTM                               = rast_template(),
    studyArea                         = create_study_area(),
    destinationPath                   = tempdir(),
    DisturbanceRate                   = NULL,
    totalDisturbanceRate              = NULL,
    disturbanceRateRelatesToBufferedArea = FALSE,
    maskOutLinesFromPolys               = FALSE,
    aggregateSameDisturbances           = FALSE,
    archiveNEW                          = NULL,
    targetFileNEW                       = NULL,
    urlNEW                              = NULL,
    archiveOLD                          = NULL,
    targetFileOLD                       = NULL,
    urlOLD                              = NULL
  )
  
  studyV <- terra::aggregate(create_study_area())
  totalArea <- terra::expanse(studyV, unit = "km", transform = FALSE)
  expected_rate <- ((5 - 0)/5)/totalArea*100
  expect_equal(result$disturbanceRate, expected_rate, tolerance = 1e-8)
  
  # Clean up
  unlink(fname)
})

# 11.
test_that("totalDisturbanceRate distributes rates correctly", {
  # Mock disturbanceInfoFromECCC to return proportion table
  stub_data <- list(
    AD_changed = data.table(Class = "C1", yearOLD = 1, yearNEW = 2),
    proportionTable = data.table(dataClass = "orig1", proportionOfTotalDisturbance = 0.5)
  )
  original_fn <- disturbanceInfoFromECCC
  assign("disturbanceInfoFromECCC", function(...) stub_data, envir = .GlobalEnv)
  on.exit(assign("disturbanceInfoFromECCC", original_fn, envir = .GlobalEnv), add = TRUE)
  
  params <- params_template()
  
  expect_warning(
    result <- calculateRate(
      disturbanceParameters           = params,
      disturbanceDT                   = disturbanceDT,
      disturbanceList                 = list(dn1 = list(orig1 = dummy_layer)),
      whichToUpdate                   = 1L,
      RTM                             = rast_template(),
      studyArea                       = create_study_area(),
      destinationPath                 = tempdir(),
      DisturbanceRate                 = NULL,
      totalDisturbanceRate            = 10, # 10% total rate
      diffYears                       = "2010_2015",
      disturbanceRateRelatesToBufferedArea = FALSE,
      maskOutLinesFromPolys           = FALSE,
      aggregateSameDisturbances       = FALSE,
      archiveNEW                      = NULL,
      targetFileNEW                   = NULL,
      urlNEW                          = NULL,
      archiveOLD                      = NULL,
      targetFileOLD                   = NULL,
      urlOLD                          = NULL
    ), "The totalDisturbanceRate was supplied as 10"
  )
  
  # Expect 50% of total rate (0.5 * 10 = 5)
  expect_equal(result$disturbanceRate, 5)
})

# 12.
test_that("disturbance is dropped when missing from proportion table", {
  # Proportion table doesn't include our disturbance
  stub_data <- list(
    proportionTable = data.table(dataClass = "other_disturbance", proportionOfTotalDisturbance = 1)
  )
  original_fn <- disturbanceInfoFromECCC
  assign("disturbanceInfoFromECCC", function(...) list(
    AD_changed = data.table(Class=character(), yearOLD=integer(), yearNEW=integer()),
    proportionTable = data.table(dataClass="other", proportionOfTotalDisturbance=1)
  ), envir = .GlobalEnv)
  on.exit(assign("disturbanceInfoFromECCC", original_fn, envir = .GlobalEnv), add = TRUE)
  
  params <- params_template()
  
  result <- suppressWarnings(calculateRate(
    disturbanceParameters           = params,
    disturbanceDT                   = disturbanceDT,
    disturbanceList                 = list(dn1 = list(orig1 = dummy_layer)),
    whichToUpdate                   = 1L,
    RTM                             = rast_template(),
    studyArea                       = create_study_area(),
    destinationPath                 = tempdir(),
    DisturbanceRate                 = NULL,
    totalDisturbanceRate            = 10,
    diffYears                       = "2010_2015",
    disturbanceRateRelatesToBufferedArea = FALSE,
    maskOutLinesFromPolys           = FALSE,
    aggregateSameDisturbances       = FALSE,
    archiveNEW                      = NULL,
    targetFileNEW                   = NULL,
    urlNEW                          = NULL,
    archiveOLD                      = NULL,
    targetFileOLD                   = NULL,
    urlOLD                          = NULL
  ))
  
  key <- params[1, .(dataName, dataClass, disturbanceType, disturbanceOrigin)]
  updated <- result[key, on=names(key)]
  expect_equal(updated$disturbanceRate, 0)
  expect_warning(calculateRate(
    disturbanceParameters           = params,
    disturbanceDT                   = disturbanceDT,
    disturbanceList                 = list(dn1 = list(orig1 = dummy_layer)),
    whichToUpdate                   = 1L,
    RTM                             = rast_template(),
    studyArea                       = create_study_area(),
    destinationPath                 = tempdir(),
    DisturbanceRate                 = NULL,
    totalDisturbanceRate            = 10,
    diffYears                       = "2010_2015",
    disturbanceRateRelatesToBufferedArea = FALSE,
    maskOutLinesFromPolys           = FALSE,
    aggregateSameDisturbances       = FALSE,
    archiveNEW                      = NULL,
    targetFileNEW                   = NULL,
    urlNEW                          = NULL,
    archiveOLD                      = NULL,
    targetFileOLD                   = NULL,
    urlOLD                          = NULL
  ), "No historical ECCC disturbance for these classes → setting rate to zero", ignore.case = TRUE)
})

# 13.
test_that("handles multiple whichToUpdate rows", {
  # Create two disturbances to update
  params <- rbind(
    params_template(),
    data.table(
      dataName          = "dn2",
      dataClass         = "orig2",
      disturbanceType   = "Generating",
      disturbanceOrigin = "orig2",
      disturbanceRate   = NA_real_
    )
  )
  whichToUpdate <- c(TRUE, TRUE)
  
  # Extend disturbanceDT
  disturbanceDT2 <- rbind(
    disturbanceDT,
    data.table(
      fieldToSearch = "Class",
      classToSearch = "C2",
      dataName      = "dn2",
      dataClass     = "orig2"
    )
  )
  
  # Mock ECCC data for both classes
  stub_data <- list(
    AD_changed = data.table(
      Class = c("C1", "C2"),
      yearOLD = c(1, 2),
      yearNEW = c(3, 4)
    ))
  original_fn <- disturbanceInfoFromECCC
  assign("disturbanceInfoFromECCC", function(...) stub_data, envir = .GlobalEnv)
  on.exit(assign("disturbanceInfoFromECCC", original_fn, envir = .GlobalEnv), add = TRUE)
  
  result <- suppressMessages(
    calculateRate(
      disturbanceParameters           = params,
      disturbanceDT                   = disturbanceDT2,
      disturbanceList                 = list(
        dn1 = list(orig1 = dummy_layer),
        dn2 = list(orig2 = dummy_layer)
      ),
      whichToUpdate                   = whichToUpdate,
      RTM                             = rast_template(),
      studyArea                       = create_study_area(),
      destinationPath                 = tempdir(),
      DisturbanceRate                 = NULL,
      totalDisturbanceRate            = NULL,
      diffYears                       = "2010_2015",
      disturbanceRateRelatesToBufferedArea = FALSE,
      maskOutLinesFromPolys           = FALSE,
      aggregateSameDisturbances       = FALSE,
      archiveNEW                      = NULL,
      targetFileNEW                   = NULL,
      urlNEW                          = NULL,
      archiveOLD                      = NULL,
      targetFileOLD                   = NULL,
      urlOLD                          = NULL
    )
  )
  
  expect_equal(nrow(result), 2)
  expect_true(all(result$disturbanceRate > 0))
})

# 14.  
test_that("final rate message includes correct total", {
  params <- params_template()
  stub_data <- list(AD_changed = data.table(Class = "C1", yearOLD = 0, yearNEW = 5))
  original_fn <- disturbanceInfoFromECCC
  assign("disturbanceInfoFromECCC", function(...) stub_data, envir = .GlobalEnv)
  on.exit(assign("disturbanceInfoFromECCC", original_fn, envir = .GlobalEnv), add = TRUE)
  
  studyArea <- create_study_area()
  studyV <- terra::aggregate(studyArea)
  totalArea <- terra::expanse(studyV, unit = "km", transform = FALSE)
  expected_rate <- ((5-0)/5/totalArea)*100
  
  expect_message(
    calculateRate(
      disturbanceParameters           = params,
      disturbanceDT                   = disturbanceDT,
      disturbanceList                 = list(dn1 = list(orig1 = dummy_layer)),
      whichToUpdate                   = 1L,
      RTM                             = rast_template(),
      studyArea                       = studyArea,
      destinationPath                 = tempdir(),
      DisturbanceRate                 = NULL,
      totalDisturbanceRate            = NULL,
      diffYears                       = "2010_2015",
      disturbanceRateRelatesToBufferedArea = FALSE,
      maskOutLinesFromPolys           = FALSE,
      aggregateSameDisturbances       = FALSE,
      archiveNEW                      = NULL,
      targetFileNEW                   = NULL,
      urlNEW                          = NULL,
      archiveOLD                      = NULL,
      targetFileOLD                   = NULL,
      urlOLD                          = NULL
    ),
    paste0("Total expected yearly rate.*", round(expected_rate, 4), "%")
  )
})

# 15. DisturbanceRate provided but layer absent (should drop disturbance)
test_that("drops disturbance when layer absent but DisturbanceRate provided", {
  params <- params_template()
  override <- data.table(
    dataName          = "dn1",
    dataClass         = "orig1",
    disturbanceType   = "Generating",
    disturbanceOrigin = "orig1",
    disturbanceRate   = 7.5
  )
  
  result <- calculateRate(
    disturbanceParameters           = params,
    disturbanceList                 = list(dn1 = list(orig1 = NULL)),  # Layer absent
    DisturbanceRate                 = override,
    disturbanceDT                   = disturbanceDT,
    whichToUpdate                   = 1L,
    RTM                             = rast_template(),
    studyArea                       = create_study_area(),
    destinationPath                 = tempdir(),
    totalDisturbanceRate            = NULL, 
    diffYears                       = "2010_2015",
    disturbanceRateRelatesToBufferedArea = FALSE,
    maskOutLinesFromPolys           = FALSE,
    aggregateSameDisturbances       = FALSE,
    archiveNEW                      = NULL,
    targetFileNEW                   = NULL,
    urlNEW                          = NULL,
    archiveOLD                      = NULL,
    targetFileOLD                   = NULL,
    urlOLD                          = NULL   
  )
  expect_equal(nrow(result), 0)  # Should be dropped
})

# 16. DisturbanceRate table missing matching row (should error)
test_that("errors when DisturbanceRate has no matching row", {
  params <- params_template()
  override <- data.table(
    dataName="invalid",
    dataClass="invalid",
    disturbanceType="Generating",
    disturbanceOrigin="orig1",
    disturbanceRate=7.5
  )
  
  expect_warning(
    calculateRate(
      disturbanceParameters           = params,
      disturbanceList                 = list(dn1 = list(orig1 = dummy_layer)),
      DisturbanceRate                 = override,
      disturbanceDT                   = disturbanceDT,
      whichToUpdate                   = 1L,
      RTM                             = rast_template(),
      studyArea                       = create_study_area(),
      destinationPath                 = tempdir(),
      totalDisturbanceRate            = NULL, 
      diffYears                       = "2010_2015",
      disturbanceRateRelatesToBufferedArea = FALSE,
      maskOutLinesFromPolys           = FALSE,
      aggregateSameDisturbances       = FALSE,
      archiveNEW                      = NULL,
      targetFileNEW                   = NULL,
      urlNEW                          = NULL,
      archiveOLD                      = NULL,
      targetFileOLD                   = NULL,
      urlOLD                          = NULL
    ),
    "No matching row in DisturbanceRate for dn1/orig1; setting rate to zero"  # Check error message
  )
})

# 17. Missing DisturbanceRate row
test_that("errors cleanly when DisturbanceRate has no matching row", {
  params <- params_template()
  override <- data.table(
    dataName="invalid",
    dataClass="invalid",
    disturbanceType="Generating",
    disturbanceOrigin="orig1",
    disturbanceRate=7.5
  )
  
  expect_warning(
    result <- calculateRate(
      disturbanceParameters = params,
      disturbanceDT = disturbanceDT,
      disturbanceList = list(dn1 = list(orig1 = dummy_layer)),
      DisturbanceRate = override,
      whichToUpdate                   = 1L,
      RTM                             = rast_template(),
      studyArea                       = create_study_area(),
      destinationPath                 = tempdir(),
      totalDisturbanceRate            = NULL, 
      diffYears                       = "2010_2015",
      disturbanceRateRelatesToBufferedArea = FALSE,
      maskOutLinesFromPolys           = FALSE,
      aggregateSameDisturbances       = FALSE,
      archiveNEW                      = NULL,
      targetFileNEW                   = NULL,
      urlNEW                          = NULL,
      archiveOLD                      = NULL,
      targetFileOLD                   = NULL,
      urlOLD                          = NULL
    ),
    "No matching row in DisturbanceRate"
  )
  
  key <- params[1, .(dataName, dataClass, disturbanceType, disturbanceOrigin)]
  updated <- result[key, on = names(key)]
  expect_equal(nrow(updated), 1L)
  expect_equal(updated$disturbanceRate, 0)
  
  # Optional: whole table preserved
  expect_equal(nrow(result), nrow(params))
})

# 18. Missing class in ECCC data
test_that("skips disturbance when class missing in ECCC data", {
  stub_data <- data.table(Class = "OTHER_CLASS", yearOLD = 5, yearNEW = 10)
  original_fn <- disturbanceInfoFromECCC
  assign("disturbanceInfoFromECCC", function(...) list(AD_changed = stub_data), 
         envir = .GlobalEnv)
  on.exit(assign("disturbanceInfoFromECCC", original_fn, envir = .GlobalEnv), add = TRUE)
  
  params <- params_template()
  
  result <- suppressWarnings(calculateRate(
    disturbanceParameters             = params,
    disturbanceDT                     = disturbanceDT,
    disturbanceList                   = list(dn1 = list(orig1 = dummy_layer)),
    whichToUpdate                     = 1L,
    RTM                               = rast_template(),
    studyArea                         = create_study_area(),
    destinationPath                   = tempdir(),
    DisturbanceRate                   = NULL,
    totalDisturbanceRate              = NULL,
    disturbanceRateRelatesToBufferedArea = FALSE,
    maskOutLinesFromPolys               = FALSE,
    aggregateSameDisturbances           = FALSE,
    archiveNEW                          = NULL,
    targetFileNEW                       = NULL,
    urlNEW                              = NULL,
    archiveOLD                          = NULL,
    targetFileOLD                       = NULL,
    urlOLD                              = NULL
  ))
  
  key <- params[1, .(dataName, dataClass, disturbanceType, disturbanceOrigin)]
  updated <- result[key, on=names(key)]
  expect_equal(updated$disturbanceRate, 0)
  
  expect_warning(calculateRate(
    disturbanceParameters             = params,
    disturbanceDT                     = disturbanceDT,
    disturbanceList                   = list(dn1 = list(orig1 = dummy_layer)),
    whichToUpdate                     = 1L,
    RTM                               = rast_template(),
    studyArea                         = create_study_area(),
    destinationPath                   = tempdir(),
    DisturbanceRate                   = NULL,
    totalDisturbanceRate              = NULL,
    disturbanceRateRelatesToBufferedArea = FALSE,
    maskOutLinesFromPolys               = FALSE,
    aggregateSameDisturbances           = FALSE,
    archiveNEW                          = NULL,
    targetFileNEW                       = NULL,
    urlNEW                              = NULL,
    archiveOLD                          = NULL,
    targetFileOLD                       = NULL,
    urlOLD                              = NULL
  ), "No historical ECCC disturbance for these classes → setting rate to zero", ignore.case = TRUE)
})

# 19. Multiple DisturbanceRate matches
test_that("warns and uses first row when multiple DisturbanceRate matches", {
  params <- params_template()
  override <- rbind(
    data.table(
      dataName          = "dn1",
      dataClass         = "orig1",
      disturbanceType   = "Generating",
      disturbanceOrigin = "orig1",
      disturbanceRate   = 7.5
    ),
    data.table(
      dataName          = "dn1",
      dataClass         = "orig1",
      disturbanceType   = "Generating",
      disturbanceOrigin = "orig1",
      disturbanceRate   = 10.5
    )
  )
  
  expect_warning(
    result <- calculateRate(
      disturbanceParameters = params,
      disturbanceDT = disturbanceDT,
      disturbanceList = list(dn1 = list(orig1 = dummy_layer)),
      DisturbanceRate = override,
      whichToUpdate                   = 1L,
      RTM                             = rast_template(),
      studyArea                       = create_study_area(),
      destinationPath                 = tempdir(),
      totalDisturbanceRate            = NULL, 
      diffYears                       = "2010_2015",
      disturbanceRateRelatesToBufferedArea = FALSE,
      maskOutLinesFromPolys           = FALSE,
      aggregateSameDisturbances       = FALSE,
      archiveNEW                      = NULL,
      targetFileNEW                   = NULL,
      urlNEW                          = NULL,
      archiveOLD                      = NULL,
      targetFileOLD                   = NULL,
      urlOLD                          = NULL
    ),
    "Multiple matches"
  )
  expect_equal(result$disturbanceRate, 7.5)
})

# 20. Test that buffered/masking/aggregation flags are forwarded
test_that("buffer, maskOutLinesFromPolys and aggregateSameDisturbances are passed to disturbanceInfoFromECCC", {
  original_fn <- get("disturbanceInfoFromECCC", envir = .GlobalEnv)
  captured_args <- NULL
  assign("disturbanceInfoFromECCC", function(...){
    captured_args <<- list(...)
    list(AD_changed = data.table(Class = "C1", yearOLD = 0, yearNEW = 0))
  }, envir = .GlobalEnv)
  on.exit(assign("disturbanceInfoFromECCC", original_fn, envir = .GlobalEnv), add = TRUE)
  
  params <- params_template()
  suppressMessages(
    calculateRate(
      disturbanceParameters             = params,
      disturbanceDT                     = disturbanceDT,
      disturbanceList                   = list(dn1 = list(orig1 = dummy_layer)),
      whichToUpdate                     = 1L,
      RTM                               = rast_template(),
      studyArea                         = create_study_area(),
      destinationPath                   = tempdir(),
      DisturbanceRate                   = NULL,
      totalDisturbanceRate              = NULL,
      disturbanceRateRelatesToBufferedArea = TRUE,
      maskOutLinesFromPolys               = TRUE,
      aggregateSameDisturbances           = TRUE,
      archiveNEW                          = NULL,
      targetFileNEW                       = NULL,
      urlNEW                              = NULL,
      archiveOLD                          = NULL,
      targetFileOLD                       = NULL,
      urlOLD                              = NULL
    )
  )
  
  expect_true(captured_args$bufferedDisturbances,      info = "buffer flag not forwarded")
  expect_true(captured_args$maskOutLinesFromPolys,     info = "maskOutLinesFromPolys not forwarded")
  expect_true(captured_args$aggregateSameDisturbances, info = "aggregateSameDisturbances not forwarded")
})

# -------------------------------------------------------------------
# 21. Test empty whichToUpdate → returns original unmodified
test_that("no-op when whichToUpdate is empty", {
  params <- params_template()
  result <- calculateRate(
    disturbanceParameters           = params,
    disturbanceDT                   = disturbanceDT,
    disturbanceList                 = list(dn1 = list(orig1 = dummy_layer)),
    whichToUpdate                   = integer(0),
    RTM                             = rast_template(),
    studyArea                       = create_study_area(),
    destinationPath                 = tempdir(),
    DisturbanceRate                 = NULL,
    totalDisturbanceRate            = NULL,
    disturbanceRateRelatesToBufferedArea = FALSE,
    maskOutLinesFromPolys               = FALSE,
    aggregateSameDisturbances           = FALSE,
    archiveNEW                          = NULL,
    targetFileNEW                       = NULL,
    urlNEW                              = NULL,
    archiveOLD                          = NULL,
    targetFileOLD                       = NULL,
    urlOLD                              = NULL
  )
  expect_equal(result, params)
})

# -------------------------------------------------------------------
# 22. Test override for an “Enlarging” disturbanceType
test_that("applies user‐supplied DisturbanceRate for Enlarging classes", {
  params <- data.table(
    dataName          = "dn1",
    dataClass         = "orig1",
    disturbanceType   = "Enlarging",
    disturbanceOrigin = "orig1",
    disturbanceRate   = NA_real_
  )
  override <- data.table(
    dataName          = "dn1",
    dataClass         = "orig1",
    disturbanceType   = "Enlarging",
    disturbanceOrigin = "orig1",
    disturbanceRate   = 9.1
  )
  result <- calculateRate(
    disturbanceParameters           = params,
    disturbanceDT                   = disturbanceDT,
    disturbanceList                 = list(dn1 = list(orig1 = dummy_layer)),
    whichToUpdate                   = 1L,
    RTM                             = rast_template(),
    studyArea                       = create_study_area(),
    destinationPath                 = tempdir(),
    DisturbanceRate                 = override,
    totalDisturbanceRate            = NULL,
    disturbanceRateRelatesToBufferedArea = FALSE,
    maskOutLinesFromPolys               = FALSE,
    aggregateSameDisturbances           = FALSE,
    archiveNEW                          = NULL,
    targetFileNEW                       = NULL,
    urlNEW                              = NULL,
    archiveOLD                          = NULL,
    targetFileOLD                       = NULL,
    urlOLD                              = NULL
  )
  expect_equal(nrow(result), 1L)
  expect_equal(result$disturbanceRate, 9.1)
})

test_that("no rows to update returns original table unchanged", {
  skip_on_cran()
  
  # --- Minimal fixtures ---
  rtm <- rast(nrows=2, ncols=2, xmin=0, xmax=20, ymin=0, ymax=20, vals=1); crs(rtm) <- "EPSG:3857"
  sa  <- as.polygons(ext(rtm)); crs(sa) <- crs(rtm)
  
  disturbanceParameters <- data.table(
    dataName="forestry", dataClass="cutblocks", disturbanceType="Generating",
    disturbanceOrigin="cutblocks", disturbanceRate=NA_real_
  )
  
  disturbanceDT <- data.table(
    fieldToSearch="Class", classToSearch="CUTBLOCKS", dataName="forestry", dataClass="cutblocks"
  )
  
  disturbanceList <- list(forestry = list(cutblocks = vect(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                                                  ncol=2, byrow=TRUE),
                                                           type="polygons", crs=crs(rtm))))
  
  out <- calculateRate(
    disturbanceParameters = copy(disturbanceParameters),
    disturbanceDT = disturbanceDT,
    disturbanceList = disturbanceList,
    whichToUpdate = integer(0),
    RTM = rtm,
    studyArea = sa,
    DisturbanceRate = NULL,
    totalDisturbanceRate = NULL,
    diffYears = "2010_2015",
    destinationPath = tempdir(),
    disturbanceRateRelatesToBufferedArea = FALSE,
    maskOutLinesFromPolys = FALSE,
    aggregateSameDisturbances = FALSE,
    archiveNEW = NULL, targetFileNEW = NULL, urlNEW = NULL,
    archiveOLD = NULL, targetFileOLD = NULL, urlOLD = NULL
  )
  
  expect_identical(out, disturbanceParameters)
})

test_that("DisturbanceRate: no matching row -> warning and rate set to 0 (not NA)", {
  skip_on_cran()
  
  rtm <- rast(nrows=2, ncols=2, xmin=0, xmax=20, ymin=0, ymax=20, vals=1); crs(rtm) <- "EPSG:3857"
  sa  <- as.polygons(ext(rtm)); crs(sa) <- crs(rtm)
  
  dParams <- data.table(
    dataName="forestry", dataClass="cutblocks", disturbanceType="Generating",
    disturbanceOrigin="cutblocks", disturbanceRate=NA_real_
  )
  
  dDT <- data.table(fieldToSearch="Class", classToSearch="CUTBLOCKS",
                    dataName="forestry", dataClass="cutblocks")
  
  # Provide a non-null lay so we reach the DisturbanceRate branch:
  dl <- list(forestry = list(cutblocks = vect(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                                     ncol=2, byrow=TRUE),
                                              type="polygons", crs="EPSG:3857")))
  
  DR <- data.table( # deliberately NOT matching (different class)
    dataName="forestry", dataClass="roads", disturbanceType="Generating",
    disturbanceOrigin="roads", disturbanceRate=1.23
  )
  
  expect_warning(
    out <- calculateRate(dParams, dDT, dl, whichToUpdate = 1L, RTM = rtm, diffYears="2010_2015",
                         destinationPath=tempdir(), studyArea=sa,
                         DisturbanceRate = DR, totalDisturbanceRate = NULL,
                         disturbanceRateRelatesToBufferedArea=FALSE,
                         maskOutLinesFromPolys=FALSE, aggregateSameDisturbances=FALSE,
                         archiveNEW=NULL, targetFileNEW=NULL, urlNEW=NULL,
                         archiveOLD=NULL, targetFileOLD=NULL, urlOLD=NULL),
    regexp = "No matching row in DisturbanceRate", ignore.case = TRUE
  )
  expect_equal(out[1, disturbanceRate], 0)
})

test_that("DisturbanceRate: multiple matches -> warning and first row is used", {
  skip_on_cran()
  
  rtm <- rast(nrows=2, ncols=2, xmin=0, xmax=20, ymin=0, ymax=20, vals=1); crs(rtm) <- "EPSG:3857"
  sa  <- as.polygons(ext(rtm)); crs(sa) <- crs(rtm)
  
  dParams <- data.table(
    dataName="forestry", dataClass="cutblocks", disturbanceType="Generating",
    disturbanceOrigin="cutblocks", disturbanceRate=NA_real_
  )
  dDT <- data.table(fieldToSearch="Class", classToSearch="CUTBLOCKS",
                    dataName="forestry", dataClass="cutblocks")
  
  dl <- list(forestry = list(cutblocks = vect(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                                     ncol=2, byrow=TRUE),
                                              type="polygons", crs="EPSG:3857")))
  
  DR <- data.table(
    dataName="forestry", dataClass="cutblocks", disturbanceType="Generating",
    disturbanceOrigin="cutblocks", disturbanceRate=c(0.5, 1.5) # first should win
  )
  
  expect_warning(
    out <- calculateRate(dParams, dDT, dl, whichToUpdate = 1L, RTM = rtm, diffYears="2010_2015",
                         destinationPath=tempdir(), studyArea=sa,
                         DisturbanceRate = DR, totalDisturbanceRate = NULL,
                         disturbanceRateRelatesToBufferedArea=FALSE,
                         maskOutLinesFromPolys=FALSE, aggregateSameDisturbances=FALSE,
                         archiveNEW=NULL, targetFileNEW=NULL, urlNEW=NULL,
                         archiveOLD=NULL, targetFileOLD=NULL, urlOLD=NULL),
    regexp = "Multiple matches", ignore.case = TRUE
  )
  expect_equal(out[1, disturbanceRate], 0.5)
})

test_that("ECCC: no rows for requested classes -> warns and sets all targeted rates to 0", {
  skip_on_cran()
  
  # Define stub_data: empty table for all requested classes
  stub_data <- data.table(
    Class = character(),
    yearOLD = numeric(),
    yearNEW = numeric()
  )
  
  # Save the original function
  original_fn <- disturbanceInfoFromECCC
  # Assign the stub
  assign("disturbanceInfoFromECCC", function(...) list(AD_changed = stub_data), envir = .GlobalEnv)
  # Ensure the original function is restored after the test
  on.exit(assign("disturbanceInfoFromECCC", original_fn, envir = .GlobalEnv), add = TRUE)
  
  rtm <- rast(nrows=2, ncols=2, xmin=0, xmax=20, ymin=0, ymax=20, vals=1); crs(rtm) <- "EPSG:3857"
  sa  <- as.polygons(ext(rtm)); crs(sa) <- crs(rtm)
  
  dParams <- data.table(
    dataName=c("forestry","mining"),
    dataClass=c("cutblocks","mines"),
    disturbanceType="Generating",
    disturbanceOrigin=c("cutblocks","mines"),
    disturbanceRate=NA_real_
  )
  dDT <- data.table(
    fieldToSearch="Class",
    classToSearch=c("CUTBLOCKS","MINES"),
    dataName=c("forestry","mining"),
    dataClass=c("cutblocks","mines")
  )
  dl <- list(
    forestry = list(cutblocks = vect(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), 2, byrow=TRUE),
                                     type="polygons", crs="EPSG:3857")),
    mining   = list(mines     = vect(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), 2, byrow=TRUE),
                                     type="polygons", crs="EPSG:3857"))
  )
  
  expect_warning(
    out <- calculateRate(
      disturbanceParameters = dParams,
      disturbanceDT = dDT,
      disturbanceList = dl,
      whichToUpdate = 1:2,
      RTM = rtm,
      diffYears="2010_2015",
      destinationPath=tempdir(),
      studyArea=sa,
      DisturbanceRate = NULL,
      totalDisturbanceRate = NULL,
      disturbanceRateRelatesToBufferedArea=FALSE,
      maskOutLinesFromPolys=FALSE,
      aggregateSameDisturbances=FALSE,
      archiveNEW=NULL, targetFileNEW=NULL, urlNEW=NULL,
      archiveOLD=NULL, targetFileOLD=NULL, urlOLD=NULL
    ),
    regexp = "No historical ECCC disturbance", ignore.case = TRUE
  )
  expect_equal(out[dataName=="forestry", disturbanceRate], 0)
  expect_equal(out[dataName=="mining",   disturbanceRate], 0)
})

test_that("ECCC: per-class missing row -> WARNING + sets only that row to 0 (document desired fix)", {
  skip_on_cran()
  
  # Define stub_data: only "CUTBLOCKS" present, "MINES" missing
  stub_data <- data.table(
    Class = c("CUTBLOCKS"),
    yearOLD = c(0),
    yearNEW = c(0)
  )
  
  # Save the original function
  original_fn <- disturbanceInfoFromECCC
  # Assign the stub
  assign("disturbanceInfoFromECCC", function(...) list(AD_changed = stub_data), envir = .GlobalEnv)
  # Ensure the original function is restored after the test
  on.exit(assign("disturbanceInfoFromECCC", original_fn, envir = .GlobalEnv), add = TRUE)
  
  rtm <- rast(nrows=2, ncols=2, xmin=0, xmax=20, ymin=0, ymax=20, vals=1); crs(rtm) <- "EPSG:3857"
  sa  <- as.polygons(ext(rtm)); crs(sa) <- crs(rtm)
  
  dParams <- data.table(
    dataName=c("forestry","mining"),
    dataClass=c("cutblocks","mines"),
    disturbanceType="Generating",
    disturbanceOrigin=c("cutblocks","mines"),
    disturbanceRate=NA_real_
  )
  dDT <- data.table(
    fieldToSearch="Class",
    classToSearch=c("CUTBLOCKS","MINES"),
    dataName=c("forestry","mining"),
    dataClass=c("cutblocks","mines")
  )
  dl <- list(
    forestry = list(cutblocks = vect(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), 2, byrow=TRUE),
                                     type="polygons", crs="EPSG:3857")),
    mining   = list(mines     = vect(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), 2, byrow=TRUE),
                                     type="polygons", crs="EPSG:3857"))
  )
  
  expect_warning(
    out <- calculateRate(
      disturbanceParameters = dParams,
      disturbanceDT = dDT,
      disturbanceList = dl,
      whichToUpdate = 1:2,
      RTM = rtm,
      diffYears="2010_2015",
      destinationPath=tempdir(),
      studyArea=sa,
      DisturbanceRate = NULL,
      totalDisturbanceRate = NULL,
      disturbanceRateRelatesToBufferedArea=FALSE,
      maskOutLinesFromPolys=FALSE,
      aggregateSameDisturbances=FALSE,
      archiveNEW=NULL, targetFileNEW=NULL, urlNEW=NULL,
      archiveOLD=NULL, targetFileOLD=NULL, urlOLD=NULL
    ),
    regexp = "No ECCC-derived rate|Using yearly disturbance rate|Repeated data", ignore.case = TRUE
  )
  expect_true(out[dataName=="forestry", disturbanceRate] >= 0)
  expect_equal(out[dataName=="mining",  disturbanceRate], 0)
})


# -------------------------------------------------------------------
# 23. Test that real‐fetch flags get to prepInputs (integration stub)
#test_that("passes fetch flags through to disturbanceInfoFromECCC", {
#  originalECCC <- get("disturbanceInfoFromECCC", envir = .GlobalEnv)
#  captured <- NULL
#  assign("disturbanceInfoFromECCC", function(...){
#    captured <<- list(...)
#    # return minimal AD_changed so calculateRate_new continues
#    list(AD_changed = data.table(Class="C1", yearOLD=0, yearNEW=0))
#  }, envir = .GlobalEnv)
#  on.exit(assign("disturbanceInfoFromECCC", originalECCC, envir = .GlobalEnv), add = TRUE)
#  
#  calculateRate(
#    disturbanceParameters           = params_template(),
#    disturbanceDT                   = disturbanceDT,
#    disturbanceList                 = list(dn1 = list(orig1 = dummy_layer)),
#    whichToUpdate                   = 1L,
#    RTM                             = rast_template(),
#    studyArea                       = create_study_area(),
#    destinationPath                 = tempdir(),
#    overwriteDisturbanceLayersNEW   = TRUE,
#    overwriteDisturbanceLayersOLD   = FALSE,
#    diffYears                       = "2010_2015",
#    archiveNEW                      = "archiveNEW.zip",
#    targetFileNEW                   = "new.shp",
#    urlNEW                          = "http://new",
#    archiveOLD                      = "archiveOLD.zip",
#    targetFileOLD                   = "old.shp",
#    urlOLD                          = "http://old",
#    DisturbanceRate                 = NULL,
#    totalDisturbanceRate            = NULL,
#    disturbanceRateRelatesToBufferedArea = FALSE,
#    maskOutLinesFromPolys               = FALSE,
#    aggregateSameDisturbances           = FALSE
#  )
#  
#  expect_true(captured$overwriteDisturbanceLayersNEW, info="overwriteDisturbanceLayersNEW not forwarded")
#  expect_equal(captured$archiveNEW,     "archiveNEW.zip")
#  expect_equal(captured$targetFileNEW,  "new.shp")
#  expect_equal(captured$urlNEW,         "http://new")
#  # similarly for OLD if you want
#})


# Further tests for totalDisturbanceRate and empirical ECCC fallback would require mocking
# disturbanceInfoFromECCC() or pre-populating the cache files.
# These are left as integration tests rather than pure unit tests.
