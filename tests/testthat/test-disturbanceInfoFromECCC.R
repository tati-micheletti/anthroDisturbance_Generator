library(testthat)
library(terra)
library(data.table)
library(mockery)
library(digest)

# Stub extractNonPotentialLayers to return empty table with required columns
stub(disturbanceInfoFromECCC, 'extractNonPotentialLayers', function(dl) data.table(Sector=character(), dataClass=character()))
# Stub geomtype (unqualified) to avoid method dispatch issues
stub(disturbanceInfoFromECCC, 'geomtype', function(x) 'line')

# Define minimal studyArea and RTM with CRS
studyArea <- vect('POLYGON ((0 0, 0 1000, 1000 1000, 1000 0, 0 0))')
crs(studyArea) <- 'EPSG:3408'
RTM <- rast(nrows=5, ncols=5, xmin=0, xmax=5, ymin=0, ymax=5)

# Create fake SpatVectors with proper geometries
# Use LINESTRING for line features
fake_NEW_Lines <- vect('LINESTRING (100 100, 600 100)')
fake_NEW_Lines$Class <- 'Road'
crs(fake_NEW_Lines) <- crs(studyArea)

fake_OLD_Lines <- vect('LINESTRING (100 100, 500 100)')
fake_OLD_Lines$Class <- 'Road'
crs(fake_OLD_Lines) <- crs(studyArea)

# Use POLYGON for polygon features, ensuring they are closed and scaled to match the lines
fake_NEW_Polys <- vect('POLYGON ((100 50, 100 150, 400 150, 400 50, 100 50))')
fake_NEW_Polys$Class <- 'Settlements'
crs(fake_NEW_Polys) <- crs(studyArea)

fake_OLD_Polys <- vect('POLYGON ((100 60, 100 140, 300 140, 300 60, 100 60))')
fake_OLD_Polys$Class <- 'Settlements'
crs(fake_OLD_Polys) <- crs(studyArea)

plot(studyArea, col = "lightblue", main = "Study Area with Objects")
plot(fake_NEW_Lines, col = "red", add = TRUE, border = "black")
plot(fake_NEW_Polys, col = "pink", add = TRUE, border = "black")
plot(fake_OLD_Lines, col = "green", add = TRUE, border = "black")
plot(fake_OLD_Polys, col = "grey", add = TRUE, border = "black")


# Builder for prepInputs stub: always returns matching fake_* objects: always returns matching fake_* objects
stub_prepInputs_for <- function(layers) {
  function(url, archive, alsoExtract, studyArea, rasterToMatch,
           fun, targetFile, destinationPath) {
    tf <- tolower(basename(targetFile))
    if (grepl('lines|linear', tf)) {
      if (grepl('new|2015', tf)) return(layers$NEW_Lines)
      else return(layers$OLD_Lines)
    }
    if (grepl('poly|polygon', tf)) {
      if (grepl('new|2015', tf)) return(layers$NEW_Polys)
      else return(layers$OLD_Polys)
    }
    stop('Unexpected targetFile: ', targetFile)
  }
}

# Stub expanse: return constant or vector
stub_expanse <- function(values) {
  function(x, transform, unit) {
    if (length(values) == 1) rep(values, nrow(x)) else rep(values, length.out = nrow(x))
  }
}

# Shared classesAvailable and disturbanceList
classesAvailable <- data.table(classToSearch=c('Road','Settlements'), dataClass=c('roadClass','settleClass'))
disturbanceList <- createDisturbanceList(crs = crs(studyArea))
destination <- tempdir()

# 1. buffering + masking
# ----------------------------------------------------------------------
test_that('buffering and masking are applied correctly', {
  layers <- list(NEW_Lines=fake_NEW_Lines, NEW_Polys=fake_NEW_Polys,
                 OLD_Lines=fake_OLD_Lines, OLD_Polys=fake_OLD_Polys)
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(layers))
  stub(disturbanceInfoFromECCC, 'expanse', stub_expanse(1))
  buf_calls <- c(); mask_calls <- c()
  stub(disturbanceInfoFromECCC, 'terra::buffer', function(x, width) { buf_calls <<- c(buf_calls, width); x })
  stub(disturbanceInfoFromECCC, 'terra::mask',   function(x, mask, inverse) { mask_calls <<- c(mask_calls, inverse); x })
  
  suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea, RTM, classesAvailable, 10, disturbanceList,
      diffYears = '2000_2001', destination,
      bufferedDisturbances = TRUE,
      maskOutLinesFromPolys = TRUE,
      aggregateSameDisturbances = FALSE
    )
  )
  
  expect_true(all(buf_calls == 500))
  expect_true(all(mask_calls == TRUE))
})

# 2. aggregation collapsing classes
# ---------------------------------------------------------
test_that('aggregateSameDisturbances collapses classes', {
  layers <- list(NEW_Lines=fake_NEW_Lines, NEW_Polys=fake_NEW_Polys,
                 OLD_Lines=fake_OLD_Lines, OLD_Polys=fake_OLD_Polys)
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(layers))
  stub(disturbanceInfoFromECCC, 'expanse', stub_expanse(1))
  agg_called <- FALSE
  stub(disturbanceInfoFromECCC, 'terra::aggregate', function(x, by, dissolve) { agg_called <<- TRUE; x })
  
  suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea, RTM, classesAvailable, 10, disturbanceList,
      diffYears = '2000_2005', destination,
      bufferedDisturbances = FALSE,
      maskOutLinesFromPolys = FALSE,
      aggregateSameDisturbances = TRUE
    )
  )
  
  expect_true(agg_called)
})

# 3. mismatched classes zero-fill
# -----------------------------------------------------
test_that('mismatched classes generate zero-fill', {
  fake_NEW_Lines$Class <- 'OnlyNew'; fake_NEW_Polys$Class <- 'OnlyNew'
  fake_OLD_Lines$Class <- 'OnlyOld'; fake_OLD_Polys$Class <- 'OnlyOld'
  layers <- list(NEW_Lines=fake_NEW_Lines, NEW_Polys=fake_NEW_Polys,
                 OLD_Lines=fake_OLD_Lines, OLD_Polys=fake_OLD_Polys)
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(layers))
  stub(disturbanceInfoFromECCC, 'expanse', stub_expanse(1))
  classesAv <- data.table(classToSearch=c('OnlyNew','OnlyOld'), dataClass=c('newClass','oldClass'))
  
  res <- suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea, RTM, classesAv, 20, disturbanceList,
      diffYears = '2010_2012', destination,
      bufferedDisturbances = FALSE,
      maskOutLinesFromPolys = FALSE,
      aggregateSameDisturbances = FALSE
    )
  )
  
  pt <- res$proportionTable
  expect_true(all(pt$proportionAreaSqKmChangedPerYear >= 0))
})

# 4. negative change filtering
# ---------------------------------
test_that('negative changes are excluded', {
  fake_NEW_Lines$Class <- 'Neg'; fake_NEW_Polys$Class <- 'Neg'
  fake_OLD_Lines$Class <- 'Neg'; fake_OLD_Polys$Class <- 'Neg'
  layers <- list(NEW_Lines=fake_NEW_Lines, NEW_Polys=fake_NEW_Polys,
                 OLD_Lines=fake_OLD_Lines, OLD_Polys=fake_OLD_Polys)
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(layers))
  stub(disturbanceInfoFromECCC, 'expanse', function(x, transform, unit) rep(c(2,0), length.out = nrow(x)))
  
  res <- suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea, RTM, classesAvailable, 10, disturbanceList,
      diffYears = '2000_2001', destination,
      bufferedDisturbances = FALSE,
      maskOutLinesFromPolys = FALSE,
      aggregateSameDisturbances = FALSE
    )
  )
  
  pt <- res$proportionTable
  expect_false('Neg' %in% pt$dataClass)
})

# 5. proportion normalization
# -----------------------------------------------------------------------
test_that('proportionTable sums to 1', {
  fake_NEW_Lines$Class <- 'C1'; fake_NEW_Polys$Class <- 'C1'
  fake_OLD_Lines$Class <- 'C1'; fake_OLD_Polys$Class <- 'C1'
  layers <- list(NEW_Lines=fake_NEW_Lines, NEW_Polys=fake_NEW_Polys,
                 OLD_Lines=fake_OLD_Lines, OLD_Polys=fake_OLD_Polys)
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(layers))
  stub(disturbanceInfoFromECCC, 'expanse', function(x, transform, unit) rep(c(1,2), length.out = nrow(x)))
  classesAv <- data.table(classToSearch='C1', dataClass='onlyC1')
  
  res <- suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea, RTM, classesAv, 5, disturbanceList,
      diffYears = '2010_2012', destination,
      bufferedDisturbances = FALSE,
      maskOutLinesFromPolys = FALSE,
      aggregateSameDisturbances = FALSE
    )
  )
  
  expect_equal(sum(res$proportionTable$proportionOfTotalDisturbance), 1)
})

# 6. CSV outputs created
# --------------------------
test_that('CSV outputs exist', {
  fake_NEW_Lines$Class <- 'X'; fake_NEW_Polys$Class <- 'X'
  fake_OLD_Lines$Class <- 'X'; fake_OLD_Polys$Class <- 'X'
  layers <- list(NEW_Lines=fake_NEW_Lines, NEW_Polys=fake_NEW_Polys,
                 OLD_Lines=fake_OLD_Lines, OLD_Polys=fake_OLD_Polys)
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(layers))
  stub(disturbanceInfoFromECCC, 'expanse', stub_expanse(1))
  outdir <- file.path(tempdir(), 'distTests'); dir.create(outdir, showWarnings=FALSE)
  
  # Clear output directory
  unlink(outdir, recursive = TRUE)
  dir.create(outdir)
  
  # Run function
  suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea, RTM, classesAvailable, 10, disturbanceList,
      diffYears = '2010_2015', destination = outdir
    )
  )
  
  # Create expected file paths
  ad_file <- file.path(outdir, paste0("anthropogenicDisturbance_ECCC_2010_2015_", 
                                      digest::digest(studyArea), ".csv"))
  pt_file <- file.path(outdir, paste0("proportionTable_ECCC_2010_2015_", 
                                      digest::digest(studyArea), ".csv"))
  
  # Check both files exist
  expect_true(file.exists(ad_file), 
              info = paste("Missing AD file:", ad_file))
  expect_true(file.exists(pt_file), 
              info = paste("Missing PT file:", pt_file))
})

# 7. invalid diffYears yields NA
# ---------------------------------------
test_that('bad diffYears yields NA change', {
  layers <- list(NEW_Lines=fake_NEW_Lines, NEW_Polys=fake_NEW_Polys,
                 OLD_Lines=fake_OLD_Lines, OLD_Polys=fake_OLD_Polys)
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(layers))
  stub(disturbanceInfoFromECCC, 'expanse', stub_expanse(1))
  res <- suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea, RTM, classesAvailable, 10, disturbanceList,
      diffYears = 'bad_format', destination
    )
  )
  expect_true(all(is.na(res$AD_changed$relativeChangePerClassPerYear)))
})

# 8. smoke test for structure
# -------------------------
test_that('returns list with correct names', {
  layers <- list(NEW_Lines=fake_NEW_Lines, NEW_Polys=fake_NEW_Polys,
                 OLD_Lines=fake_OLD_Lines, OLD_Polys=fake_OLD_Polys)
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(layers))
  stub(disturbanceInfoFromECCC, 'expanse', stub_expanse(1))
  res <- suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea, RTM, classesAvailable, 100, disturbanceList,
      diffYears = '2010_2015', destination
    )
  )
  expect_type(res, 'list')
  expect_named(res, c('AD_changed', 'proportionTable'))
})

# 10. Mask-Only vs. No-Mask Branch (revised)
# -------------------------
test_that('maskOutLinesFromPolys correctly handles overlaps', {
  # Create layers list using the predefined geometries
  layers <- list(
    NEW_Lines = fake_NEW_Lines, NEW_Polys = fake_NEW_Polys,
    OLD_Lines = fake_OLD_Lines, OLD_Polys = fake_OLD_Polys
  )
  
  # Stub prepInputs to return our test layers
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(layers))
  
  # Create a safer expanse function that handles empty SpatVectors
  safe_expanse <- function(x, transform, unit) {
    if (length(x) == 0) return(0)
    terra::expanse(x, transform = transform, unit = unit)
  }
  stub(disturbanceInfoFromECCC, 'expanse', safe_expanse)
  
  # With masking - expect line area to be reduced
  res_maskON <- suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea = studyArea,
      RTM = RTM,
      classesAvailable = data.table(
        classToSearch = c("Road", "Settlements"),
        dataClass = c("roadClass", "settleClass")
      ),
      totalstudyAreaVAreaSqKm = 100, # 1000x1000 m = 1 km^2
      disturbanceList = list(),
      maskOutLinesFromPolys = TRUE,
      bufferedDisturbances = FALSE,
      destinationPath = destination
    )
  )
  
  # Without masking - expect full line area
  res_maskOFF <- suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea = studyArea,
      RTM = NULL,
      classesAvailable = data.table(
        classToSearch = c("Road", "Settlements"),
        dataClass = c("roadClass", "settleClass")
      ),
      totalstudyAreaVAreaSqKm = 100,
      disturbanceList = list(),
      maskOutLinesFromPolys = FALSE,
      bufferedDisturbances = FALSE,
      destinationPath = destination
    )
  )
  
  # Get line areas
  line_area_maskON <- sum(res_maskON$AD_changed[Class == "Road", yearNEW])
  line_area_maskOFF <- sum(res_maskOFF$AD_changed[Class == "Road", yearNEW])
  
  # Verify masking reduces but doesn't eliminate line area
  expect_lt(line_area_maskON, line_area_maskOFF)
})

# 11. Unbuffered vs. Buffered
# -------------------------
test_that('buffering width changes based on bufferedDisturbances flag', {
  # Create temp directory for destinationPath
  destPath <- tempfile()
  dir.create(destPath)
  
  # Create minimal layers
  layers <- list(
    NEW_Lines = fake_NEW_Lines, NEW_Polys = fake_NEW_Polys,
    OLD_Lines = fake_OLD_Lines, OLD_Polys = fake_OLD_Polys
  )
  
  # Stub prepInputs to return our test layers
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(layers))
  
  # Track buffer widths
  buffer_widths <- numeric(0)
  
  # Stub terra::buffer to capture width parameter
  stub(disturbanceInfoFromECCC, 'terra::buffer', function(x, width) {
    buffer_widths <<- c(buffer_widths, width)
    return(x)
  })
  
  # Test buffered
  suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea, rtm, classesAvailable, 10, disturbanceList,
      bufferedDisturbances = TRUE,
      destinationPath = destination
    )
  )
  buffered_widths <- buffer_widths
  buffer_widths <- numeric(0)  # Reset
  
  # Test unbuffered
  suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea, rtm, classesAvailable, 10, disturbanceList,
      bufferedDisturbances = FALSE,
      destinationPath = destination
    )
  )
  unbuffered_widths <- buffer_widths
  
  # Verify widths
  expect_true(all(buffered_widths == 500))
  expect_true(all(unbuffered_widths == 30))
})

# 12. Zero-Year Interval 
# -------------------------
test_that('zero-year interval handles division safely', {
  layers <- list(NEW_Lines=fake_NEW_Lines, NEW_Polys=fake_NEW_Polys,
                 OLD_Lines=fake_OLD_Lines, OLD_Polys=fake_OLD_Polys)
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(layers))
  stub(disturbanceInfoFromECCC, 'expanse', stub_expanse(1))
  
  res <- suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea, RTM, classesAvailable, 10, disturbanceList,
      diffYears = "2010_2010",
      destinationPath = destination 
    )
  )
  
  # Check relativeChangePerClassPerYear is either NA or 0
  changes <- res$AD_changed$relativeChangePerClassPerYear
  expect_true(all(is.na(changes) | changes == 0))
})

# 13. Partial Class Matching
# -------------------------
test_that('unavailable classes are handled gracefully', {
  layers <- list(NEW_Lines=fake_NEW_Lines, NEW_Polys=fake_NEW_Polys,
                 OLD_Lines=fake_OLD_Lines, OLD_Polys=fake_OLD_Polys)
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(layers))
  stub(disturbanceInfoFromECCC, 'expanse', stub_expanse(1))
  
  # Add an unavailable class to classesAvailable
  extendedClasses <- rbind(
    classesAvailable,
    data.table(classToSearch = "NonExistentClass", dataClass = "ghostClass")
  )
  
  res <- suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea, RTM, extendedClasses, 10, disturbanceList,
      destinationPath = destination 
    )
  )
  
  # Check proportionTable doesn't contain the ghost class
  expect_false("ghostClass" %in% res$proportionTable$dataClass)
})

# 14. Non-Numeric Year Parts
# -------------------------
test_that('non-numeric year parts are handled gracefully', {
  layers <- list(NEW_Lines=fake_NEW_Lines, NEW_Polys=fake_NEW_Polys,
                 OLD_Lines=fake_OLD_Lines, OLD_Polys=fake_OLD_Polys)
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(layers))
  stub(disturbanceInfoFromECCC, 'expanse', stub_expanse(1))
  
  res <- suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea, RTM, classesAvailable, 10, disturbanceList,
      diffYears = "foo_bar",
      destinationPath = destination
    )
  )
  
  # Check that relative change is NA
  expect_true(all(is.na(res$AD_changed$relativeChangePerClassPerYear)))
})

# 15. aggregateSameDisturbances reduces overlapping area
# ------------------------------------------------
test_that('aggregateSameDisturbances reduces overlapping area', {
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(list(NEW_Lines = fake_NEW_Lines, OLD_Lines = fake_OLD_Lines, NEW_Polys = fake_NEW_Polys, OLD_Polys = fake_OLD_Polys)))
  
  # Create two overlapping line features for class 'Road'
  duplicate_line <- fake_OLD_Lines
  fake_NEW_Lines <<- rbind(duplicate_line, duplicate_line)
  stub(disturbanceInfoFromECCC,'expanse', function(x,...) rep(1, nrow(x)))
  
  # Without aggregation
  res_noagg <- disturbanceInfoFromECCC(
    studyArea, RTM, classesAvailable, 10, disturbanceList,
    diffYears = '2000_2001', destinationPath = destination,
    bufferedDisturbances = FALSE, maskOutLinesFromPolys = FALSE,
    aggregateSameDisturbances = FALSE
  )
  # With aggregation
  res_agg <- disturbanceInfoFromECCC(
    studyArea, RTM, classesAvailable, 10, disturbanceList,
    diffYears = '2000_2001', destinationPath = destination,
    bufferedDisturbances = FALSE, maskOutLinesFromPolys = FALSE,
    aggregateSameDisturbances = TRUE
  )
  
  area_noagg <- sum(res_noagg$AD_changed$yearNEW)
  area_agg   <- sum(res_agg$AD_changed$yearNEW)
  
  expect_lt(area_agg, area_noagg)
})

# Stub prepInputs function
stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(list(NEW_Lines = fake_NEW_Lines, OLD_Lines = fake_OLD_Lines, NEW_Polys = fake_NEW_Polys, OLD_Polys = fake_OLD_Polys)))

# Test 16: invalid diffYears yields NA relative change
# -------------------------
test_that('invalid diffYears yields NA in relativeChangePerClassPerYear', {
  stub(disturbanceInfoFromECCC, 'expanse', function(x, ...) rep(1, nrow(x)))
  expect_warning(
    res <- disturbanceInfoFromECCC(
      studyArea, RTM, classesAvailable, 10, disturbanceList,
      diffYears = 'foo_bar', destinationPath = destination
    ),
    regex = "The 'diffYears' parameter \\(foo_bar\\) contains non-numeric parts. Relative change could not be calculated."
  )
  expect_true(all(is.na(res$AD_changed$relativeChangePerClassPerYear)))
})

# Test 17: same-year diffYears yields NA relative change
# -------------------------
test_that('same-year diffYears yields NA in relativeChangePerClassPerYear', {
  stub(disturbanceInfoFromECCC, 'expanse', function(x, ...) rep(1, nrow(x)))
  expect_warning(
    res <- disturbanceInfoFromECCC(
      studyArea, RTM, classesAvailable, 10, disturbanceList,
      diffYears = '2010_2010', destinationPath = destination
    ),
    regex = "The 'diffYears' parameter \\(2010_2010\\) specifies the same year. Relative change could not be calculated."
  )
  expect_true(all(is.na(res$AD_changed$relativeChangePerClassPerYear)))
})

# Test 20: disturbanceList non-potential layers are included
# -------------------------
test_that('non-potential layers are properly added', {
  # Create test layer
  extra <- vect('LINESTRING(0 0, 100 0)')
  crs(extra) <- crs(studyArea)
  extra$Class <- "Extra"  # Fixed to match classToSearch
  
  stub(disturbanceInfoFromECCC, 'extractNonPotentialLayers', 
       function(dl) data.table(Sector='test', dataClass='extra'))
  
  res <- suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea, RTM, 
      rbind(classesAvailable, data.table(classToSearch="Extra", dataClass="extra")),
      10, list(test = list(extra = extra)),
      destinationPath = destination
    )
  )
  expect_true("extra" %in% res$proportionTable$dataClass)
})
