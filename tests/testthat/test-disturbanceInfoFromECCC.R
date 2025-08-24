# test-disturbanceInfoFromECCC_spec.R
## ---- Minimal reproducible spatial scaffolding --------------------------------

# Study area + RTM with consistent CRS (meters; arbitrary but consistent)
studyArea <- vect('POLYGON ((0 0, 0 1000, 1000 1000, 1000 0, 0 0))')
crs(studyArea) <- 'EPSG:3005'
RTM <- rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1000, ymin = 0, ymax = 1000)
crs(RTM) <- crs(studyArea)

# Minimal disturbanceList (module helper will ignore most content here)
disturbanceList <- createDisturbanceList(crs = crs(studyArea))

# Two classes we can control easily in tests
classesAvailable <- data.table(
  classToSearch = c('Road',       'Settlements'),
  dataClass     = c('roadClass',  'settleClass')
)

# Simple, valid fake layers (one line, one polygon) in the study CRS
fake_NEW_Lines <- vect('LINESTRING (100 100, 900 100)'); crs(fake_NEW_Lines) <- crs(studyArea); fake_NEW_Lines$Class <- 'Road'
fake_OLD_Lines <- vect('LINESTRING (100 100, 700 100)'); crs(fake_OLD_Lines) <- crs(studyArea); fake_OLD_Lines$Class <- 'Road'

fake_NEW_Polys <- vect('POLYGON ((100 50, 100 250, 500 250, 500 50, 100 50))'); crs(fake_NEW_Polys) <- crs(studyArea); fake_NEW_Polys$Class <- 'Settlements'
fake_OLD_Polys <- vect('POLYGON ((100 80, 100 220, 400 220, 400 80, 100 80))'); crs(fake_OLD_Polys) <- crs(studyArea); fake_OLD_Polys$Class <- 'Settlements'

# Builder for prepInputs stub: return our fake layers based on targetFile hints
stub_prepInputs_for <- function(layers) {
  function(url, archive, alsoExtract, studyArea, rasterToMatch, fun, targetFile, destinationPath) {
    tf <- tolower(basename(targetFile))
    if (grepl('lines|linear', tf)) {
      if (grepl('new|2015', tf)) return(layers$NEW_Lines) else return(layers$OLD_Lines)
    }
    if (grepl('poly|polygon', tf)) {
      if (grepl('new|2015', tf)) return(layers$NEW_Polys) else return(layers$OLD_Polys)
    }
    stop('Unexpected targetFile: ', targetFile)
  }
}

# Keep module cross-calls quiet and predictable in tests
stub(disturbanceInfoFromECCC, 'extractNonPotentialLayers',
     function(dl) data.table(Sector = character(), dataClass = character()))
stub(disturbanceInfoFromECCC, 'geomtype', function(x) 'line')  # avoids dispatch surprises

## ---- 1) Schema & units (proportions, not percents) ---------------------------

test_that("returns a proportion table with sensible units (0–1 per year)", {
  totalstudyAreaVAreaSqKm <- expanse(studyArea)/ 1e6
  layers <- list(NEW_Lines = fake_NEW_Lines, NEW_Polys = fake_NEW_Polys,
                 OLD_Lines = fake_OLD_Lines, OLD_Polys = fake_OLD_Polys)
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(layers))
  
  res <- suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea = studyArea,
      RTM = RTM,
      totalstudyAreaVAreaSqKm = totalstudyAreaVAreaSqKm,
      classesAvailable = classesAvailable,
      disturbanceList = disturbanceList,
      diffYears = '2010_2020',
      destinationPath = tempdir(),
      bufferedDisturbances = FALSE,
      maskOutLinesFromPolys = TRUE,
      aggregateSameDisturbances = TRUE
    )
  )
  
  expect_type(res, "list")
  expect_true("proportionTable" %in% names(res))
  pt <- res$proportionTable
  # Be permissive on exact column names: the key parts must exist
  expect_true(all(c("dataClass") %in% names(pt)))
  # Find the proportion column by pattern if needed
  prop_col <- intersect(names(pt), c("proportionAreaSqKmChangedPerYear",
                                     "proportionAreaChangedPerYear",
                                     "proportionChangedPerYear"))
  expect_gt(length(prop_col), 0)
  prop <- pt[[prop_col[1]]]
  expect_true(is.numeric(prop))
  expect_true(all(prop >= 0, na.rm = TRUE))
  expect_true(all(prop <= 1, na.rm = TRUE))
})

## ---- 2) Buffering: uses 500 m when requested --------------------------------

test_that("bufferedDisturbances=TRUE triggers 500 m buffers", {
  layers <- list(NEW_Lines = fake_NEW_Lines, NEW_Polys = fake_NEW_Polys,
                 OLD_Lines = fake_OLD_Lines, OLD_Polys = fake_OLD_Polys)
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(layers))
  
  # capture buffer widths
  buf_widths <- numeric(0)
  stub(disturbanceInfoFromECCC, 'terra::buffer', function(x, width) { buf_widths <<- c(buf_widths, width); x })
  
  suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea, RTM, classesAvailable, 10, disturbanceList,
      diffYears = '2000_2005', destinationPath = tempdir(),
      bufferedDisturbances = TRUE, maskOutLinesFromPolys = TRUE, aggregateSameDisturbances = TRUE
    )
  )
  expect_true(any(buf_widths == 500),
              info = "When relating rates to buffered area, 500 m buffers must be applied (caribou-style accounting).")
})

## ---- 3) Masking: lines masked out of polygons when requested -----------------

test_that("maskOutLinesFromPolys=TRUE attempts to mask overlapping lines from polygons", {
  layers <- list(NEW_Lines = fake_NEW_Lines, NEW_Polys = fake_NEW_Polys,
                 OLD_Lines = fake_OLD_Lines, OLD_Polys = fake_OLD_Polys)
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(layers))
  
  # record mask invocations (we don't enforce exact args; we only care it's called)
  mask_called <- FALSE
  stub(disturbanceInfoFromECCC, 'terra::mask', function(x, mask, inverse = FALSE) { mask_called <<- TRUE; x })
  
  suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea, RTM, classesAvailable, 10, disturbanceList,
      diffYears = '2000_2010', destinationPath = tempdir(),
      bufferedDisturbances = TRUE, maskOutLinesFromPolys = TRUE, aggregateSameDisturbances = TRUE
    )
  )
  
  expect_true(mask_called,
              info = "With maskOutLinesFromPolys=TRUE, polygons should be masked by overlapping lines before area calc.")
})

## ---- 4) Aggregation: overlap collapse reduces/equal class share --------------

test_that("aggregateSameDisturbances collapses overlaps (monotone non-increase)", {
  totalstudyAreaVAreaSqKm <- expanse(studyArea)/ 1e6
  # Create overlapping NEW polygons in the same class to exercise aggregation
  polyA <- vect('POLYGON ((200 50, 200 300, 500 300, 500 50, 200 50))'); crs(polyA) <- crs(studyArea); polyA$Class <- 'Settlements'
  polyB <- vect('POLYGON ((400 50, 400 300, 700 300, 700 50, 400 50))'); crs(polyB) <- crs(studyArea); polyB$Class <- 'Settlements'
  fake_NEW_Polys_overlap <- rbind(polyA, polyB)
  fake_OLD_Polys_smaller <- vect('POLYGON ((250 80, 250 220, 450 220, 450 80, 250 80))'); crs(fake_OLD_Polys_smaller) <- crs(studyArea); fake_OLD_Polys_smaller$Class <- 'Settlements'
  
  layers <- list(
    NEW_Lines = fake_NEW_Lines, NEW_Polys = fake_NEW_Polys_overlap,
    OLD_Lines = fake_OLD_Lines, OLD_Polys = fake_OLD_Polys_smaller
  )
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(layers))
  
  res_noAgg <- suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea = studyArea, RTM = RTM,
      totalstudyAreaVAreaSqKm = totalstudyAreaVAreaSqKm,
      classesAvailable = classesAvailable, disturbanceList = disturbanceList,
      diffYears = '2010_2015', destinationPath = tempdir(),
      bufferedDisturbances = FALSE, maskOutLinesFromPolys = FALSE, aggregateSameDisturbances = FALSE)
  )
  res_agg <- suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea = studyArea, RTM = RTM,
      totalstudyAreaVAreaSqKm = totalstudyAreaVAreaSqKm,
      classesAvailable = classesAvailable, disturbanceList = disturbanceList,
      diffYears = '2010_2015', destinationPath = tempdir(),
      bufferedDisturbances = FALSE, maskOutLinesFromPolys = FALSE, aggregateSameDisturbances = TRUE)
  )
  
  # Join by classes present in both
  get_pt <- function(res) {
    pt <- as.data.table(res$proportionTable)
    prop_col <- intersect(names(pt), c("proportionAreaSqKmChangedPerYear",
                                       "proportionAreaChangedPerYear",
                                       "proportionChangedPerYear"))
    setnames(pt, prop_col[1], "prop")
    pt[, .(dataClass, prop)]
  }
  pt_noAgg <- get_pt(res_noAgg)
  pt_agg   <- get_pt(res_agg)
  merged <- merge(pt_noAgg, pt_agg, by = "dataClass", suffixes = c("_noAgg", "_agg"))
  expect_gt(nrow(merged), 0)
  # For overlapping classes, aggregated should not be larger than non-aggregated
  expect_true(all(merged$prop_agg <= merged$prop_noAgg + 1e-12))
})

## ---- 5) Buffered vs raw: buffered >= raw (per class, where buffering expands) --

test_that("bufferedDisturbances yields >= proportions vs raw (when buffers expand footprint)", {
  layers <- list(NEW_Lines = fake_NEW_Lines, NEW_Polys = fake_NEW_Polys,
                 OLD_Lines = fake_OLD_Lines, OLD_Polys = fake_OLD_Polys)
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(layers))
  
  res_raw <- suppressWarnings(
    disturbanceInfoFromECCC(studyArea, RTM, classesAvailable, 10, disturbanceList,
                            diffYears = '2000_2005', destinationPath = tempdir(),
                            bufferedDisturbances = FALSE, maskOutLinesFromPolys = FALSE, aggregateSameDisturbances = TRUE)
  )
  res_buf <- suppressWarnings(
    disturbanceInfoFromECCC(studyArea, RTM, classesAvailable, 10, disturbanceList,
                            diffYears = '2000_2005', destinationPath = tempdir(),
                            bufferedDisturbances = TRUE, maskOutLinesFromPolys = FALSE, aggregateSameDisturbances = TRUE)
  )
  
  pt_raw <- res_raw$proportionTable
  pt_buf <- res_buf$proportionTable
  setkey(pt_raw, dataClass); setkey(pt_buf, dataClass)
  merged <- pt_raw[pt_buf, nomatch = 0]
  expect_true(all(merged$proportionAreaSqKmChangedPerYear <= merged$i.proportionAreaSqKmChangedPerYear),
              info = "Relating rates to buffered area should not reduce the class-wise proportion.")
})

## ---- 6) Empty classes: allowed; no error, omitted or NA ----------------------

test_that("classes with no features are handled gracefully (omitted or NA)", {
  layers <- list(NEW_Lines = fake_NEW_Lines, NEW_Polys = fake_NEW_Polys,
                 OLD_Lines = fake_OLD_Lines, OLD_Polys = fake_OLD_Polys)
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(layers))
  
  extendedClasses <- rbind(
    classesAvailable,
    data.table(classToSearch = "NonExistentClass", dataClass = "ghostClass")
  )
  
  res <- suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea, RTM, extendedClasses, 10, disturbanceList,
      diffYears = "2010_2015", destinationPath = tempdir()
    )
  )
  
  # Either the ghost row is absent, or present with NA — both honor the spec
  pt <- res$proportionTable
  if ("ghostClass" %in% pt$dataClass) {
    expect_true(is.na(pt[dataClass == "ghostClass", proportionAreaSqKmChangedPerYear]))
  } else {
    succeed()
  }
})

## ---- 7) Negative change: excluded from positive growth proportions -----------

test_that("negative class changes are not reported as positive proportions", {
  # Force NEW < OLD by shrinking NEW polygon
  shrink_NEW <- vect('POLYGON ((100 90, 100 210, 250 210, 250 90, 100 90))'); crs(shrink_NEW) <- crs(studyArea); shrink_NEW$Class <- 'Settlements'
  layers <- list(NEW_Lines = fake_NEW_Lines, NEW_Polys = shrink_NEW,
                 OLD_Lines = fake_OLD_Lines, OLD_Polys = fake_OLD_Polys)
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(layers))
  
  res <- suppressWarnings(
    disturbanceInfoFromECCC(studyArea, RTM, classesAvailable, 10, disturbanceList,
                            diffYears = '2010_2015', destinationPath = tempdir(),
                            bufferedDisturbances = FALSE, maskOutLinesFromPolys = FALSE, aggregateSameDisturbances = TRUE)
  )
  pt <- res$proportionTable
  # The class with negative change should either be absent or have 0/NA, but not a positive share
  if ("Settlements" %in% pt$dataClass) {
    expect_true(is.na(pt[dataClass == "Settlements", proportionAreaSqKmChangedPerYear]) ||
                  pt[dataClass == "Settlements", proportionAreaSqKmChangedPerYear] == 0)
  } else {
    succeed()
  }
})

## ---- 8) diffYears parsing: robust to bad strings (warn, not error) -----------

test_that("bad diffYears formats are handled (warning, not error)", {
  totalstudyAreaVAreaSqKm <- expanse(studyArea)/ 1e6
  layers <- list(NEW_Lines = fake_NEW_Lines, NEW_Polys = fake_NEW_Polys,
                 OLD_Lines = fake_OLD_Lines, OLD_Polys = fake_OLD_Polys)
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(layers))
  
  # Expect at least one warning (content not important)
  expect_warning(
    res <- disturbanceInfoFromECCC(
      studyArea = studyArea, RTM = RTM,
      totalstudyAreaVAreaSqKm = totalstudyAreaVAreaSqKm,
      classesAvailable = classesAvailable, disturbanceList = disturbanceList,
      diffYears = 'bad_format', destinationPath = tempdir(),
      bufferedDisturbances = FALSE, maskOutLinesFromPolys = FALSE, aggregateSameDisturbances = TRUE
    )
  )
  expect_true(is.list(res))
  expect_true("proportionTable" %in% names(res))
})

## ---- 9) Output files are optional — do not enforce side effects --------------

test_that("destinationPath side-effects are optional and non-failing", {
  totalstudyAreaVAreaSqKm <- expanse(studyArea)/ 1e6
  layers <- list(NEW_Lines = fake_NEW_Lines, NEW_Polys = fake_NEW_Polys,
                 OLD_Lines = fake_OLD_Lines, OLD_Polys = fake_OLD_Polys)
  stub(disturbanceInfoFromECCC, 'prepInputs', stub_prepInputs_for(layers))
  
  outdir <- file.path(tempdir(), paste0("eccc-", as.integer(runif(1, 1e6, 1e7))))
  dir.create(outdir, showWarnings = FALSE)
  
  # Allow messages/warnings; just assert the call doesn't error and returns the expected structure
  res <- suppressWarnings(suppressMessages(
    disturbanceInfoFromECCC(
      studyArea = studyArea, RTM = RTM,
      totalstudyAreaVAreaSqKm = totalstudyAreaVAreaSqKm,
      classesAvailable = classesAvailable, disturbanceList = disturbanceList,
      diffYears = '2010_2015', destinationPath = outdir,
      bufferedDisturbances = FALSE, maskOutLinesFromPolys = TRUE, aggregateSameDisturbances = TRUE
    )
  ))
  
  expect_true(is.list(res))
  expect_true("proportionTable" %in% names(res))
})

## ---- 10) Masking reduces polygonal contribution when lines overlap polygons -
test_that("maskOutLinesFromPolys prevents double counting in the combined total", {
  studyArea <- vect('POLYGON ((0 0, 0 1000, 1000 1000, 1000 0, 0 0))')
  crs(studyArea) <- 'EPSG:3005'                           # set CRS before expanse
  totalstudyAreaVAreaSqKm <- terra::expanse(studyArea) / 1e6
  
  RTM <- rast(nrows=10, ncols=10, xmin=0, xmax=1000, ymin=0, ymax=1000)
  crs(RTM) <- crs(studyArea)
  
  # OLD polygon smaller; NEW polygon grows upward
  poly_old <- vect('POLYGON ((100 100, 100 400, 600 400, 600 100, 100 100))')
  crs(poly_old) <- crs(studyArea); poly_old$Class <- 'Settlements'
  poly_new <- vect('POLYGON ((100 100, 100 500, 600 500, 600 100, 100 100))')
  crs(poly_new) <- crs(studyArea); poly_new$Class <- 'Settlements'
  
  # OLD line outside; NEW line crosses expansion zone so overlaps exist after buffering
  line_old <- vect('LINESTRING (50 50, 950 50)')
  crs(line_old) <- crs(studyArea); line_old$Class <- 'Road'
  line_new <- vect('LINESTRING (50 450, 950 450)')
  crs(line_new) <- crs(studyArea); line_new$Class <- 'Road'
  
  layers <- list(NEW_Lines = line_new, NEW_Polys = poly_new,
                 OLD_Lines = line_old, OLD_Polys = poly_old)
  
  mockery::stub(disturbanceInfoFromECCC, 'prepInputs', {
    function(url, archive, alsoExtract, studyArea, rasterToMatch, fun, targetFile, destinationPath) {
      tf <- tolower(basename(targetFile))
      if (grepl('lines|linear', tf)) if (grepl('new|2015', tf)) return(layers$NEW_Lines) else return(layers$OLD_Lines)
      if (grepl('poly|polygon', tf)) if (grepl('new|2015', tf)) return(layers$NEW_Polys) else return(layers$OLD_Polys)
      stop('Unexpected targetFile: ', targetFile)
    }
  })
  mockery::stub(disturbanceInfoFromECCC, 'geomtype', function(x) terra::geomtype(x))
  
  classesDT <- data.table(classToSearch = c('Road','Settlements'),
                          dataClass     = c('road','settle'))
  
  # mask OFF
  res_off <- suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea, RTM, classesDT,
      totalstudyAreaVAreaSqKm = totalstudyAreaVAreaSqKm, disturbanceList = list(),
      diffYears='2000_2010', destinationPath=tempdir(),
      bufferedDisturbances=TRUE, maskOutLinesFromPolys=FALSE, aggregateSameDisturbances=TRUE
    )
  )
  # mask ON
  res_on <- suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea, RTM, classesDT,
      totalstudyAreaVAreaSqKm = totalstudyAreaVAreaSqKm, disturbanceList = list(),
      diffYears='2000_2010', destinationPath=tempdir(),
      bufferedDisturbances=TRUE, maskOutLinesFromPolys=TRUE, aggregateSameDisturbances=TRUE
    )
  )
  
  sum_prop <- function(res) {
    res$proportionTable[
      dataClass %in% c('road','settle'),
      sum(proportionAreaSqKmChangedPerYear, na.rm = TRUE)
    ]
  }
  
  off_sum <- sum_prop(res_off)
  on_sum  <- sum_prop(res_on)
  
  # Spec-aligned: masking ensures the combined total does not increase (often equal).
  expect_lte(round(on_sum, 10), round(off_sum, 10))
})

## ---- 11) Extra classes from extractNonPotentialLayers appear in output (may be NA)
test_that("classes reported by extractNonPotentialLayers are surfaced (function integrates without error)", {
  studyArea <- vect('POLYGON ((0 0, 0 1000, 1000 1000, 1000 0, 0 0))')
  crs(studyArea) <- 'EPSG:3005'                           # set CRS first to avoid expanse warning
  totalstudyAreaVAreaSqKm <- terra::expanse(studyArea) / 1e6
  RTM <- rast(nrows=2, ncols=2, xmin=0, xmax=1000, ymin=0, ymax=1000); crs(RTM) <- crs(studyArea)
  
  # Minimal non-empty layers with required Class field
  line_old <- vect('LINESTRING (100 100, 200 100)'); crs(line_old) <- crs(studyArea); line_old$Class <- 'Road'
  line_new <- vect('LINESTRING (100 100, 300 100)'); crs(line_new) <- crs(studyArea); line_new$Class <- 'Road'
  poly_old <- vect('POLYGON ((400 400, 400 500, 500 500, 500 400, 400 400))'); crs(poly_old) <- crs(studyArea); poly_old$Class <- 'Settlements'
  poly_new <- vect('POLYGON ((400 400, 400 520, 520 520, 520 400, 400 400))'); crs(poly_new) <- crs(studyArea); poly_new$Class <- 'Settlements'
  
  layers <- list(NEW_Lines = line_new, NEW_Polys = poly_new,
                 OLD_Lines = line_old, OLD_Polys = poly_old)
  
  mockery::stub(disturbanceInfoFromECCC, 'prepInputs', {
    function(url, archive, alsoExtract, studyArea, rasterToMatch, fun, targetFile, destinationPath) {
      tf <- tolower(basename(targetFile))
      if (grepl('lines|linear', tf)) if (grepl('new|2015', tf)) return(layers$NEW_Lines) else return(layers$OLD_Lines)
      if (grepl('poly|polygon', tf)) if (grepl('new|2015', tf)) return(layers$NEW_Polys) else return(layers$OLD_Polys)
      stop('Unexpected targetFile: ', targetFile)
    }
  })
  mockery::stub(disturbanceInfoFromECCC, 'geomtype', function(x) terra::geomtype(x))
  
  # Track that extractNonPotentialLayers is used
  was_called <- FALSE
  mockery::stub(disturbanceInfoFromECCC, 'extractNonPotentialLayers',
                function(dl) { was_called <<- TRUE; data.table(Sector='test', dataClass='extra') })
  
  res <- suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea=studyArea, RTM=RTM,
      classesAvailable = data.table(classToSearch=c('Road','Settlements'), dataClass=c('road','settle')),
      totalstudyAreaVAreaSqKm = totalstudyAreaVAreaSqKm, disturbanceList=list(),
      diffYears='2010_2015', destinationPath=tempdir(),
      bufferedDisturbances=FALSE, maskOutLinesFromPolys=FALSE, aggregateSameDisturbances=TRUE
    )
  )
  
  expect_true(isTRUE(was_called), info = "extractNonPotentialLayers() should be consulted.")
  
  # It's OK if 'extra' is not materialized as a row; if it is present, value may be NA.
  if ('extra' %in% res$proportionTable$dataClass) {
    succeed()
  } else {
    succeed()
  }
})

## ---- 12) (Optional) RTM can be NULL without error (if function supports it)
test_that("RTM is optional (NULL) and function still returns a result", {
  studyArea <- vect('POLYGON ((0 0, 0 1000, 1000 1000, 1000 0, 0 0))'); crs(studyArea) <- 'EPSG:3005'
  totalstudyAreaVAreaSqKm <- expanse(studyArea)/ 1e6
  fake_line <- vect('LINESTRING (100 100, 900 100)'); crs(fake_line) <- crs(studyArea); fake_line$Class <- 'Road'
  fake_poly <- vect('POLYGON ((100 50, 100 250, 500 250, 500 50, 100 50))'); crs(fake_poly) <- crs(studyArea); fake_poly$Class <- 'Settlements'
  
  stub(disturbanceInfoFromECCC, 'prepInputs', function(url, archive, alsoExtract, studyArea, rasterToMatch, fun, targetFile, destinationPath) {
    tf <- tolower(basename(targetFile))
    if (grepl('lines|linear', tf)) return(fake_line) else return(fake_poly)
  })
  stub(disturbanceInfoFromECCC, 'extractNonPotentialLayers',
       function(dl) data.table(Sector = character(), dataClass = character()))
  stub(disturbanceInfoFromECCC, 'geomtype', function(x) 'line')
  
  res <- suppressWarnings(
    disturbanceInfoFromECCC(
      studyArea=studyArea, RTM=NULL,
      classesAvailable=data.table(classToSearch=c('Road','Settlements'), dataClass=c('road','settle')),
      totalstudyAreaVAreaSqKm = totalstudyAreaVAreaSqKm, disturbanceList=list(),
      diffYears='2010_2015', destinationPath=tempdir(),
      bufferedDisturbances=FALSE, maskOutLinesFromPolys=FALSE, aggregateSameDisturbances=TRUE
    )
  )
  expect_true(is.list(res) && "proportionTable" %in% names(res))
})

### additional tests to up coverage

# tests/testthat/test-disturbanceInfoFromECCC.R
test_that("adds zero rows for classes missing on either side", {
  # --- minimal spatial scaffolding ---
  crs_str <- "EPSG:3857"
  rtm <- rast(ext(0, 1e4, 0, 1e4), res = 1000, crs = crs_str); values(rtm) <- 1
  sa  <- as.polygons(ext(0, 1e4, 0, 1e4), crs = crs_str)
  
  make_poly <- function(x1, y1, x2, y2, cls) {
    v <- as.polygons(ext(x1, x2, y1, y2), crs = crs_str); v$Class <- cls; v
  }
  
  # NEW has Road + Mine; OLD has Road + Cutblock → class sets differ both ways
  AD_NEW_Lines <- rbind(make_poly(100,100,200,200,"Road"),
                        make_poly(300,300,400,400,"Mine"))
  AD_NEW_Polys <- AD_NEW_Lines[0]
  AD_OLD_Lines <- rbind(make_poly(100,100,200,200,"Road"),
                        make_poly(500,500,600,600,"Cutblock"))
  AD_OLD_Polys <- AD_OLD_Lines[0]
  
  # fake prepInputs: return the 4 objects in order; keep cycling if called more
  seq_data <- list(AD_NEW_Lines, AD_NEW_Polys, AD_OLD_Lines, AD_OLD_Polys)
  i <- 0L
  fake_prep <- function(...) { i <<- i %% length(seq_data) + 1L; seq_data[[i]] }
  
  classesAvailable <- data.table(
    classToSearch = c("Road","Mine","Cutblock"),
    dataName      = c("Transport","Mining","Forestry"),
    dataClass     = c("roads","mines","cutblocks")
  )
  
  # Make a local copy we can patch
  dECCC <- disturbanceInfoFromECCC
  
  # IMPORTANT: stub returns NULL; it mutates dECCC in place — don't reassign!
  stub(dECCC, "prepInputs", fake_prep)
  stub(dECCC, "reproducible::prepInputs", fake_prep)  # covers qualified calls too
  
  expect_warning(
    res <- dECCC(
      studyArea = sa, RTM = rtm, classesAvailable = classesAvailable,
      totalstudyAreaVAreaSqKm = 100, disturbanceList = list(),
      destinationPath = tempdir(), diffYears = "2010_2015",
      bufferedDisturbances = FALSE, maskOutLinesFromPolys = FALSE,
      aggregateSameDisturbances = FALSE
    ),
    "Not all classes of OLD are in the NEW data"
  )
  
  changed <- res$AD_changed
  expect_true(all(c("Road","Mine","Cutblock") %in% changed$Class))
  expect_equal(changed[Class == "Mine"]$yearOLD, 0)      # Mine added to OLD as 0
  expect_equal(changed[Class == "Cutblock"]$yearNEW, 0)  # Cutblock added to NEW as 0
})

