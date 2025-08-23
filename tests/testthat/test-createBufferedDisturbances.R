# test-createBufferedDisturbances.R

testthat::test_that("setup: required packages are available", {
  testthat::skip_if_not_installed("terra")
  testthat::skip_if_not_installed("raster")
  testthat::skip_if_not_installed("sf")
  testthat::skip_if_not_installed("fasterize")
})

suppressPackageStartupMessages({
  library(testthat)
  library(terra)
  library(raster)
  library(sf)
})

# --- helpers ---------------------------------------------------------------

make_template <- function(n = 100L, res_m = 10, crs_str = "EPSG:3005") {
  # square raster, side = n * res_m meters
  r <- rast(nrows = n, ncols = n, xmin = 0, xmax = n * res_m, ymin = 0, ymax = n * res_m, vals = 1)
  crs(r) <- crs_str
  r
}

make_study_area <- function(r) {
  sa <- as.polygons(ext(r))
  crs(sa) <- crs(r)
  sa
}

pt_sv <- function(x, y, crs_str) {
  v <- vect(cbind(x, y), type = "points", crs = crs_str)
  v
}

# Create a small "non-potential" polygon feature
make_nonpot_poly <- function(r, x_frac = 0.3, y_frac = 0.3, inner_buf = 15) {
  crs_str <- crs(r)
  ext_vals <- as.vector(ext(r))
  x <- ext_vals[1] + x_frac * (ext_vals[2] - ext_vals[1])
  y <- ext_vals[3] + y_frac * (ext_vals[4] - ext_vals[3])
  p <- pt_sv(x, y, crs_str)
  # small circle
  terra::buffer(p, width = inner_buf)
}

# Create a small "potential" polygon feature (should be ignored)
make_potential_poly <- function(r, x_frac = 0.7, y_frac = 0.7, inner_buf = 15) {
  make_nonpot_poly(r, x_frac, y_frac, inner_buf)
}

# Create a binary 1/NA SpatRaster disturbance block
make_binary_raster_block <- function(r, x_idx = 20:30, y_idx = 20:30) {
  # clamp indices to raster dimensions to avoid subscript out-of-bounds
  x_idx <- unique(pmax(1, pmin(ncol(r), x_idx)))
  y_idx <- unique(pmax(1, pmin(nrow(r), y_idx)))
  
  r_block <- rast(r)
  values(r_block) <- NA
  
  # set requested row/col cells to 1 using cell indices (avoids transpose confusion)
  if (length(x_idx) > 0 && length(y_idx) > 0) {
    rc <- expand.grid(row = y_idx, col = x_idx)
    cells <- terra::cellFromRowCol(r_block, rc$row, rc$col)
    r_block[cells] <- 1
  }
  r_block[r_block != 1] <- NA
  r_block
}

# Build a disturbanceList with nested sectors/classes
make_disturbance_list <- function(r) {
  crs_str <- crs(r)
  # non-potential vector
  poly_np1 <- make_nonpot_poly(r, 0.3, 0.3, inner_buf = 12)
  poly_np2 <- make_nonpot_poly(r, 0.7, 0.3, inner_buf = 12)
  # potential vector (ignored)
  poly_pot <- make_potential_poly(r, 0.7, 0.7, inner_buf = 20)
  
  # non-potential raster (1-coded)
  rast_np <- make_binary_raster_block(r, 40:50, 40:50)
  # potential raster (ignored; indices may exceed dims for small templates -> clamped)
  rast_pot <- make_binary_raster_block(r, 60:65, 60:65)
  
  list(
    forestry = list(
      cutblocks = poly_np1,
      potentialCutblocks = poly_pot
    ),
    mining = list(
      mining = poly_np2
    ),
    roads = list(
      roads = rast_np,
      potentialRoads = rast_pot
    )
  )
}

# convenience: area in m^2 for SpatVector polygons
area_m2 <- function(v) {
  if (is.null(v) || nrow(v) == 0) return(0)
  as.numeric(terra::expanse(v, unit = "m"))
}

# --- tests -----------------------------------------------------------------

test_that("Vector mode: returns dissolved SpatVector and excludes potential layers", {
  skip_if_not_installed("fasterize") # used inside the function even in vector workflow
  r <- make_template(n = 100, res_m = 10)
  sa <- make_study_area(r)
  
  dl_with_pot <- make_disturbance_list(r)
  dl_no_pot <- dl_with_pot
  # Remove potential entries explicitly to create a reference output
  dl_no_pot$forestry$potentialCutblocks <- NULL
  dl_no_pot$roads$potentialRoads <- NULL
  
  # Use a modest external buffer
  buf <- 40
  
  out_with <- createBufferedDisturbances(
    disturbanceList = dl_with_pot,
    bufferSize = buf,
    rasterToMatch = r,
    studyArea = sa,
    currentTime = 1L,
    convertToRaster = FALSE
  )
  
  out_no   <- createBufferedDisturbances(
    disturbanceList = dl_no_pot,
    bufferSize = buf,
    rasterToMatch = r,
    studyArea = sa,
    currentTime = 1L,
    convertToRaster = FALSE
  )
  
  # Class & structure
  expect_s4_class(out_with, "SpatVector")
  expect_equal(nrow(out_with), 1L) # dissolved into single (multi)feature
  
  # Exclusion of potential layers -> identical geometry to manual "no potential" list
  expect_true(all.equal(area_m2(out_with), area_m2(out_no), tolerance = 1e-6) == TRUE)
})

test_that("Vector mode: empty raster inputs are ignored (NULL dropped) and union still works", {
  skip_if_not_installed("fasterize")
  
  r <- make_template(n = 60, res_m = 10)
  sa <- make_study_area(r)
  dl <- make_disturbance_list(r)
  
  # Add an empty raster layer (all NA or all 0 -> will be treated as empty)
  empty_r <- rast(r); values(empty_r) <- 0
  dl$seismic <- list(seismic = empty_r)
  
  out <- suppressWarnings(createBufferedDisturbances(
    disturbanceList = dl,
    bufferSize = 30,
    rasterToMatch = r,
    studyArea = sa,
    currentTime = 1L,
    convertToRaster = FALSE
  ))
  
  expect_s4_class(out, "SpatVector")
  expect_true(area_m2(out) > 0)
})

test_that("Raster mode (SpatRaster template): 1/0/NA semantics and class parity", {
  skip_if_not_installed("fasterize")
  
  r <- make_template(n = 80, res_m = 10)
  sa <- make_study_area(r)
  dl <- make_disturbance_list(r)
  
  out <- createBufferedDisturbances(
    disturbanceList = dl,
    bufferSize = 25,
    rasterToMatch = r,          # SpatRaster template
    studyArea = sa,
    currentTime = 1L,
    convertToRaster = TRUE
  )
  
  expect_s4_class(out, "SpatRaster")
  
  vals <- unique(values(out))
  vals <- vals[!is.na(vals)]
  expect_true(all(vals %in% c(0, 1)))
  expect_true(any(values(out) == 1, na.rm = TRUE)) # some disturbance exists
})

test_that("Raster mode (RasterLayer template): class tracks template class", {
  skip_if_not_installed("fasterize")
  
  r <- make_template(n = 80, res_m = 10)
  sa <- make_study_area(r)
  dl <- make_disturbance_list(r)
  
  r_template <- raster::raster(r) # RasterLayer template
  
  out <- createBufferedDisturbances(
    disturbanceList = dl,
    bufferSize = 25,
    rasterToMatch = r_template,
    studyArea = sa,
    currentTime = 1L,
    convertToRaster = TRUE
  )
  
  expect_true(inherits(out, "RasterLayer"))  # class parity with template
  # still 1/0/NA
  vals <- raster::unique(out)
  vals <- vals[!is.na(vals)]
  expect_true(all(vals %in% c(0, 1)))
})

test_that("Raster mode: potential layers do not change the result", {
  skip_if_not_installed("fasterize")
  
  r <- make_template(n = 60, res_m = 10)
  sa <- make_study_area(r)
  
  dl_with_pot <- make_disturbance_list(r)
  dl_no_pot <- dl_with_pot
  dl_no_pot$forestry$potentialCutblocks <- NULL
  dl_no_pot$roads$potentialRoads <- NULL
  
  out_with <- suppressWarnings(createBufferedDisturbances(
    disturbanceList = dl_with_pot,
    bufferSize = 20,
    rasterToMatch = r,
    studyArea = sa,
    currentTime = 1L,
    convertToRaster = TRUE
  ))
  
  out_no <- suppressWarnings(createBufferedDisturbances(
    disturbanceList = dl_no_pot,
    bufferSize = 20,
    rasterToMatch = r,
    studyArea = sa,
    currentTime = 1L,
    convertToRaster = TRUE
  ))
  
  # Compare number of disturbed (==1) cells; identical when potentials are ignored
  n1_with <- sum(values(out_with) == 1, na.rm = TRUE)
  n1_no   <- sum(values(out_no) == 1,   na.rm = TRUE)
  expect_equal(n1_with, n1_no)
})

test_that("Inputs: raster disturbance layers must be 1-coded to be kept", {
  skip_if_not_installed("fasterize")
  
  r <- make_template(n = 40, res_m = 10)
  sa <- make_study_area(r)
  
  # non-potential 1-coded block
  r1 <- make_binary_raster_block(r, 10:15, 10:15)
  
  # non-potential but NOT 1-coded -> should be dropped by the function
  r2 <- rast(r); values(r2) <- NA
  m2 <- matrix(0, nrow = nrow(r), ncol = ncol(r))
  m2[20:25, 20:25] <- 2L           # code 2 (not 1)
  values(r2) <- as.vector(t(m2))
  
  dl <- list(
    mines = list(mining = r1, mining2 = r2) # neither name contains "potential"; only r1 should survive
  )
  
  out <- createBufferedDisturbances(
    disturbanceList = dl,
    bufferSize = 15,
    rasterToMatch = r,
    studyArea = sa,
    currentTime = 1L,
    convertToRaster = TRUE
  )
  
  # We expect only the 1-coded block to contribute
  n1 <- sum(values(out) == 1, na.rm = TRUE)
  expect_true(n1 > 0)
  
  # crude sanity: removing r1 should clear all 1s (empty -> internal empty polygons emit warnings)
  dl_only_bad <- list(mines = list(mining2 = r2))
  out_bad <- suppressWarnings(createBufferedDisturbances(
    disturbanceList = dl_only_bad, bufferSize = 15,
    rasterToMatch = r, studyArea = sa, currentTime = 1L,
    convertToRaster = TRUE
  ))
  n1_bad <- sum(values(out_bad) == 1, na.rm = TRUE)
  expect_equal(n1_bad, 0)
})
