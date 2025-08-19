library(testthat)
library(terra)
library(sf)
library(truncnorm)

# Helper: create simple lines as a SpatVector via WKT
make_lines <- function(wkt_vec, cluster_id, crs_def = 
                         "+proj=aeqd +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs") {
  if (length(wkt_vec) == 0) {
    empty <- vect()  # Create empty SpatVector first
    crs(empty) <- crs_def  # Then assign CRS
    empty$Pot_Clus <- character(0)
    return(empty)
  }
  sv <- vect(wkt_vec, crs = crs_def)
  sv$Pot_Clus <- cluster_id
  return(sv)
}

normalize_angle <- function(a) {
  a <- a %% 180
  if(a > 90) a - 180 else a
}


# Test that simulateLines returns a SpatVector with the correct number of lines
test_that("simulateLines returns correct number of lines", {
  wkt <- rep("LINESTRING (0 0, 1 0)", 3)
  lines <- do.call(rbind, lapply(seq_along(wkt), function(i) make_lines(wkt[i], 1)))
  
  result <- simulateLines(lines, distThreshold = 10, distNewLinesFact = 1, refinedStructure = FALSE)
  
  expect_s4_class(result, "SpatVector")
  expect_equal(nrow(result), 3)
})

# Test that a single-line cluster preserves the exact length
test_that("single-line cluster preserves length", {
  wkt <- "LINESTRING (0 0, 0 2)"
  single <- make_lines(wkt, 42)
  out <- simulateLines(single, distThreshold = 5, distNewLinesFact = 1, refinedStructure = FALSE)
  
  expect_equal(out$lineLength, terra::perim(single), tolerance = 1e-6)
})

# Test that randomized structure (refinedStructure=FALSE) draws lengths within bounds
test_that("random structure samples lengths within original range", {
  wkts <- c("LINESTRING (0 0, 3 0)",
            "LINESTRING (0 0, 4 0)",
            "LINESTRING (0 0, 5 0)",
            "LINESTRING (0 0, 6 0)",
            "LINESTRING (0 0, 7 0)")
  lines <- do.call(rbind, lapply(seq_along(wkts), function(i) make_lines(wkts[i], 100)))
  
  out <- simulateLines(lines, distThreshold = 20, distNewLinesFact = 2, refinedStructure = FALSE)
  
  sims <- out$lineLength
  eps  <- 1e-6
  L    <- terra::perim(lines)
  
  expect_true(all(sims >= min(L) - eps))
  expect_true(all(sims <= max(L) + eps))
})

# Test refinedStructure generates perpendicular pairs for a 2-line perpendicular cluster
test_that("refinedStructure preserves perpendicular orientation for perpendicular inputs", {
  wkt_x <- "LINESTRING (5 5, 105 5)"  # Horizontal
  wkt_y <- "LINESTRING (10 10, 10 100)"  # Vertical
  ex <- make_lines(wkt_x, 5)
  ey <- make_lines(wkt_y, 5)
  lines <- rbind(ex, ey)
  
  set.seed(1123)
  out <- simulateLines(lines, distThreshold = 10, distNewLinesFact = 1, refinedStructure = TRUE)
  sim_angles <- sapply(1:nrow(out), function(i) calculateLineAngle(out[i, ]))
  ang_diff <- abs(sim_angles[1] - sim_angles[2])
  ang_diff <- min(ang_diff, 180 - ang_diff)
  expect_true(all(abs(ang_diff - 90) <= 5))
})

# Test that no NA lengths are produced in multi-line clusters
test_that("no NA lengths for multi-line clusters", {
  wkts <- rep("LINESTRING (0 0, 2 2)", 4)
  lines <- do.call(rbind, lapply(seq_along(wkts), function(i) make_lines(wkts[i], 7)))
  
  out <- simulateLines(lines, distThreshold = 15, distNewLinesFact = 1, refinedStructure = TRUE)
  
  expect_false(any(is.na(out$lineLength)))
})

# Test behavior with different CRS inputs
# 1. Geographic CRS: should run without error, preserve count
# 2. Missing CRS: should warn about unknown CRS and still produce output
# 3. Projected CRS (e.g. UTM): should run silently

test_that("simulateLines works with geographic CRS", {
  wkts <- rep("LINESTRING (0 0, 1 0)", 3)
  geo_lines <- do.call(rbind, lapply(seq_along(wkts), function(i) make_lines(wkts[i], 1, "+proj=longlat +datum=WGS84")))
  expect_silent({
    res <- simulateLines(geo_lines, distThreshold = 10, distNewLinesFact = 1, refinedStructure = FALSE)
  })
  expect_equal(nrow(res), 3)
})

test_that("simulateLines warns when CRS is missing", {
  wkt <- "LINESTRING (0 0, 0 1)"
  no_crs <- vect(wkt)
  no_crs$Pot_Clus <- 9
  expect_warning(res <- simulateLines(no_crs, distThreshold = 5, distNewLinesFact = 1, refinedStructure = FALSE), "unknown CRS")
  expect_s4_class(res, "SpatVector")
  expect_equal(nrow(res), 1)
})

test_that("simulateLines silent with projected CRS", {
  wkt <- rep("LINESTRING (0 0, 2 0)", 2)
  utm_lines <- do.call(rbind, lapply(seq_along(wkt), function(i) make_lines(wkt[i], 3, "EPSG:32633")))
  expect_silent({
    out <- simulateLines(utm_lines, distThreshold = 8, distNewLinesFact = 1, refinedStructure = FALSE)
  })
  expect_equal(nrow(out), 2)
})

# Edge-case tests
test_that("constant-length clusters yield identical simulated lengths", {
  wkt <- rep("LINESTRING (0 0, 10 0)", 4)
  lines <- do.call(rbind, lapply(seq_along(wkt), function(i) make_lines(wkt[i], 1)))
  orig_len <- terra::perim(lines)[1]
  set.seed(42)
  out <- simulateLines(lines, distThreshold = 10, distNewLinesFact = 1, refinedStructure = FALSE)
  
  expect_equal(out$calculatedLength, rep(orig_len, 4), tolerance = 1e-6)
})

test_that("errors on empty input", {
  # Create an empty spatial object with a specified CRS
  empty <- st_sf(geometry = st_sfc(), crs = "+proj=aeqd +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs")
  empty$Pot_Clus <- character(0)
  expect_error(simulateLines(empty), "no input features")
})

test_that("simulateLines preserves row count for odd-number clusters with refinedStructure", {
  wkts <- c("LINESTRING (0 0, 1 0)",
            "LINESTRING (0 0, 0 1)",
            "LINESTRING (0 0, 2 0)")
  lines <- do.call(rbind, lapply(seq_along(wkts), function(i) make_lines(wkts[i], 2)))
  set.seed(123)
  out <- simulateLines(lines, distThreshold = 5, distNewLinesFact = 1, refinedStructure = TRUE)
  expect_equal(nrow(out), 3)
})

test_that("simulateLines errors when Pot_Clus column is missing", {
  wkt <- "LINESTRING (0 0, 1 0)"
  no_pc <- vect(wkt, crs = "+proj=aeqd +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs")
  expect_error(simulateLines(no_pc), "Pot_Clus")
})

# Test for 2-line parallel cluster
test_that("refinedStructure preserves parallel orientation for 2-line parallel inputs", {
  wkt1 <- "LINESTRING (0 100, 100 200)"  # ~45°
  wkt2 <- "LINESTRING (0 0, 200 200)"  # ~45°
  ex <- make_lines(wkt1, 100)
  ey <- make_lines(wkt2, 100)
  lines <- rbind(ex, ey)
  
  set.seed(101)
  out <- simulateLines(lines, distThreshold = 5, distNewLinesFact = 1, refinedStructure = TRUE)
  sim_angles <- sapply(seq_len(nrow(out)), function(i) calculateLineAngle(out[i, ]))
  ang_diff <- abs(normalize_angle(sim_angles[1]) - normalize_angle(sim_angles[2]))
  ang_diff <- min(ang_diff, 180 - ang_diff)
  expect_true(ang_diff <= 5)
})

# Test for multi-line refined structure >2: multiple parallel pairs
test_that("refinedStructure correctly simulates multiple parallel pairs in multi-line cluster", {
  wkts <- c(
    "LINESTRING (0 0, 1 0)",
    "LINESTRING (0 0, 1 0)",
    "LINESTRING (0 0, 0 1)",
    "LINESTRING (0 0, 0 1)"
  )
  lines <- do.call(rbind, lapply(seq_along(wkts), function(i) make_lines(wkts[i], 200)))
  
  set.seed(2025)
  out <- simulateLines(lines, distThreshold = 10, distNewLinesFact = 1, refinedStructure = TRUE)
  sim_angles <- sapply(seq_len(nrow(out)), function(i) calculateLineAngle(out[i, ]))
  # Expect first two to be parallel, next two also parallel
  expect_true(abs(sim_angles[1] - sim_angles[2]) <= 5)
  expect_true(abs(sim_angles[3] - sim_angles[4]) <= 5)
})

# 1. Spatial extent within buffered centroid envelope
test_that("simulated lines lie within expected buffered extent", {
  wkts <- rep("LINESTRING (0 0, 0 1)", 2)
  lines <- do.call(rbind, lapply(seq_along(wkts), function(i) make_lines(wkts[i], 3)))
  distThreshold <- 10; fact <- 2; buffDist <- distThreshold * fact
  
  set.seed(222)
  out <- simulateLines(lines, distThreshold, fact, FALSE)
  
  ext_orig <- terra::ext(lines)
  ext_expected <- c(
    ext_orig[1] - buffDist,
    ext_orig[2] + buffDist,
    ext_orig[3] - buffDist,
    ext_orig[4] + buffDist
  )
  ext_out <- terra::ext(out)
  
  expect_true(ext_out[1] >= ext_expected[1] - 1e-6)
  expect_true(ext_out[2] <= ext_expected[2] + 1e-6)
  expect_true(ext_out[3] >= ext_expected[3] - 1e-6)
  expect_true(ext_out[4] <= ext_expected[4] + 1e-6)
})

# 2. Leftover lines have non-paired orientations
test_that("leftover lines have non-paired orientations", {
  wkts <- c(
    "LINESTRING (0 0, 1 0)",
    "LINESTRING (0 0, 1 0)",
    "LINESTRING (0 0, 0 1)",
    "LINESTRING (0 0, 0 1)",
    "LINESTRING (0 0, 2 0)"
  )
  lines <- do.call(rbind, lapply(seq_along(wkts), function(i) make_lines(wkts[i], 300)))
  
  set.seed(789)
  out <- simulateLines(lines, 5, 1, TRUE)
  sim_angles <- sapply(seq_len(nrow(out)), function(i) calculateLineAngle(out[i, ]))
  tol <- 5
  is_parallel <- function(a, b) abs(a - b) <= tol || abs(abs(a - b) - 180) <= tol
  is_perp_to_12 <- abs(abs(sim_angles[5] - sim_angles[1]) - 90) <= tol
  is_perp_to_34 <- abs(abs(sim_angles[5] - sim_angles[3]) - 90) <= tol
  
  expect_false(
    is_parallel(sim_angles[5], sim_angles[1]) ||
      is_parallel(sim_angles[5], sim_angles[3]) ||
      is_perp_to_12 || is_perp_to_34
  )
})

# 3. Mixed parallel and perpendicular patterns
test_that("simulateLines maintains complex spatial relationships", {
  wkts <- c(
    "LINESTRING (0 0, 1 0)",  # parallel group
    "LINESTRING (0 0, 1 0)",
    "LINESTRING (0 0, 0 1)",  # perpendicular group
    "LINESTRING (0 0, 1 1)",  # 45°
    "LINESTRING (0 0, -1 1)"  # other
  )
  lines <- do.call(rbind, lapply(seq_along(wkts), function(i) make_lines(wkts[i], 1)))
  set.seed(123)
  out <- simulateLines(lines, distThreshold = 10, distNewLinesFact = 1, refinedStructure = TRUE)
  angles <- sapply(seq_len(nrow(out)), function(i) calculateLineAngle(out[i, ]))
  tol <- 5
  # First two remain parallel
  expect_true(abs(angles[1] - angles[2]) <= tol)
})

# 4. Angle tolerance near-zero degrees treated as parallel
test_that("simulateLines correctly handles angle tolerance", {
  wkt1 <- "LINESTRING (0 0, 1 0)"
  wkt2 <- "LINESTRING (0 0, 1 0.087)"  # ~5° difference
  ex <- make_lines(wkt1, 1)
  ey <- make_lines(wkt2, 1)
  lines <- rbind(ex, ey)
  set.seed(123)
  out <- simulateLines(lines, distThreshold = 10, distNewLinesFact = 1, refinedStructure = TRUE)
  angles <- sapply(seq_len(nrow(out)), function(i) calculateLineAngle(out[i, ]))
  expect_true(abs(angles[1] - angles[2]) <= 5)
})
