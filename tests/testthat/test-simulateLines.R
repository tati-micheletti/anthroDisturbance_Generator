library(testthat)
library(terra)
library(truncnorm)

# Helper: create simple lines as a SpatVector via WKT
make_lines <- function(wkt_vec, cluster_id, crs_def = "+proj=aeqd +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs") {
  sv <- vect(wkt_vec, crs = crs_def)
  sv$Pot_Clus <- cluster_id
  return(sv)
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
  expect_equal(out$lineLength, lengths(single), tolerance = 1e-6)
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
  eps <- 1e-6
  expect_true(all(sims >= min(lengths(lines)) - eps))
  expect_true(all(sims <= max(lengths(lines)) + eps))
})

# Test refinedStructure generates perpendicular pairs for a 2-line perpendicular cluster
test_that("refinedStructure preserves perpendicular orientation for perpendicular inputs", {
  wkt_x <- "LINESTRING (0 0, 1 0)"
  wkt_y <- "LINESTRING (0 0, 0 1)"
  ex <- make_lines(wkt_x, 5)
  ey <- make_lines(wkt_y, 5)
  lines <- rbind(ex, ey)
  
  set.seed(123)
  out <- simulateLines(lines, distThreshold = 10, distNewLinesFact = 1, refinedStructure = TRUE)
  
  sim_angles <- sapply(1:nrow(out), function(i) calculateLineAngle(out[i, ]))
  ang_diff <- abs(diff(sim_angles))
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
  wkt <- rep("LINESTRING (0 0, 1 0)", 4)
  lines <- do.call(rbind, lapply(seq_along(wkt), function(i) make_lines(wkt[i], 1)))
  orig_len <- lengths(lines)[1]
  set.seed(42)
  out <- simulateLines(lines, distThreshold = 10, distNewLinesFact = 1, refinedStructure = FALSE)
  expect_equal(out$lineLength, rep(orig_len, 4), tolerance = 1e-6)
})

test_that("simulateLines errors on completely empty input", {
  # Create an empty SpatVector with proper CRS
  empty <- make_lines(character(0), cluster_id = integer(0),
                      crs_def = "+proj=aeqd +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs")
  expect_error(simulateLines(empty), "subscript out of bounds")
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
  wkt1 <- "LINESTRING (0 0, 1 1)"  # ~45°
  wkt2 <- "LINESTRING (0 0, 2 2)"  # ~45°
  ex <- make_lines(wkt1, 100)
  ey <- make_lines(wkt2, 100)
  lines <- rbind(ex, ey)
  
  set.seed(2025)
  out <- simulateLines(lines, distThreshold = 5, distNewLinesFact = 1, refinedStructure = TRUE)
  sim_angles <- sapply(seq_len(nrow(out)), function(i) calculateLineAngle(out[i, ]))
  expect_true(abs(diff(sim_angles)) <= 5)
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
