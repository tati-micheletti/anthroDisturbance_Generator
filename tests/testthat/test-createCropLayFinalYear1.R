library(testthat)
library(terra)
library(reproducible)
library(data.table)
library(tictoc)
library(dplyr)

# Disable reproducible caching for speed
options(reproducible.useCache = FALSE)
assignInNamespace("Cache",
                  function(fun, ...) fun(...),
                  ns = "reproducible"
)

# Helper to disable caching between tests
reset_cache <- function() {
  options(reproducible.cachePath = tempfile())
}

# Helper to unpack function result
unpack_res <- function(res) {
  if (inherits(res, "list")) return(res$lines)
  res
}

# 1. Cropping
test_that("Lines outside the potential polygon are dropped (cropping)", {
  reset_cache()
  pot <- vect("POLYGON ((0 0, 0 1000, 1000 1000, 1000 0, 0 0))", crs = "EPSG:32633")
  pot$Potential <- 1L
  line_inside  <- vect("LINESTRING (100 100, 500 500)", crs = "EPSG:32633")
  line_outside <- vect("LINESTRING (2000 2000, 3000 3000)", crs = "EPSG:32633")
  line_crossing <-vect("LINESTRING (750 750, 2000 750)", crs = "EPSG:32633")
  Lay <- rbind(line_inside, line_outside, line_crossing)
  
  plot(Lay)
  plot(pot, add=TRUE)
  
  suppressMessages({
    res <- createCropLayFinalYear1(
      Lay                    = Lay,
      potLayTopValid         = pot,
      runClusteringInParallel= FALSE,
      clusterDistance        = 5,
      studyAreaHash          = "test_crop"
    )
  })
  lines <- unpack_res(res)
  
  plot(lines, add=TRUE, col="red")
  
  expect_s4_class(lines, "SpatVector")
  expect_equal(nrow(lines), 2)
  rel <- terra::relate(lines, pot, relation = "T********", pairs = FALSE)
  expect_true(any(rel))
})

# 2. Deduplication
test_that("Overlapping collinear lines get deduplicated (keep the longer one)", {
  reset_cache()
  pot <- vect("POLYGON ((0 0, 0 1000, 1000 1000, 1000 0, 0 0))", crs = "EPSG:32631")
  pot$Potential <- 1L
  line1 <- vect("LINESTRING (1 1, 8 1)", crs = "EPSG:32631")
  line2 <- vect("LINESTRING (1 1, 5 1)", crs = "EPSG:32631")
  Lay <- rbind(line1, line2)
  
  suppressMessages({
    res <- createCropLayFinalYear1(
      Lay                     = Lay,
      potLayTopValid          = pot,
      runClusteringInParallel = FALSE,
      clusterDistance         = 5,
      studyAreaHash           = "test_dedupe"
    )
  })
  lines <- unpack_res(res)
  
  expect_s4_class(lines, "SpatVector")
  expect_equal(nrow(lines), 1)
  total_length <- terra::perim(lines)
  expect_equal(total_length, 7, tolerance = 1e-6)
})

# 3. Angle exemption
test_that("Intersecting lines at different angles are not removed", {
  reset_cache()
  pot <- vect("POLYGON ((0 0, 0 1000, 1000 1000, 1000 0, 0 0))", crs = "EPSG:32633")
  pot$Potential <- 1L
  horizontal <- vect("LINESTRING (2 5, 8 5)", crs = "EPSG:32633")
  vertical   <- vect("LINESTRING (5 2, 5 8)", crs = "EPSG:32633")
  Lay <- rbind(horizontal, vertical)
  
  suppressMessages({
    res <- createCropLayFinalYear1(
      Lay                     = Lay,
      potLayTopValid          = pot,
      runClusteringInParallel = FALSE,
      clusterDistance         = 5,
      studyAreaHash           = "test_angles"
    )
  })
  lines <- unpack_res(res)
  
  expect_s4_class(lines, "SpatVector")
  expect_equal(nrow(lines), 2)
  angles <- vapply(seq_len(nrow(lines)), function(i) calculateLineAngle(lines[i, ]), numeric(1))
  expect_true(any(abs(angles) < 1))
  expect_true(any(abs(abs(angles) - 90) < 1))
})

# 4. Clustering
test_that("Clustering assigns separate clusters to distant line groups", {
  reset_cache()
  pot <- vect("POLYGON ((0 0, 0 100, 100 100, 100 0, 0 0))", crs = "EPSG:32633")
  pot$Potential <- 1L
  near_line1 <- vect("LINESTRING (10 10, 20 10)", crs = "EPSG:32633")
  near_line2 <- vect("LINESTRING (15 15, 25 15)", crs = "EPSG:32633")
  far_line   <- vect("LINESTRING (80 80, 90 80)", crs = "EPSG:32633")
  Lay <- rbind(near_line1, near_line2, far_line)
  
  suppressMessages({
    res <- createCropLayFinalYear1(
      Lay                     = Lay,
      potLayTopValid          = pot,
      runClusteringInParallel = FALSE,
      clusterDistance         = 10,
      studyAreaHash           = "test_cluster"
    )
  })
  lines <- unpack_res(res)
  
  expect_s4_class(lines, "SpatVector")
  expect_equal(nrow(lines), 3)
  clusters <- unique(lines$cluster)
  expect_equal(length(clusters), 2)
  sizes <- as.integer(table(lines$cluster))
  expect_true(any(sizes == 2))
  expect_true(any(sizes == 1))
})

# Edge cases

# Empty input
test_that("Empty Lay yields empty result", {
  reset_cache()
  tmp <- vect("LINESTRING (0 0, 0 0)", crs = "EPSG:32633")
  Lay <- tmp[0, ]
  pot <- vect("POLYGON ((0 0, 0 1000, 1000 1000, 1000 0, 0 0))", crs = "EPSG:32633")
  pot$Potential <- 1L
  
  suppressWarnings({
    res <- createCropLayFinalYear1(
      Lay                     = Lay,
      potLayTopValid          = pot,
      runClusteringInParallel = FALSE,
      clusterDistance         = 5,
      studyAreaHash           = "edge_empty"
    )
  })
  lines <- unpack_res(res)
  expect_true(nrow(lines) == 0)
})

# All outside
test_that("All lines outside potential yield empty result", {
  reset_cache()
  line_out <- vect("LINESTRING (2000 2000, 3000 3000)", crs = "EPSG:32633")
  Lay <- rbind(line_out)
  pot <- vect("POLYGON ((0 0, 0 1000, 1000 1000, 1000 0, 0 0))", crs = "EPSG:32633")
  pot$Potential <- 1L
  
  suppressWarnings({
    res <- createCropLayFinalYear1(
      Lay                     = Lay,
      potLayTopValid          = pot,
      runClusteringInParallel = FALSE,
      clusterDistance         = 5,
      studyAreaHash           = "edge_all_out"
    )
  })
  lines <- unpack_res(res)
  expect_true(nrow(lines) == 0)
})

# Zero distance clustering
test_that("Zero clusterDistance still clusters into singletons", {
  reset_cache()
  pot <- vect("POLYGON ((0 0, 0 1000, 1000 1000, 1000 0, 0 0))", crs = "EPSG:32633")
  pot$Potential <- 1L
  line1 <- vect("LINESTRING (1 1, 2 2)", crs = "EPSG:32633")
  line2 <- vect("LINESTRING (3 3, 4 4)", crs = "EPSG:32633")
  Lay <- rbind(line1, line2)
  
  suppressMessages({
    res <- createCropLayFinalYear1(
      Lay, pot,
      runClusteringInParallel = FALSE,
      clusterDistance         = 0,
      studyAreaHash           = "test_dist0"
    )
  })
  lines <- unpack_res(res)
  
  expect_s4_class(lines, "SpatVector")
  expect_equal(nrow(lines), 2)
  expect_true("cluster" %in% names(lines))
  clusters <- unique(lines$cluster)
  expect_equal(length(clusters), 2)
})

# Opposite-direction lines
test_that("Collinear lines in opposite directions are kept", {
  reset_cache()
  pot <- vect("POLYGON ((0 0, 0 1000, 1000 1000, 1000 0, 0 0))", crs = "EPSG:32631")
  pot$Potential <- 1L
  line1 <- vect("LINESTRING (1 1, 5 1)", crs = "EPSG:32631")  # → 0°
  line2 <- vect("LINESTRING (5 1, 1 1)", crs = "EPSG:32631")  # ← 180°
  Lay <- rbind(line1, line2)
  
  suppressMessages({
    res <- createCropLayFinalYear1(
      Lay, pot,
      runClusteringInParallel = FALSE,
      clusterDistance         = 5,
      studyAreaHash           = "test_opposite_angles"
    )
  })
  lines <- unpack_res(res)
  
  expect_s4_class(lines, "SpatVector")
  expect_equal(nrow(lines), 2)
})

# Errors on bad clusterDistance
test_that("Negative clusterDistance throws an error", {
  reset_cache()
  pot <- vect("POLYGON ((0 0, 0 1000, 1000 1000, 1000 0, 0 0))", crs = "EPSG:32633")
  pot$Potential <- 1L
  Lay <- vect("LINESTRING (1 1, 2 2)", crs = "EPSG:32633")
  
  expect_error(
    createCropLayFinalYear1(
      Lay,
      pot,
      runClusteringInParallel = FALSE,
      clusterDistance        = -10,
      studyAreaHash          = "test_neg_dist"
    ),
    "`clusterDistance` must be non\u2010negative"
  )
})

test_that("Non-numeric clusterDistance throws an error", {
  reset_cache()
  pot <- vect("POLYGON ((0 0, 0 1000, 1000 1000, 1000 0, 0 0))", crs = "EPSG:32633")
  pot$Potential <- 1L
  Lay <- vect("LINESTRING (1 1, 2 2)", crs = "EPSG:32633")
  
  expect_error(
    createCropLayFinalYear1(
      Lay,
      pot,
      runClusteringInParallel = FALSE,
      clusterDistance        = "foo",
      studyAreaHash          = "test_non_numeric_dist"
    ),
    "clusterDistance.*numeric"
  )
})

# Erase carve-out test
test_that("Erase carves out a 50m buffer from the potential polygon", {
  reset_cache()
  pot <- vect("POLYGON ((0 0, 0 1000, 1000 1000, 1000 0, 0 0))", crs = "EPSG:32633")
  pot$Potential <- 1L
  line <- vect("LINESTRING (150 150, 60 50)", crs = "EPSG:32633")
  Lay <- rbind(line)
  
  suppressMessages({
    res <- createCropLayFinalYear1(
      Lay,
      pot,
      runClusteringInParallel = FALSE,
      clusterDistance         = 5,
      studyAreaHash           = "test_erase"
    )
  })
  area <- res$availableArea
  
  buf <- buffer(line, width = 50)
  expect_false(any(terra::relate(area, buf, relation = "T********", pairs = FALSE)))
})
