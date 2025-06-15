library(testthat)
library(sf)
library(terra)

# Create a simple horizontal line from (0,0) to (10,0)
line_sf <- st_as_sf(st_sfc(st_linestring(rbind(c(0,0), c(10,0)))), crs = 4326)
# Center point directly above the start of the line
center_pt <- vect(matrix(c(0, 10), ncol = 2), type = "points", crs = crs(vect(line_sf)))

# Common parameters for tests
distanceBetween <- 5
fraction <- 0.0
linesLengthExpr <- "perim(Line)"  # keep original length

# 1. Input validation

test_that("invalid Line input throws error", {
  expect_error(
    seismicGrids(Line = "not a line", 
                 centerPoint = center_pt,
                 originalLayer = NULL,
                 linesLength = linesLengthExpr,
                 distanceBetweenLines = distanceBetween,
                 howMany = c(1,1),
                 existingLine = TRUE),
    "Line musst be sf or SpatVector"
  )
})

test_that("invalid centerPoint input throws error", {
  expect_error(
    seismicGrids(Line = line_sf, 
                 centerPoint = "not a point", 
                 originalLayer = NULL,
                 linesLength = linesLengthExpr,
                 distanceBetweenLines = distanceBetween,
                 howMany = c(1,1),
                 existingLine = TRUE),
    "centerPoint musst be a SpatVector"
  )
})

# 2. Basic functionality: correct number of lines

test_that("returns expected number of parallel and crossing lines", {
  # We ask for 2 parallels (excluding original) and 2 crossings
  out <- seismicGrids(
    Line                 = line_sf,
    centerPoint          = center_pt,
    originalLayer        = NULL,
    linesLength          = linesLengthExpr,
    distanceBetweenLines = distanceBetween,
    howMany              = c(2, 2),
    existingLine         = FALSE
  )
  # Expect a SpatVector result
  expect_true(is(out, "SpatVector"))
  # Expect total lines = parallels + crossings (no random here)
  expect_equal(nrow(out), 4)
})

# 3. Geometry: spacing between parallels matches distanceBetweenLines

test_that("parallel lines are spaced correctly", {
  out <- seismicGrids(
    Line                 = line_sf,
    centerPoint          = center_pt,
    originalLayer        = NULL,
    linesLength          = linesLengthExpr,
    distanceBetweenLines = distanceBetween,
    howMany              = c(1, 0),  # only one parallel, no crossings
    existingLine         = TRUE
  )
  # Extract parallel line and original: distances should be exactly 5 units
  original_coords <- st_coordinates(line_sf)[1,]
  parallel_coords <- st_coordinates(st_as_sf(out)[1,])[1,]
  dist <- sqrt((parallel_coords["X"] - original_coords["X"])^2 + 
                 (parallel_coords["Y"] - original_coords["Y"])^2)
  expect_equal(dist, distanceBetween)
})

# 4. Perpendicular orientation

test_that("crossing lines are approximately perpendicular to base line", {
  out <- seismicGrids(
    Line                 = line_sf,
    centerPoint          = center_pt,
    originalLayer        = NULL,
    linesLength          = linesLengthExpr,
    distanceBetweenLines = distanceBetween,
    howMany              = c(1, 1),
    existingLine         = FALSE
  )
  # Extract the crossing line
  crossing <- st_as_sf(out)[nrow(out),]
  coords <- st_coordinates(crossing)
  # Compute angle of segment
  delta <- coords[2,] - coords[1,]
  angle <- atan2(delta[2], delta[1])
  # Base line is horizontal, so crossing should be vertical (~pi/2)
  expect_true(abs(abs(angle) - pi/2) < 1e-6)
})
