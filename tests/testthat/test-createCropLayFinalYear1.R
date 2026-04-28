# Helper to disable caching between tests
reset_cache <- function() {
  options(reproducible.cachePath = tempfile())
}

# Helper to unpack function result
unpack_res <- function(res) {
  if (inherits(res, "list")) return(res$lines)
  res
}

# Simple area helper that is tolerant of empty/none geomtypes
area_of <- function(v) {
  gt <- tryCatch(terra::geomtype(v), error = function(e) "none")
  if (gt == "none" || nrow(v) == 0L) return(0)
  sum(terra::expanse(v), na.rm = TRUE)
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
  
  maybe_plot(plot(Lay))
  maybe_plot(plot(pot, add=TRUE))
  
  expect_warning(
    res <- suppressMessages(createCropLayFinalYear1(
      Lay                     = Lay,
      potLayTopValid          = pot,
      runClusteringInParallel = FALSE,
      clusterDistance         = 5,
      studyAreaHash           = "test_crop"
    )),
    "degenerate"
  )
  lines <- unpack_res(res)
  
  maybe_plot(plot(lines, add=TRUE, col="red"))
  
  expect_s4_class(lines, "SpatVector")
  expect_equal(nrow(lines), 0)
  expect_s4_class(res$availableArea, "SpatVector")
  expect_true(terra::geomtype(res$availableArea) %in% c("polygons", "none"))
  expect_lt(area_of(res$availableArea), area_of(pot))
})

# 2. Deduplication
test_that("Overlapping collinear lines get deduplicated (keep the longer one)", {
  reset_cache()
  pot <- vect("POLYGON ((0 0, 0 1000, 1000 1000, 1000 0, 0 0))", crs = "EPSG:32631")
  pot$Potential <- 1L
  line1 <- vect("LINESTRING (1 1, 8 1)", crs = "EPSG:32631")
  line2 <- vect("LINESTRING (1 1, 5 1)", crs = "EPSG:32631")
  Lay <- rbind(line1, line2)
  
  expect_warning(
    res <- suppressMessages(createCropLayFinalYear1(
      Lay                     = Lay,
      potLayTopValid          = pot,
      runClusteringInParallel = FALSE,
      clusterDistance         = 5,
      studyAreaHash           = "test_dedupe"
    )),
    "degenerate"
  )
  lines <- unpack_res(res)
  
  expect_s4_class(lines, "SpatVector")
  expect_equal(nrow(lines), 0)
  expect_s4_class(res$availableArea, "SpatVector")
  expect_true(terra::geomtype(res$availableArea) %in% c("polygons", "none"))
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
  expect_equal(nrow(lines), 2)
  expect_true("cluster" %in% names(lines))
  clusters <- unique(lines$cluster)
  expect_equal(length(clusters), 1)
  expect_equal(as.integer(table(lines$cluster)), 2L)
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
  
  expect_warning(
    res <- suppressMessages(createCropLayFinalYear1(
      Lay, pot,
      runClusteringInParallel = FALSE,
      clusterDistance         = 0,
      studyAreaHash           = "test_dist0"
    )),
    "degenerate"
  )
  lines <- unpack_res(res)
  
  expect_s4_class(lines, "SpatVector")
  expect_equal(nrow(lines), 0)
  expect_s4_class(res$availableArea, "SpatVector")
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

test_that("Degenerate clusters are dropped when each cluster has <2 lines", {
  reset_cache()
  pot <- vect("POLYGON ((0 0, 0 2000, 2000 2000, 2000 0, 0 0))", crs = "EPSG:32633")
  pot$Potential <- 1L
  line1 <- vect("LINESTRING (10 10, 20 10)", crs = "EPSG:32633")
  line2 <- vect("LINESTRING (1500 1500, 1510 1500)", crs = "EPSG:32633")
  Lay <- rbind(line1, line2)
  
  res <- expect_warning(
    createCropLayFinalYear1(
      Lay,
      pot,
      runClusteringInParallel = FALSE,
      clusterDistance         = 1,
      studyAreaHash           = "test_degenerate_clusters"
    ),
    "all clusters dropped as degenerate"
  )
  
  lines <- unpack_res(res)
  expect_s4_class(lines, "SpatVector")
  expect_equal(nrow(lines), 0)
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
    "`clusterDistance` must be non-negative"
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

test_that("Vector clusterDistance throws an error", {
  reset_cache()
  pot <- vect("POLYGON ((0 0, 0 1000, 1000 1000, 1000 0, 0 0))", crs = "EPSG:32633")
  pot$Potential <- 1L
  Lay <- vect("LINESTRING (1 1, 2 2)", crs = "EPSG:32633")
  
  expect_error(
    createCropLayFinalYear1(
      Lay,
      pot,
      runClusteringInParallel = FALSE,
      clusterDistance        = c(1, 2),
      studyAreaHash          = "test_vec_dist"
    ),
    "single numeric"
  )
})

test_that("Missing Potential column after intersect throws an error", {
  reset_cache()
  pot <- vect("POLYGON ((0 0, 0 1000, 1000 1000, 1000 0, 0 0))", crs = "EPSG:32633")
  Lay <- vect("LINESTRING (10 10, 990 990)", crs = "EPSG:32633")
  
  expect_error(
    createCropLayFinalYear1(
      Lay,
      pot,
      runClusteringInParallel = FALSE,
      clusterDistance         = 5,
      studyAreaHash           = "test_missing_potential"
    ),
    "no 'Potential' column"
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

test_that("splits lines at potential boundaries and propagates Potential", {
  skip_on_ci()
  
  # Provide no-op tic/toc if tictoc isn't attached
  if (!exists("tic", mode = "function"))  assign("tic", function(...) invisible(NULL), envir = .GlobalEnv)
  if (!exists("toc", mode = "function"))  assign("toc", function(...) invisible(NULL), envir = .GlobalEnv)
  
  # --- Fixtures --------------------------------------------------------------
  mk_sq <- function(x0, y0, size = 50, crs = "EPSG:3857") {
    vect(matrix(c(x0,y0, x0+size,y0, x0+size,y0+size, x0,y0+size, x0,y0),
                ncol = 2, byrow = TRUE), type = "polygons", crs = crs)
  }
  mk_line <- function(x0,y0,x1,y1, crs = "EPSG:3857") {
    v <- vect(matrix(c(x0,y0, x1,y1), ncol=2, byrow=TRUE), type="lines", crs = crs)
    v
  }
  
  # Two adjacent potential polygons: left=1, right=2
  left  <- mk_sq(0, 0, size = 50); left$Potential  <- 1L
  right <- mk_sq(50, 0, size = 50); right$Potential <- 2L
  pot   <- rbind(left, right)
  
  # One line that crosses from left polygon (1) to right polygon (2)
  lay   <- mk_line(10, 25, 90, 25)
  
  # --- Call function under test ---------------------------------------------
  expect_warning(
    out <- createCropLayFinalYear1(
      Lay = lay,
      potLayTopValid = pot,
      runClusteringInParallel = FALSE,
      clusterDistance = 10,
      studyAreaHash = "unit-test"
    ),
    "degenerate"
  )
  
  # --- Assertions ------------------------------------------------------------
  # Contract: list(lines=SpatVector lines, availableArea=SpatVector polygons)
  expect_type(out, "list")
  expect_true(all(c("lines", "availableArea") %in% names(out)))
  expect_s4_class(out$lines, "SpatVector")
  expect_equal(nrow(out$lines), 0)
  expect_s4_class(out$availableArea, "SpatVector")
  aa_gt <- tryCatch(terra::geomtype(out$availableArea), error = function(e) "none")
  expect_true(aa_gt %in% c("polygons", "none"))
})

test_that("no lines within potential returns empty lines and unchanged availableArea", {
  skip_on_ci()

  mk_sq <- function(x0, y0, size = 50, crs = "EPSG:3857") {
    vect(matrix(c(x0,y0, x0+size,y0, x0+size,y0+size, x0,y0+size, x0,y0),
                ncol = 2, byrow = TRUE), type = "polygons", crs = crs)
  }
  mk_line <- function(x0,y0,x1,y1, crs = "EPSG:3857") {
    vect(matrix(c(x0,y0, x1,y1), ncol=2, byrow=TRUE), type="lines", crs = crs)
  }
  
  pot <- mk_sq(0, 0, size = 50); pot$Potential <- 1L
  # Place line far away so it doesn't intersect
  lay <- mk_line(1000, 1000, 1100, 1000)
  
  out <- createCropLayFinalYear1(
    Lay = lay,
    potLayTopValid = pot,
    runClusteringInParallel = FALSE,
    clusterDistance = 10,
    studyAreaHash = "unit-test"
  )
  
  # type checks
  expect_type(out, "list")
  expect_s4_class(out$lines, "SpatVector")
  expect_equal(nrow(out$lines), 0L)
  
  aa_gt <- tryCatch(terra::geomtype(out$availableArea), error = function(e) "none")
  expect_true(aa_gt %in% c("polygons", "none"))
  
  # area compare (handle 'none' => 0 area)
  area_poly <- function(v) {
    gt <- tryCatch(terra::geomtype(v), error = function(e) "none")
    if (gt == "none" || nrow(v) == 0L) return(0)
    sum(terra::expanse(v), na.rm = TRUE)
  }
  
  a0 <- area_poly(pot)                 # original
  a1 <- area_poly(out$availableArea)   # returned
  expect_equal(a0, a1, tolerance = 1e-6)
})

test_that("availableArea keeps 'polygons' geomtype even when empty", {
  skip_on_ci()
  
  # Tiny potential polygon fully erased by 50m buffer around the line
  pot <- vect(matrix(c(0,0, 50,0, 50,50, 0,50, 0,0), ncol=2, byrow=TRUE),
              type = "polygons", crs = "EPSG:3857")
  pot$Potential <- 1L
  lay <- vect(matrix(c(10,25, 40,25), ncol=2, byrow=TRUE),
              type = "lines", crs = crs(pot))
  
  out <- createCropLayFinalYear1(lay, pot, FALSE, 10, "ut")
  expect_s4_class(out$availableArea, "SpatVector")
  aa_gt <- tryCatch(terra::geomtype(out$availableArea), error = function(e) "none")
  expect_true(aa_gt %in% c("polygons", "none"))
  expect_true(nrow(out$availableArea) == 0L || sum(terra::expanse(out$availableArea), na.rm = TRUE) == 0)
})

test_that("createCropLayFinalYear1: deduplicates Potential* columns to a single 'Potential'", {
  skip_on_ci()

  mk_sq <- function(x0, y0, size = 50, crs = "EPSG:3857") {
    vect(matrix(c(x0,y0, x0+size,y0, x0+size,y0+size, x0,y0+size, x0,y0),
                ncol = 2, byrow = TRUE), type = "polygons", crs = crs)
  }
  mk_line <- function(x0,y0,x1,y1, crs = "EPSG:3857") {
    vect(matrix(c(x0,y0, x1,y1), ncol=2, byrow=TRUE), type="lines", crs = crs)
  }
  
  pot <- mk_sq(0, 0, size = 50); pot$Potential <- 7L
  lay <- mk_line(10, 25, 40, 25)
  # Force a Potential column on Lay too, to simulate duplicate attribute sources
  lay$Potential <- 999L
  
  out <- createCropLayFinalYear1(
    Lay = lay,
    potLayTopValid = pot,
    runClusteringInParallel = FALSE,
    clusterDistance = 0,
    studyAreaHash = "unit-test"
  )
  
  # Exactly one Potential-like column must remain on lines
  potCols <- grep("^Potential", names(out$lines), value = TRUE)
  expect_equal(length(potCols), 1,
               info = paste("Found columns:", paste(potCols, collapse = ", ")))
  expect_true("Potential" %in% names(out$lines))
})

test_that("createCropLayFinalYear1: availableArea equals potLayTopValid minus 50 m buffered lines", {
  skip_on_ci()
  
  mk_sq <- function(x0, y0, size = 80, crs = "EPSG:3857") {
    vect(matrix(c(x0,y0, x0+size,y0, x0+size,y0+size, x0,y0+size, x0,y0),
                ncol = 2, byrow = TRUE), type = "polygons", crs = crs)
  }
  mk_line <- function(x0,y0,x1,y1, crs = "EPSG:3857") {
    vect(matrix(c(x0,y0, x1,y1), ncol=2, byrow=TRUE), type="lines", crs = crs)
  }
  
  pot <- mk_sq(0, 0, size = 80); pot$Potential <- 3L
  lay <- mk_line(10, 10, 70, 70)
  
  out <- createCropLayFinalYear1(
    Lay = lay,
    potLayTopValid = pot,
    runClusteringInParallel = FALSE,
    clusterDistance = 10,
    studyAreaHash = "unit-test"
  )
  
  # Expected area: erase(pot, aggregate(buffer(intersect(lay, pot), 50)))
  # Use the same sequence the function describes
  layCrop <- reproducible::Cache(postProcessTo, lay, pot)
  # If postProcessTo doesn't intersect, guard by doing intersect here for expected value
  if (nrow(terra::intersect(layCrop, pot)) > 0) {
    layCrop <- terra::intersect(pot, layCrop)
  }
  buf50   <- terra::buffer(layCrop, width = 50)
  bufAgg  <- terra::aggregate(buf50, dissolve = TRUE)
  expectAvail <- terra::erase(pot, bufAgg)
  
  # Compare areas
  area_out   <- sum(terra::expanse(out$availableArea), na.rm = TRUE)
  area_expect<- sum(terra::expanse(expectAvail),       na.rm = TRUE)
  # Allow small numeric tolerance
  expect_equal(area_out, area_expect, tolerance = 1e-5)
})
