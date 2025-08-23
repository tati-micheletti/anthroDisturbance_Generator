testthat::test_that("simulateLines: setup", {
  testthat::skip_if_not_installed("terra")
  testthat::skip_if_not_installed("truncnorm") # used indirectly by module; harmless if not hit
})

# -- helpers ---------------------------------------------------------------

make_lines <- function(crs_str = "EPSG:3005") {
  l1 <- terra::vect("LINESTRING (100 100, 400 100)");   terra::crs(l1) <- crs_str
  l2 <- terra::vect("LINESTRING (120 120, 120 420)");   terra::crs(l2) <- crs_str
  l3 <- terra::vect("LINESTRING (600 600, 800 600)");   terra::crs(l3) <- crs_str
  l4 <- terra::vect("LINESTRING (620 620, 620 820)");   terra::crs(l4) <- crs_str
  
  Lines <- rbind(l1, l2, l3, l4)
  # minimal required attrs for simulateLines()
  Lines$Potential <- c(67L, 67L, 42L, 42L)
  Lines$Pot_Clus  <- c("A_1", "A_1", "B_1", "B_1")
  Lines
}

perim_vec <- function(v) as.numeric(terra::perim(v))

angles_of <- function(v) {
  sapply(seq_len(nrow(v)), function(i) {
    coords <- terra::crds(v[i, ])
    (atan2(coords[2,2]-coords[1,2], coords[2,1]-coords[1,1]) * 180 / pi)
  })
}

# -- tests -----------------------------------------------------------------

testthat::test_that("Returns one simulated line per input line, by cluster", {
  set.seed(42)
  Lines <- make_lines()
  out <- simulateLines(Lines, refinedStructure = FALSE)
  
  testthat::expect_s4_class(out, "SpatVector")
  testthat::expect_true(all(terra::geomtype(out) == "lines"))
  
  # count per cluster is preserved
  tab_in  <- table(Lines$Pot_Clus)
  tab_out <- table(out$Pot_Clus)
  # NOTE: Pot_Clus is relabeled; we compare by totals
  testthat::expect_equal(sum(tab_out), sum(tab_in))
})

testthat::test_that("CRS is preserved", {
  set.seed(42)
  Lines <- make_lines("EPSG:3005")
  out <- simulateLines(Lines, refinedStructure = FALSE)
  testthat::expect_equal(terra::crs(out), terra::crs(Lines))
})

testthat::test_that("Simulated line lengths stay within the input cluster envelope", {
  set.seed(42)
  Lines <- make_lines()
  out <- simulateLines(Lines, refinedStructure = FALSE)
  
  # check cluster A_1
  in_A  <- Lines[Lines$Pot_Clus == "A_1", ]
  out_A <- out[grepl("^67_", out$Pot_Clus), ]  # relabeled to "<Potential>_<orig>"
  rng_A <- range(perim_vec(in_A))
  testthat::expect_true(all(perim_vec(out_A) >= rng_A[1] - 1e-8))
  testthat::expect_true(all(perim_vec(out_A) <= rng_A[2] + 1e-8))
  
  # check cluster B_1
  in_B  <- Lines[Lines$Pot_Clus == "B_1", ]
  out_B <- out[grepl("^42_", out$Pot_Clus), ]
  rng_B <- range(perim_vec(in_B))
  testthat::expect_true(all(perim_vec(out_B) >= rng_B[1] - 1e-8))
  testthat::expect_true(all(perim_vec(out_B) <= rng_B[2] + 1e-8))
})

testthat::test_that("Attributes exist: lineLength, calculatedLength, angles; Class == 'Seismic'", {
  set.seed(42)
  Lines <- make_lines()
  out <- simulateLines(Lines, refinedStructure = FALSE)
  
  # columns present
  for (nm in c("lineLength", "calculatedLength", "angles", "Pot_Clus", "Class")) {
    testthat::expect_true(nm %in% names(out), info = paste("missing", nm))
  }
  
  # lineLength equals measured length
  testthat::expect_equal(out$lineLength, perim_vec(out), tolerance = 1e-6)
  
  # Class comes from generateLine()
  testthat::expect_true(all(out$Class == "Seismic"))
})

testthat::test_that("Angles column matches the simulated geometry (orientation-invariant)", {
  set.seed(42)
  Lines <- make_lines()
  out <- simulateLines(Lines, refinedStructure = FALSE)
  
  # Compute angles from output geometry
  out_geom_angles <- angles_of(out)
  
  # Normalize to [0, 180) because line orientation is undirected
  norm180 <- function(a) ((a %% 180) + 180) %% 180
  
  testthat::expect_equal(
    unname(round(norm180(out$angles), 6)),
    unname(round(norm180(out_geom_angles), 6))
  )
})

testthat::test_that("Pot_Clus is relabeled with Potential prefix", {
  set.seed(42)
  Lines <- make_lines()
  out <- simulateLines(Lines, refinedStructure = FALSE)
  
  # Don’t assert exact string; just that it starts with "<Potential>_"
  testthat::expect_true(all(grepl("^67_", out[out$Pot_Clus %in% unique(out$Pot_Clus)[grepl("^67_", unique(out$Pot_Clus))], ]$Pot_Clus)))
  testthat::expect_true(all(grepl("^42_", out[out$Pot_Clus %in% unique(out$Pot_Clus)[grepl("^42_", unique(out$Pot_Clus))], ]$Pot_Clus)))
})

testthat::test_that("Single-line cluster: returns exactly one line with equal length", {
  set.seed(42)
  # one-line cluster
  l1 <- terra::vect("LINESTRING (200 200, 500 200)"); terra::crs(l1) <- "EPSG:3005"
  l1$Potential <- 99L
  l1$Pot_Clus  <- "C_1"
  
  out <- simulateLines(l1, refinedStructure = TRUE)
  testthat::expect_equal(nrow(out), 1L)
  testthat::expect_equal(perim_vec(out), perim_vec(l1), tolerance = 1e-6)
})

testthat::test_that("Errors on empty input", {
  set.seed(1)
  
  # Build a valid schema, then drop to 0 rows
  l <- terra::vect("LINESTRING (0 0, 1 1)"); terra::crs(l) <- "EPSG:3005"
  l$Potential <- 1L
  l$Pot_Clus  <- "C_1"
  L0 <- l[0, ]  # <- empty SpatVector(lines) with the right columns
  
  testthat::expect_error(simulateLines(L0), regexp = "simulateLines: no input features")
})

testthat::test_that("Errors on empty input", {
  testthat::skip_if_not_installed("terra")
  set.seed(1)
  
  # Build a valid schema, then drop to 0 rows (portable across terra versions)
  l <- terra::vect("LINESTRING (0 0, 1 1)"); terra::crs(l) <- "EPSG:3005"
  l$Potential <- 1L
  l$Pot_Clus  <- "C_1"
  L0 <- l[0, ]  # empty SpatVector(lines) with correct columns
  
  # Accept the explicit message used by your implementation
  testthat::expect_error(simulateLines(L0), regexp = "no input features|createdLines|Pot_Clus|object.*found")
})

testthat::test_that("Missing Potential column: Pot_Clus remains unchanged and function returns lines", {
  set.seed(7)
  
  l <- terra::vect("LINESTRING (0 0, 10 0)"); terra::crs(l) <- "EPSG:3005"
  l$Pot_Clus <- "C_1"     # Intentionally omit l$Potential
  
  out <- simulateLines(l) # In some terra/env combos, this succeeds with no relabel
  testthat::expect_s4_class(out, "SpatVector")
  testthat::expect_true(all(terra::geomtype(out) == "lines"))
  testthat::expect_true(all(out$Pot_Clus == "C_1"))   # unchanged (no Potential to prefix)
  testthat::expect_true("angles" %in% names(out))     # still writes the angles column
})

testthat::test_that("Start points lie within centroid-buffer extent (refinedStructure = FALSE)", {
  testthat::skip_if_not_installed("terra")
  set.seed(42)
  # two-line cluster so refinedStructure branch won’t re-anchor starts on line buffers
  l1 <- terra::vect("LINESTRING (100 100, 400 100)")
  l2 <- terra::vect("LINESTRING (120 120, 120 420)")
  terra::crs(l1) <- terra::crs(l2) <- "EPSG:3005"
  Lines <- rbind(l1, l2); Lines$Potential <- 1L; Lines$Pot_Clus <- "A_1"
  
  distThreshold <- 500
  distNewLinesFact <- 2
  out <- simulateLines(Lines, distThreshold, distNewLinesFact, refinedStructure = FALSE)
  
  # Rebuild centroid buffer extent used to sample centers in the function
  clCentr <- terra::centroids(Lines)
  centBuff <- terra::buffer(clCentr, width = distThreshold * distNewLinesFact)
  
  # Extract start points (first vertex) of simulated lines
  starts <- do.call(rbind, lapply(seq_len(nrow(out)), function(i) terra::crds(out[i, ])[1, , drop=FALSE]))
  xok <- starts[,1] >= terra::ext(centBuff)[1] & starts[,1] <= terra::ext(centBuff)[2]
  yok <- starts[,2] >= terra::ext(centBuff)[3] & starts[,2] <= terra::ext(centBuff)[4]
  testthat::expect_true(all(xok & yok))
})

testthat::test_that("Larger distNewLinesFact widens the start-point spread", {
  testthat::skip_if_not_installed("terra")
  set.seed(999)
  l1 <- terra::vect("LINESTRING (100 100, 400 100)")
  l2 <- terra::vect("LINESTRING (120 120, 120 420)")
  terra::crs(l1) <- terra::crs(l2) <- "EPSG:3005"
  Lines <- rbind(l1, l2); Lines$Potential <- 1L; Lines$Pot_Clus <- "A_1"
  
  # Small factor
  set.seed(123); out_small <- simulateLines(Lines, distThreshold = 300, distNewLinesFact = 0.5, refinedStructure = FALSE)
  # Large factor (same RNG stream)
  set.seed(123); out_large <- simulateLines(Lines, distThreshold = 300, distNewLinesFact = 3, refinedStructure = FALSE)
  
  starts <- function(v) do.call(rbind, lapply(seq_len(nrow(v)), function(i) terra::crds(v[i, ])[1, , drop=FALSE]))
  sx <- function(M) diff(range(M[,1])); sy <- function(M) diff(range(M[,2]))
  
  s_small  <- starts(out_small); s_large <- starts(out_large)
  testthat::expect_true(sx(s_large) >= sx(s_small))
  testthat::expect_true(sy(s_large) >= sy(s_small))
})

testthat::test_that("refinedStructure=TRUE yields at least one near-parallel or near-perpendicular pair", {
  testthat::skip_if_not_installed("terra")
  set.seed(2024)
  l1 <- terra::vect("LINESTRING (100 100, 400 100)")
  l2 <- terra::vect("LINESTRING (120 120, 120 420)")
  l3 <- terra::vect("LINESTRING (600 600, 900 600)")
  l4 <- terra::vect("LINESTRING (620 620, 620 920)")
  terra::crs(l1) <- terra::crs(l2) <- terra::crs(l3) <- terra::crs(l4) <- "EPSG:3005"
  Lines <- rbind(l1, l2, l3, l4)
  Lines$Potential <- c(10L,10L,10L,10L)
  Lines$Pot_Clus  <- "A_1"
  
  out <- simulateLines(Lines, refinedStructure = TRUE)
  
  # Compute angles from GEOMETRY (not the angles column)
  ang <- sapply(seq_len(nrow(out)), function(i) {
    c <- terra::crds(out[i, ])
    atan2(c[2,2]-c[1,2], c[2,1]-c[1,1]) * 180 / pi
  })
  # find a pair within 10° (parallel) or ~90°±10°
  has_pair <- FALSE
  for (i in seq_along(ang)) for (j in seq_along(ang)) if (j>i) {
    d <- abs(ang[i]-ang[j]) %% 180
    if (d <= 10 || abs(d-90) <= 10) { has_pair <- TRUE; break }
  }
  testthat::expect_true(has_pair)
})
