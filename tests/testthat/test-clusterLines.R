# Unit tests for clusterLines() using testthat

library(testthat)
library(terra)
library(doParallel)
library(foreach)

# 1. Test that clustering with a very large distThreshold assigns all lines to a single cluster

# Create a simple SpatVector of two distant lines
line1_matrix <- matrix(c(0, 1, 0, 0), ncol = 2, byrow = TRUE)
line2_matrix <- matrix(c(10000, 10001, 10000, 10000), ncol = 2, byrow = TRUE)

line1 <- vect(line1_matrix, type = "lines")
line2 <- vect(line2_matrix, type = "lines")
lines <- rbind(line1, line2)

lines$Potential <- c(1, 1)
lines$individualID <- c("l1", "l2")
crs(lines) <- "EPSG:4326"

# 1. Test that clustering with a very large distThreshold assigns all lines to a single cluster
res_all <- clusterLines(lines, distThreshold = 1e6, currPotential = 1, totPotential = 1, runInParallel = FALSE)

test_that("All lines in one cluster for large threshold", {
  expect_equal(length(unique(res_all$cluster)), 1)
})

# 2. Test that clustering with a tiny distThreshold assigns each line to its own cluster
res_none <- clusterLines(lines, distThreshold = 0.001, currPotential = 1, totPotential = 1, runInParallel = FALSE)

test_that("Each line in its own cluster for tiny threshold", {
  expect_equal(length(unique(res_none$cluster)), 2)
})

# 3. Test that calculatedLength matches perim() output for each line
res_len <- clusterLines(lines, distThreshold = 1e6, currPotential = 1, totPotential = 1, runInParallel = FALSE)

test_that("calculatedLength equals terra::perim output", {
  true_lengths <- sapply(seq_len(nrow(res_len)), function(i) perim(res_len[i, ]))
  expect_equal(res_len$calculatedLength, true_lengths)
})

# 4. Test that angles are numeric and finite
res_ang <- clusterLines(lines, distThreshold = 1e6, currPotential = 1, totPotential = 1, runInParallel = FALSE)

test_that("Angles are numeric and finite", {
  expect_true(is.numeric(res_ang$angles))
  expect_true(all(is.finite(res_ang$angles)))
})

# 5. Test that Pot_Clus column correctly combines Potential and cluster ID
res_pot <- clusterLines(lines, distThreshold = 1e6, currPotential = 1, totPotential = 1, runInParallel = FALSE)

unique_clusters <- unique(res_pot$cluster)
expected <- paste0(1, "_", unique_clusters)
res_vals <- unique(res_pot$Pot_Clus)

test_that("Pot_Clus matches Potential_clusterID format", {
  expect_equal(sort(res_vals), sort(expected))
})

# Edge case tests

# 6. Empty SpatVector input: should error on no geometries
empty_lines <- vect()
empty_lines$individualID <- character(0)

test_that("Empty input returns empty output without error", {
  expect_silent({
    res_empty <- clusterLines(empty_lines, distThreshold=100,
                              currPotential=1, totPotential=1,
                              runInParallel=FALSE)
  })
  expect_true(nrow(res_empty) == 0)
})

# 7. Identical-line geometries cluster together
line_a_matrix <- matrix(c(0, 0, 5, 0), ncol = 2, byrow = TRUE)
line_b_matrix <- line_a_matrix
line_a <- vect(line_a_matrix, type = "lines")
line_b <- vect(line_b_matrix, type = "lines")
lines_identical <- rbind(line_a, line_b)
lines_identical$Potential <- c(2, 2)
lines_identical$individualID <- c("a", "b")
crs(lines_identical) <- "EPSG:4326"

test_that("Duplicate lines cluster as one when threshold >= 0", {
  res_ident <- clusterLines(lines_identical, distThreshold=0,
                            currPotential=2, totPotential=2,
                            runInParallel=FALSE)
  expect_equal(length(unique(res_ident$cluster)), 1)
})

# 8. Parallel mode yields same clusters as sequential
test_that("Parallel and sequential produce identical cluster results", {
  testthat::skip_on_covr()
  skip_if(parallel::detectCores() < 2, "Need >=2 cores for parallel test")
  # create a small PSOCK cluster for reproducibility
  cl <- parallel::makeCluster(2, type="PSOCK")
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  # sequential result
  res_seq <- clusterLines(lines, distThreshold=1e6,
                          currPotential=1, totPotential=1,
                          runInParallel=FALSE)
  # parallel result
  res_par <- clusterLines(lines, distThreshold=1e6,
                          currPotential=1, totPotential=1,
                          runInParallel=TRUE)
  expect_equal(res_par$cluster, res_seq$cluster)
  expect_equal(res_par$Pot_Clus, res_seq$Pot_Clus)
})

# 9. plotClusterDiagnostic does not error
test_that("plotClusterDiagnostic runs cleanly", {
  expect_warning(
    clusterLines(lines, distThreshold=1e6,
                 currPotential=1, totPotential=1,
                 plotClusterDiagnostic=TRUE,
                 runInParallel=FALSE), NA
  )
})

