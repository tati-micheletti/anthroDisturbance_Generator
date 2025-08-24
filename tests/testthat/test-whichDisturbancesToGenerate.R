# Basic sequencing
test_that("picks correct indices for simple intervals", {
  params <- data.frame(disturbanceInterval = c(1, 2, 3))
  
  # 4-year: intervals 1 & 2 hit
  expect_equal(
    whichDisturbancesToGenerate(0, 10, 4, params),
    c(1L, 2L)
  )
  
  # 10-year: intervals 1 & 2 (3-year interval does not include 10)
  expect_equal(
    whichDisturbancesToGenerate(0, 10, 10, params),
    c(1L, 2L)
  )
  
  # 5-year: only the 1-year interval hits
  expect_equal(
    whichDisturbancesToGenerate(0, 10, 5, params),
    c(1L)
  )
  
  # 7-year: only the 1-year interval hits
  expect_equal(
    whichDisturbancesToGenerate(0, 10, 7, params),
    c(1L)
  )
})

# Edge cases at the bounds
test_that("includes startTime and endTime when on-sequence", {
  params <- data.frame(disturbanceInterval = c(2, 5))
  
  # startTime = 0: both intervals
  expect_equal(
    whichDisturbancesToGenerate(0, 10, 0, params),
    c(1L, 2L)
  )
  
  # endTime = 10: both 2- and 5-year intervals include 10
  expect_equal(
    whichDisturbancesToGenerate(0, 10, 10, params),
    c(1L, 2L)
  )
})

# Intervals larger than window
test_that("handles intervals larger than the simulation span", {
  params <- data.frame(disturbanceInterval = c(20, 15))
  
  # only startTime (0) falls
  expect_equal(
    whichDisturbancesToGenerate(0, 10, 0, params),
    c(1L, 2L)
  )
  
  # currentTime inside but no startTime => none
  expect_equal(
    whichDisturbancesToGenerate(0, 10, 10, params),
    integer(0)
  )
})

# Non-positive or NA intervals
test_that("returns empty for non-positive or NA intervals", {
  params <- data.frame(disturbanceInterval = c(0, -1, NA))
  idx <- suppressWarnings(whichDisturbancesToGenerate(0, 10, 5, params))
  expect_warning(whichDisturbancesToGenerate(0, 10, 5, params),
                 "Skipping 3 row\\(s\\) with non-positive or non-numeric disturbanceInterval: 1, 2, 3")
  expect_equal(idx, integer(0))
})

