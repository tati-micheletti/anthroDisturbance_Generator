test_that("middle point for odd dimensions", {
  # 5 rows, 7 cols → trunc(5/2)=2, trunc(7/2)=3
  r <- rast(nrows = 5, ncols = 7)
  # fill in a known sequence so we can spot-check
  values(r) <- seq_len(ncell(r))
  expect_equal(
    getRasterMiddlePoint(r),
    terra::cellFromRowCol(r, 2, 3)
  )
})

test_that("middle point for even dimensions", {
  # 6 rows, 8 cols → trunc(6/2)=3, trunc(8/2)=4
  r <- rast(nrows = 6, ncols = 8)
  values(r) <- seq_len(ncell(r))
  expect_equal(
    getRasterMiddlePoint(r),
    terra::cellFromRowCol(r, 3, 4)
  )
})

test_that("works for 1×1 raster", {
  r <- rast(nrows = 1, ncols = 1)
  values(r) <- 42
  # trunc(1/2)=0, trunc(1/2)=0; terra::cellFromRowCol treats row=0 or col=0 as NA
  expect_true(is.na(getRasterMiddlePoint(r)))
})

test_that("tiny even raster 2×2 gives [1,1]", {
  r <- rast(nrows = 2, ncols = 2)
  expect_equal(getRasterMiddlePoint(r), terra::cellFromRowCol(r, 1, 1))
})