# NA handling
test_that("generateLine handles NA inputs correctly", {
  # Test NA angle
  expect_error(generateLine(NA, 10, c(0, 10), c(0, 10), "EPSG:4326"), "angle is NA. Please debug")
  
  # Test NA xlim
  expect_error(generateLine(45, 10, NA, c(0, 10), "EPSG:4326"), "xlim is NA. Please debug")
  
  # Test NA ylim
  expect_error(generateLine(45, 10, c(0, 10), NA, "EPSG:4326"), "ylim is NA. Please debug")
  
  # Test NA crs
  expect_error(generateLine(45, 10, c(0, 10), c(0, 10), NA), "mCrs is NA or missing. Please debug")
  
  # Test NA length
  expect_error(generateLine(45, NA, c(0, 10), c(0, 10), "EPSG:4326"), "length is NA. Please debug")
})

# horizontal lines
test_that("generateLine produces correct horizontal line", {
  set.seed(123)
  mCrs <- "EPSG:4326"
  line <- generateLine(0, 10, c(0,100), c(0,100), mCrs)
  coords <- crds(line)
  n <- nrow(coords)
  
  dx <- unname(coords[n,1] - coords[1,1])
  dy <- unname(coords[n,2] - coords[1,2])
  
  # horizontal: dx ≈ 10, dy ≈ 0
  expect_equal(dx, 10, tolerance = 1e-6)
  expect_equal(dy,  0, tolerance = 1e-6)
})

# vertical lines
test_that("generateLine produces correct vertical line", {
  set.seed(456)
  mCrs <- "EPSG:4326"
  line <- generateLine(90, 20, c(0,100), c(0,100), mCrs)
  coords <- crds(line)
  n <- nrow(coords)
  
  dx <- unname(coords[n,1] - coords[1,1])
  dy <- unname(coords[n,2] - coords[1,2])
  
  # vertical: dx ≈ 0, dy ≈ 20
  expect_equal(dx,  0, tolerance = 1e-6)
  expect_equal(dy, 20, tolerance = 1e-6)
})

# diagonal tests
test_that("generateLine at 45° and 225°", {
  set.seed(1)
  mCrs <- "EPSG:4326"
  # force start-point to (0,0)
  line45 <- generateLine(45, sqrt(2), c(0,0), c(0,0), mCrs)
  coords45 <- crds(line45)
  dx45 <- unname(coords45[2,1] - coords45[1,1])
  dy45 <- unname(coords45[2,2] - coords45[1,2])
  expect_equal(dx45, 1, tolerance = 1e-6)
  expect_equal(dy45, 1, tolerance = 1e-6)
  
  line225 <- generateLine(225, sqrt(2), c(0,0), c(0,0), mCrs)
  coords225 <- crds(line225)
  dx225 <- unname(coords225[2,1] - coords225[1,1])
  dy225 <- unname(coords225[2,2] - coords225[1,2])
  expect_equal(dx225, -1, tolerance = 1e-6)
  expect_equal(dy225, -1, tolerance = 1e-6)
})

# zero-length
test_that("generateLine with length = 0 produces zero-length line", {
  mCrs <- "EPSG:4326"
  line0 <- generateLine(123, 0, c(5,5), c(7,7), mCrs)
  coords0 <- crds(line0)
  expect_equal(diff(coords0[,1]), 0, tolerance = 1e-6)
  expect_equal(diff(coords0[,2]), 0, tolerance = 1e-6)
})

# metadata checks
test_that("generateLine returns SpatVector with correct CRS and Class", {
  mCrs <- "EPSG:3857"
  line <- generateLine(0, 1, c(0,0), c(0,0), mCrs)
  
  expect_true(inherits(line, "SpatVector"))

  expect_equal(terra::crs(line, proj=TRUE),
               terra::crs(mCrs, proj=TRUE))

  expect_equal(line$Class, "Seismic")
})

# invalid-type errors
test_that("generateLine rejects non-numeric inputs", {
  mCrs <- "EPSG:4326"
  expect_error(generateLine("ninety", 10, c(0,100), c(0,100), mCrs), 
               "angle must be a single numeric value")
  expect_error(generateLine(0, 10, c("a","b"), c(0,100), mCrs), 
               "xlim must be a numeric vector of length 2")
  expect_error(generateLine(0, 10, c(0,100), c("c","d"), mCrs), 
               "ylim must be a numeric vector of length 2")
  expect_error(generateLine(0, -5, c(0,100), c(0,100), mCrs), 
               "length must be a single non-negative numeric")
})
