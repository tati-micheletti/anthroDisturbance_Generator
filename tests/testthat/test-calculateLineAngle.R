
test_that("horizontal line to the right has angle 0°", {
  pts <- matrix(c(0,0, 10,0), ncol = 2, byrow = TRUE)
  l <- vect(pts, type = "lines")
  expect_equal(unname(calculateLineAngle(l)), 0)
})

test_that("horizontal line to the left has angle 180° or -180°", {
  pts <- matrix(c(10,0, 0,0), ncol = 2, byrow = TRUE)
  l <- vect(pts, type = "lines")
  angle <- calculateLineAngle(l)
  expect_true(angle == 180 || angle == -180)
})

test_that("vertical line upward has angle 90°", {
  pts <- matrix(c(0,0, 0,10), ncol = 2, byrow = TRUE)
  l <- vect(pts, type = "lines")
  expect_equal(unname(calculateLineAngle(l)), 90)
})

test_that("vertical line downward has angle -90°", {
  pts <- matrix(c(0,10, 0,0), ncol = 2, byrow = TRUE)
  l <- vect(pts, type = "lines")
  expect_equal(unname(calculateLineAngle(l)), -90)
})

test_that("45° diagonal line", {
  pts <- matrix(c(0,0, 10,10), ncol = 2, byrow = TRUE)
  l <- vect(pts, type = "lines")
  expect_equal(unname(calculateLineAngle(l)), 45)
})

test_that("−135° diagonal (down-left)", {
  pts <- matrix(c(0,0, -10,-10), ncol = 2, byrow = TRUE)
  l <- vect(pts, type = "lines")
  expect_equal(unname(calculateLineAngle(l)), -135)
})

test_that("zero-length line returns 0° (atan2(0,0) == 0)", {
  pts <- matrix(c(5,5, 5,5), ncol = 2, byrow = TRUE)
  l <- vect(pts, type = "lines")
  expect_equal(unname(calculateLineAngle(l)), 0)
})

test_that("result is a plain numeric of length 1 with no names", {
  pts <- matrix(c(0,0, 0,1), ncol=2, byrow=TRUE)
  l <- vect(pts, type="lines")
  out <- unname(calculateLineAngle(l))
  expect_type(out, "double")
  expect_length(out, 1)
  expect_null(names(out))
})

test_that("no error when input has fewer than two vertices", {
  pts <- matrix(c(0,0), ncol = 2)
  l <- try(vect(pts, type="lines"), silent = TRUE)
  expect_silent(calculateLineAngle(l))  # accept any error for invalid input
})

test_that("line with more than two points uses only first segment's angle", {
  # First segment from (0,0) to (1, tan(30°)) → angle ≈ 30°
  # Second segment goes somewhere else (we choose 60°) but it should be ignored.
  pts <- matrix(c(
    0, 0,
    1, tan(pi/6),                # 30° segment
    1 + cos(pi/3), tan(pi/6) + sin(pi/3)  # second segment at 60°
  ), ncol = 2, byrow = TRUE)
  l <- vect(pts, type = "lines")
  expect_equal(
    unname(calculateLineAngle(l)),
    30,
    tolerance = 1e-6
  )
})

test_that("random two-point lines reproducibly give the expected angle", {
  set.seed(123)
  for(i in seq_len(5)) {
    x1 <- runif(1) 
    y1 <- runif(1)
    ang0 <- runif(1, -180, 180)
    len <- runif(1, 0.1, 5)
    x2 <- x1 + len * cos(ang0 * pi/180)
    y2 <- y1 + len * sin(ang0 * pi/180)
    l <- vect(matrix(c(x1,y1, x2,y2), ncol=2, byrow=TRUE), type="lines")
    expect_equal(unname(calculateLineAngle(l)), ang0, tolerance = 1e-6)
  }
})
