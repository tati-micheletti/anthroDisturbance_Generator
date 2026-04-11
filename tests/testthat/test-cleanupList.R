# Unit tests for cleanupList()

# Helper to clone a list deeply
copy_list <- function(x) { unserialize(serialize(x, NULL)) }

# 1. Test inner removal of NULLs
test_that("inner removes NULL elements in sublists", {
  input <- list(
    a = list(1, NULL, 2),
    b = list(NULL, 3)
  )
  result <- cleanupList(input, inner = TRUE, outter = FALSE,
                        cleanEmpty = FALSE, nullEmpty = FALSE)
  expect_equal(result$a, list(1, 2))
  expect_equal(result$b, list(3))
})

# 2. Test skipping inner when inner = FALSE
test_that("inner flag FALSE leaves sublists intact", {
  input <- list(
    a = list(1, NULL, 2),
    b = list(NULL, 3)
  )
  result <- cleanupList(input, inner = FALSE, outter = FALSE,
                        cleanEmpty = FALSE, nullEmpty = FALSE)
  expect_equal(result$a, input$a)
  expect_equal(result$b, input$b)
})

# 3. Test outter removal of NULLs
test_that("outter removes NULL top-level elements", {
  input <- list(
    a = 1,
    b = NULL,
    c = 3
  )
  result <- cleanupList(input, outter = TRUE, inner = FALSE,
                        cleanEmpty = FALSE, nullEmpty = FALSE)
  expect_false("b" %in% names(result))
  expect_equal(result$a, 1)
  expect_equal(result$c, 3)
})

# 4. Test nullEmpty converts empty lists to NULL
test_that("nullEmpty converts empty lists to NULL", {
  input <- list(
    a = list(),
    b = list(1)
  )
  result <- cleanupList(input, outter = FALSE, inner = FALSE,
                        cleanEmpty = FALSE, nullEmpty = TRUE)
  expect_null(result$a)
  expect_equal(result$b, list(1))
})

# 5. Test cleanEmpty removes empty lists
test_that("cleanEmpty removes empty list elements", {
  input <- list(
    a = list(),
    b = list(1),
    c = list()
  )
  result <- cleanupList(input, outter = FALSE, inner = FALSE,
                        cleanEmpty = TRUE, nullEmpty = FALSE)
  expect_false("a" %in% names(result))
  expect_false("c" %in% names(result))
  expect_equal(result$b, list(1))
})

# 6. Combined flags: inner and outter together
test_that("combined inner and outter flags", {
  input <- list(
    a = list(1, NULL),
    b = NULL,
    c = list()
  )
  result <- cleanupList(input,
                        inner = TRUE,
                        outter = TRUE,
                        cleanEmpty = TRUE,
                        nullEmpty = TRUE)
  # a: inner removes NULL, becomes list(1), kept; c: empty list converted to NULL then removed
  expect_equal(result$a, list(1))
  expect_false("b" %in% names(result))
  expect_false("c" %in% names(result))
})

# 7. No flags: identical copy
test_that("no flags returns identical copy", {
  input <- list(a = list(), b = NULL, c = list(5))
  result <- cleanupList(input, inner = FALSE, outter = FALSE,
                        cleanEmpty = FALSE, nullEmpty = FALSE)
  expect_identical(result, input)
})
