# Single include pattern should behave like grepl
x <- c("apple", "banana", "apricot", "orange")
test_that("single include pattern works like grepl", {
  expect_equal(
    grepMulti(x, patterns = "ap"),
    c("apple", "apricot")
  )
})

# Multiple include patterns: element must contain all patterns
x2 <- c("foobar", "barfoo", "foobarbaz")
test_that("multiple include patterns require all patterns", {
  expect_equal(
    grepMulti(x2, patterns = c("foo", "bar")),
    c("foobar", "barfoo", "foobarbaz")
  )
  expect_equal(
    grepMulti(x2, patterns = c("foo", "baz")),
    c("foobarbaz")
  )
})

# No matches returns empty character vector
test_that("no include matches returns empty", {
  expect_equal(
    grepMulti(x, patterns = "z"),
    character(0)
  )
})

# Unwanted patterns: single pattern excludes elements containing it
x3 <- c("apple", "apricot", "application")
test_that("single unwanted pattern filters out correctly", {
  expect_equal(
    grepMulti(x3, patterns = "app", unwanted = "ric"),
    c("apple", "application")
  )
})

# Multiple unwanted patterns: only elements matching all unwanted are removed
test_that("multiple unwanted patterns require all to exclude", {
  expect_equal(
    grepMulti(x2, patterns = "foo", unwanted = c("bar", "baz")),
    c("foobar", "barfoo")
  )
})

