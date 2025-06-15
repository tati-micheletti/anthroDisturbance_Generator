# Unit tests for wrapTerraList and unwrapTerraList

library(testthat)
library(terra)
library(qs)
library(stringi)
library(zip)
library(Require)
library(mockery)

# Helper: create a simple SpatVector for testing
create_dummy_vect <- function() {
  terra::vect(matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE), type = "points")
}


# ---- Basic tests ----

test_that("wrapTerraList wraps and saves SpatVectors", {
  tmp <- tempfile(); dir.create(tmp)
  on.exit(unlink(tmp, recursive=TRUE), add=TRUE)
  tv <- list(g1=list(a=create_dummy_vect()), g2=list(b=create_dummy_vect()))
  res <- suppressMessages(wrapTerraList(tv, generalPath=tmp, zipFiles=FALSE))
  
  expect_named(res, names(tv))
  for(grp in names(tv)){
    expect_named(res[[grp]], names(tv[[grp]]))
    for(feat in names(tv[[grp]])){
      pth <- res[[grp]][[feat]]
      expect_true(file.exists(pth))
      wrapped <- qread(pth)
      rec <- terra::vect(wrapped)
      expect_equal(as.matrix(rec), as.matrix(tv[[grp]][[feat]]))
    }
  }
})

test_that("unwrapTerraList reconstructs SpatVectors from paths", {
  tmp <- tempfile(); dir.create(tmp)
  on.exit(unlink(tmp, recursive=TRUE), add=TRUE)
  tv <- list(x=list(p=create_dummy_vect()))
  wrapped <- wrapTerraList(tv, generalPath=tmp, zipFiles=FALSE)
  rec <- unwrapTerraList(wrapped, generalPath=tmp)
  
  expect_s4_class(rec$x$p, "SpatVector")
  expect_equal(as.matrix(rec$x$p), as.matrix(tv$x$p))
})

# ---- Zip and Drive tests ----

test_that("wrapTerraList creates zip and manifest when zipFiles=TRUE", {
  tmp <- tempfile(); dir.create(tmp)
  on.exit(unlink(tmp, recursive=TRUE), add=TRUE)
  tv <- list(x=list(v=create_dummy_vect()))
  wrapTerraList(tv, generalPath=tmp, zipFiles=TRUE)
  expect_true(file.exists(file.path(tmp, "disturbanceList.zip")))
  expect_true(file.exists(file.path(tmp, "theList.qs")))
})

test_that("wrapTerraList uploads zip to Google Drive when uploadZip provided", {
  tmp <- tempfile(); dir.create(tmp)
  on.exit(unlink(tmp, recursive=TRUE), add=TRUE)
  tv <- list(a=list(b=create_dummy_vect()))
  # stub drive_upload
  mockery::stub(wrapTerraList, 'drive_upload', function(media, path) {
    expect_true(grepl("disturbanceList.zip$", media))
    expect_equal(path, as_id("drive123"))
    message("uploaded to drive123")
    TRUE
  })
  expect_message(
    wrapTerraList(tv, generalPath=tmp, zipFiles=TRUE, uploadZip="drive123"),
    "uploaded to drive123"
  )
})

test_that("unwrapTerraList downloads and unwraps when given a Drive link", {
  orig <- tempfile(); dir.create(orig)
  on.exit(unlink(orig, recursive=TRUE), add=TRUE)
  tv <- list(a=list(b=create_dummy_vect()))
  wrapTerraList(tv, generalPath=orig, zipFiles=TRUE)
  zipf <- file.path(orig, "disturbanceList.zip")
  
  dest <- tempfile(); dir.create(dest)
  on.exit(unlink(dest, recursive=TRUE), add=TRUE)
  
  # stub drive_download
  mockery::stub(unwrapTerraList, 'drive_download', function(file, path) {
    expect_equal(file, as_id("link123"))
    file.copy(zipf, file.path(dest, "disturbanceList.zip"))
    NULL
  })
  
  rec <- unwrapTerraList("link123", generalPath=dest)
  expect_s4_class(rec$a$b, "SpatVector")
  expect_equal(as.matrix(rec$a$b), as.matrix(tv$a$b))
})

# ---- Edge cases ----

test_that("wrapTerraList errors on empty terraList", {
  tmp <- tempfile(); dir.create(tmp)
  on.exit(unlink(tmp, recursive=TRUE), add=TRUE)
  expect_error(wrapTerraList(list(), generalPath=tmp, zipFiles=FALSE))
})

test_that("wrapTerraList errors on nested empty groups", {
  tmp <- tempfile(); dir.create(tmp)
  on.exit(unlink(tmp, recursive=TRUE), add=TRUE)
  tv <- list(g1=list(), g2=list(c=create_dummy_vect()))
  expect_error(wrapTerraList(tv, generalPath=tmp, zipFiles=FALSE))
})

test_that("unwrapTerraList errors if .qs files missing", {
  tmp <- tempfile(); dir.create(tmp)
  on.exit(unlink(tmp, recursive=TRUE), add=TRUE)
  tv <- list(a=list(b=create_dummy_vect()))
  wrapped <- wrapTerraList(tv, generalPath=tmp, zipFiles=FALSE)
  file.remove(unlist(wrapped))
  expect_error(unwrapTerraList(wrapped, generalPath=tmp))
})

test_that("wrapTerraList propagates zip errors", {
  tmp <- tempfile(); dir.create(tmp)
  on.exit(unlink(tmp, recursive=TRUE), add=TRUE)
  tv <- list(x=list(y=create_dummy_vect()))
  mockery::stub(wrapTerraList, 'zip', function(...) stop("zip failed"))
  expect_error(wrapTerraList(tv, generalPath=tmp, zipFiles=TRUE), "zip failed")
})

test_that("wrapTerraList handles special characters in generalPath", {
  tmp <- file.path(tempdir(), "sp c!@# ünicode")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive=TRUE), add=TRUE)
  tv <- list(a=list(b=create_dummy_vect()))
  res <- wrapTerraList(tv, generalPath=tmp, zipFiles=FALSE)
  expect_true(all(file.exists(unlist(res, use.names=FALSE))))
})

test_that("unwrapTerraList errors on invalid input", {
  tmp <- tempfile(); dir.create(tmp)
  on.exit(unlink(tmp, recursive=TRUE), add=TRUE)
  expect_error(unwrapTerraList(123, generalPath=tmp))
})
