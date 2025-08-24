# Disable crayon for clean message testing
withr::local_options(list(crayon.enabled = FALSE))

# 1. Precondition checks
test_that("Precondition checks work", {
  # Both sim/pathInput missing
  expect_error(
    createModObject(data = "x", sim = NULL, pathInput = NULL, currentTime = 1),
    "Either a simList or a folder containing the data need to be supplied"
  )
  
  # Invalid currentTime - requires valid pathInput setup
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))
  saveRDS("dummy", file.path(tmp, "x_001.rds"))
  
  expect_error(
    createModObject(data = "x", pathInput = tmp, currentTime = "one"),
    "Current time needs to be numeric!"
  )
  expect_error(
    createModObject(data = "x", pathInput = tmp, currentTime = NA),
    "Current time needs to be numeric!"
  )
})

# 2. simList behavior
test_that("simList handling works", {
  # Object in simList
  sim <- list(foo = "bar")
  expect_silent(
    result <- createModObject("foo", sim = sim, pathInput = "dummy", currentTime = 1)
  )
  expect_equal(result, "bar")
  
  # Prefers sim over disk
  tmp <- tempfile()
  dir.create(tmp)
  saveRDS("disk", file.path(tmp, "foo001.rds"))
  expect_silent(
    result <- createModObject("foo", sim = sim, pathInput = tmp, currentTime = 1)
  )
  expect_equal(result, "bar")
  
  # Uses sim even with invalid pathInput
  result <- createModObject("foo", sim = sim, pathInput = "invalid_path", currentTime = 1)
  expect_equal(result, "bar")
})

# 3. Empty directory handling
test_that("Empty directory handling works", {
  tmp <- tempfile()
  dir.create(tmp)
  
  # returnNULL = FALSE
  expect_error(
    createModObject("foo", pathInput = tmp, currentTime = 1),
    "Please place the data in the input folder"
  )
  
  # returnNULL = TRUE
  expect_message(
    result <- createModObject("foo", pathInput = tmp, currentTime = 1, returnNULL = TRUE),
    "The file for foo was not found. Returning NULL"
  )
  expect_null(result)
})

# 4. File loading logic
test_that("File loading works correctly", {
  tmp <- tempfile()
  dir.create(tmp)
  
  # Successful load
  obj <- list(a = 1)
  saveRDS(obj, file.path(tmp, "foo001.rds"))
  expect_message(
    result <- createModObject("foo", pathInput = tmp, currentTime = 1),
    "foo loaded from .* for time 001"
  )
  expect_equal(result, obj)
  
  # Custom loader
  writeLines("text", file.path(tmp, "foo002.txt"))
  result <- createModObject("foo", pathInput = tmp, currentTime = 2, fun = readLines)
  expect_equal(result, "text")
  
  # Recursive search
  dir.create(file.path(tmp, "nested"))
  saveRDS("nested", file.path(tmp, "nested/foo003.rds"))
  expect_equal(
    createModObject("foo", pathInput = tmp, currentTime = 3),
    "nested"
  )
  
  # Fallback message
  expect_message(
    createModObject("foo", sim = list(), pathInput = tmp, currentTime = 1),
    "not supplied by another module"
  )
})

# 5. File matching logic
test_that("File matching works correctly", {
  tmp <- tempfile()
  dir.create(tmp)
  
  # Time formatting
  saveRDS("padded", file.path(tmp, "time_001.rds"))
  saveRDS("unpadded", file.path(tmp, "time_1.rds"))
  expect_equal(
    createModObject("time", pathInput = tmp, currentTime = 1),
    "padded"
  )
  
  # Large times
  saveRDS("big", file.path(tmp, "big_1000.rds"))
  expect_equal(
    createModObject("big", pathInput = tmp, currentTime = 1000),
    "big"
  )
  
  # Partial name matching
  saveRDS("partial", file.path(tmp, "partial_001.rds"))
  expect_equal(
    createModObject("part", pathInput = tmp, currentTime = 1),
    "partial"
  )
  
  # Case sensitivity
  saveRDS("upper", file.path(tmp, "UPPER_001.rds"))
  expect_null(
    suppressMessages(
      createModObject("upper", pathInput = tmp, currentTime = 1, returnNULL = TRUE)
    )
  )
})

# 6. Multi-file handling
test_that("Multi-file selection works", {
  tmp <- tempfile()
  dir.create(tmp)
  
  # Same directory
  saveRDS("old", file.path(tmp, "multi_001.rds"))
  Sys.sleep(0.1)
  saveRDS("new", file.path(tmp, "multi_001_v2.rds"))
  
  expect_message(
    result <- createModObject("multi", pathInput = tmp, currentTime = 1),
    "Choosing most recent: multi_001_v2.rds"
  )
  expect_equal(result, "new")
  
  # Different directories
  dir.create(file.path(tmp, "dir1"))
  dir.create(file.path(tmp, "dir2"))
  saveRDS("dir1", file.path(tmp, "dir1/multi_001.rds"))
  saveRDS("dir2", file.path(tmp, "dir2/multi_001.rds"))
  Sys.setFileTime(file.path(tmp, "dir1/multi_001.rds"), Sys.time() - 1000)
  Sys.setFileTime(file.path(tmp, "dir2/multi_001.rds"), Sys.time())
  
  expect_message(
    result <- createModObject("multi", pathInput = tmp, currentTime = 1),
    "Choosing most recent: dir2/multi_001.rds"
  )
  expect_equal(result, "dir2")
})

# 7. Edge cases and failures
test_that("Edge cases and failures are handled", {
  tmp <- tempfile()
  dir.create(tmp)
  
  # No matching files
  saveRDS("data", file.path(tmp, "nomatch_002.rds"))
  expect_message(
    result <- createModObject("test", pathInput = tmp, currentTime = 1),
    "No file found for 1. Returning NULL"
  )
  expect_null(result)
  
  # Loader failure
  saveRDS("ok", file.path(tmp, "err_001.rds"))
  bad_loader <- function(path) stop("Loader error")
  expect_message(
    result <- createModObject("err", pathInput = tmp, currentTime = 1, fun = bad_loader),
    "Failed to load err_001.rds. Returning NULL"
  )
  expect_null(result)
  
  # Corrupted file
  writeBin(charToRaw("invalid"), file.path(tmp, "corrupt_001.rds"))
  expect_message(
    result <- createModObject("corrupt", pathInput = tmp, currentTime = 1),
    "Failed to load corrupt_001.rds. Returning NULL"
  )
  expect_null(result)
  
  # Hidden files
  saveRDS("visible", file.path(tmp, "visible_001.rds"))
  saveRDS("hidden", file.path(tmp, ".hidden_001.rds"))
  expect_equal(
    createModObject("visible", pathInput = tmp, currentTime = 1),
    "visible"
  )
  expect_null(
    suppressMessages(
      createModObject(".hidden", pathInput = tmp, currentTime = 1, returnNULL = TRUE)
    )
  )
})