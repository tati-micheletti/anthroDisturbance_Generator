# tests/testthat/test-disturbanceList-url.R
# (unwrapTerraList() and its tests moved to the anthroDisturbanceUtils package)

testthat::test_that("the default disturbanceList is the .rds archive on Zenodo, not Google Drive", {
  testthat::skip_on_cran()
  modFile <- testthat::test_path("..", "..", "anthroDisturbance_Generator.R")
  txt <- paste(readLines(modFile), collapse = "\n")
  block <- regmatches(txt, regexpr('(?s)expectsInput\\(objectName = "disturbanceList".*?sourceURL = "[^"]*"',
                                   txt, perl = TRUE))
  url <- sub('.*sourceURL = "([^"]*)"$', "\\1", block)
  testthat::expect_false(grepl("drive.google.com", url, fixed = TRUE))
  testthat::expect_match(url, "^https://zenodo.org/records/[0-9]+/files/disturbanceList.zip$")
})
