# Unit tests for extractNonPotentialLayers()
# Helper to compare data.tables ignoring row order
expect_dt_equal <- function(actual, expected) {
  setkey(actual, Sector, dataClass)
  setkey(expected, Sector, dataClass)
  expect_equal(actual, expected)
}

# 1. Empty disturbanceList returns empty data.table
test_that("empty list returns empty data.table", {
  empty <- list()
  result <- extractNonPotentialLayers(empty)
  expect_s3_class(result, "data.table")
  expect_true(nrow(result) == 0)
})

# 2. Only potential layers should be filtered out
test_that("only potential layers yields empty result", {
  dl <- list(
    forestry = list(
      potentialCutblocks = NULL,
      potentialRoads     = NULL
    )
  )
  result <- extractNonPotentialLayers(dl)
  expect_s3_class(result, "data.table")
  expect_true(nrow(result) == 0)
})

# 3. Single sector with mixed layers
test_that("mixed layers returns only non-potential ones", {
  dl <- list(
    energy = list(
      potentialWind = NULL,
      solar         = NULL,
      gridLines     = NULL
    )
  )
  result <- extractNonPotentialLayers(dl)
  expected <- data.table(
    Sector    = c("energy", "energy"),
    dataClass = c("solar", "gridLines")
  )
  expect_dt_equal(result, expected)
})

# 4. Multiple sectors and duplicates
test_that("multiple sectors and duplicate entries are handled correctly", {
  dl <- list(
    mining = list(
      minePits        = NULL,
      potentialShafts = NULL,
      minePits        = NULL  # duplicate name
    ),
    forestry = list(
      cutblocks       = NULL,
      thinning        = NULL,
      potentialFire   = NULL
    )
  )
  result <- extractNonPotentialLayers(dl)
  expected <- data.table(
    Sector    = c("mining", "forestry", "forestry"),
    dataClass = c("minePits", "cutblocks", "thinning")
  )
  expect_dt_equal(result, expected)
})

# 5. Names that contain 'potential' but not at start are kept
test_that("data classes containing 'potential' not at start are kept", {
  dl <- list(
    sectorX = list(
      subpotentialLayer = NULL,
      potentialLayer    = NULL,
      actualLayer       = NULL
    )
  )
  result <- extractNonPotentialLayers(dl)
  expected <- data.table(
    Sector    = c("sectorX", "sectorX"),
    dataClass = c("subpotentialLayer", "actualLayer")
  )
  expect_dt_equal(result, expected)
})

# Edge-case tests

# 6. Empty second-level lists return empty result


# 7. Unnamed sectors yield empty result
test_that("unnamed sectors yield empty result", {
  dl <- list(
    list(layer1 = NULL, layer2 = NULL)
  )
  names(dl) <- NULL
  result <- extractNonPotentialLayers(dl)
  expect_s3_class(result, "data.table")
  expect_true(nrow(result) == 0)
})

# 8. Classes with uppercase 'Potential' prefix are kept (case-sensitive filter)
test_that("uppercase 'Potential' prefix is not filtered out", {
  dl <- list(
    sectorB = list(
      PotentialLayer = NULL,
      potentialLayer = NULL,
      realLayer      = NULL
    )
  )
  result <- extractNonPotentialLayers(dl)
  expected <- data.table(
    Sector    = c("sectorB", "sectorB"),
    dataClass = c("PotentialLayer", "realLayer")
  )
  expect_dt_equal(result, expected)
})

# 9. Empty dataClass names are retained
# Construct second-level list with an empty name using setNames

test_that("empty dataClass names are retained", {
  second_level <- setNames(list(NULL, NULL, NULL), c("", "potentialXYZ", "actualXYZ"))
  dl <- list(sectorC = second_level)
  result <- extractNonPotentialLayers(dl)
  expected <- data.table(
    Sector    = c("sectorC", "sectorC"),
    dataClass = c("", "actualXYZ")
  )
  expect_dt_equal(result, expected)
})

# 10. Atomic (non-list) sector values produce empty result
test_that("atomic sector values produce empty data.table", {
  dl <- list(
    sectorD = 42
  )
  result <- extractNonPotentialLayers(dl)
  expect_s3_class(result, "data.table")
  expect_true(nrow(result) == 0)
})

# 11. NA dataClass names cause an error
test_that("NA dataClass names are filtered out", {
  second_level <- setNames(list(NULL, NULL), c(NA_character_, "foo"))
  dl <- list(x = second_level)
  result <- extractNonPotentialLayers(dl)
  expected <- data.table(
    Sector    = "x",
    dataClass = "foo"
  )
  expect_dt_equal(result, expected)
})

# 12. 
test_that("dataClass named 'NA' is kept", {
  dl <- list(sector = list(`NA` = NULL, potentialX = NULL))
  result <- extractNonPotentialLayers(dl)
  expected <- data.table(Sector = "sector", dataClass = "NA")
  expect_equal(result, expected)
})

# 13. Named vector sector processes names
test_that("named vector sector processes names", {
  dl <- list(sector = c(a = 1, b = 2, potentialX = 3))
  result <- extractNonPotentialLayers(dl)
  # Named vector is not a list, so should be skipped entirely\ n  expect_s3_class(result, "data.table")
  expect_true(nrow(result) == 0)
})

# 15
test_that("duplicate dataClass names in different sectors are kept", {
  dl <- list(
    sector1 = list(x = NULL),
    sector2 = list(x = NULL)
  )
  result <- extractNonPotentialLayers(dl)
  expected <- data.table(Sector = c("sector1", "sector2"), dataClass = c("x", "x"))
  expect_dt_equal(result, expected)
})

# 16
test_that("whitespace after 'potential' is still filtered", {
  dl <- list(sector = list("potential " = NULL, real = NULL))
  result <- extractNonPotentialLayers(dl)
  expected <- data.table(Sector = "sector", dataClass = "real")
  expect_dt_equal(result, expected)
})

# 17
test_that("nested lists are handled correctly", {
  dl <- list(
    sectorA = list(
      subSector = list(
        layer1 = NULL,
        potentialLayer = NULL
      )
    )
  )
  result <- extractNonPotentialLayers(dl)
  expected <- data.table(
    Sector = "sectorA",
    dataClass = "subSector"
  )
  expect_dt_equal(result, expected)
})

# 18
test_that("special characters in dataClass names are handled correctly", {
  dl <- list(
    sectorA = list(
      "layer@1" = NULL,
      "potential#Layer" = NULL
    )
  )
  result <- extractNonPotentialLayers(dl)
  expected <- data.table(
    Sector = "sectorA",
    dataClass = "layer@1"
  )
  expect_dt_equal(result, expected)
})

# 19
test_that("whitespace in dataClass names is handled correctly", {
  dl <- list(
    sectorA = list(
      " layer1 " = NULL,
      "potentialLayer " = NULL
    )
  )
  result <- extractNonPotentialLayers(dl)
  expected <- data.table(
    Sector = "sectorA",
    dataClass = " layer1 "
  )
  expect_dt_equal(result, expected)
})
