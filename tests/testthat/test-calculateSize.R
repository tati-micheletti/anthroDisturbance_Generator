library(testthat)
library(data.table)
library(terra)
library(crayon)
library(sf)
library(msm)

# Create helper inputs
r <- rast(nrows=10, ncols=10, xmin=0, xmax=10, ymin=0, ymax=10, vals=1)
crs(r) <- "EPSG:3005"
# Utilize provided helper functions
disturbanceList <- createDisturbanceList(crs="EPSG:3005")
disturbanceParameters <- createDisturbanceParameters(disturbanceList)

# Adjusted test for wind turbines with current behavior
test_that("potentialWindTurbines missing layer results in row removal", {
  dp <- data.table(
    dataName = "energy",
    dataClass = "potentialWindTurbines",
    disturbanceOrigin = "wt",
    disturbanceType = "Generating"
  )
  dl <- list(energy = list(wt = NULL))
  out <- calculateSize(disturbanceParameters = dp, disturbanceList = dl, whichToUpdate = 1)
  expect_equal(nrow(out), 0)  # Expecting row to be removed
})

# Updated test to check if missing non-wind layers result in row removal
test_that("calculateSize removes rows for non-wind classes with missing layers", {
  disturbanceList$oilGas$missingLayer <- NULL
  disturbanceParameters <- rbind(disturbanceParameters, data.table(
    dataName="oilGas", dataClass="missingLayer", disturbanceType="Generating",
    disturbanceRate=NA_real_, disturbanceSize=NA_real_, disturbanceOrigin="missingLayer",
    disturbanceEnd="", disturbanceInterval=1L, resolutionVector=list(NA)
  ))
  
  updatedParams <- calculateSize(disturbanceParameters,
                                 disturbanceList,
                                 whichToUpdate = which(disturbanceParameters$dataClass == "missingLayer"))
  
  expect_false("missingLayer" %in% updatedParams$dataClass)
})

# Adjusted test remains for raster input correctness
test_that("calculateSize converts raster inputs correctly", {
  disturbanceList$forestry$rasterDisturbance <- r
  disturbanceParameters <- rbind(disturbanceParameters, data.table(
    dataName="forestry", dataClass="rasterDisturbance", disturbanceType="Generating",
    disturbanceRate=NA_real_, disturbanceSize=NA_real_, disturbanceOrigin="rasterDisturbance",
    disturbanceEnd="", disturbanceInterval=1L, resolutionVector=list(NA)
  ))
  
  updatedParams <- calculateSize(disturbanceParameters,
                                 disturbanceList,
                                 whichToUpdate = which(disturbanceParameters$dataClass == "rasterDisturbance"))
  
  expect_true(!is.na(updatedParams[dataClass == "rasterDisturbance", disturbanceSize]))
  expect_true(grepl("rtnorm", updatedParams[dataClass == "rasterDisturbance", disturbanceSize]))
})

# Remaining tests (unchanged)
test_that("calculateSize computes disturbanceSize correctly for polygon inputs", {
  disturbanceParameters[dataClass == "potentialSettlements", disturbanceSize := NA]
  updatedParams <- calculateSize(disturbanceParameters,
                                 disturbanceList,
                                 whichToUpdate = which(disturbanceParameters$dataClass == "potentialSettlements"))
  
  expect_true(!is.na(updatedParams[dataClass == "potentialSettlements", disturbanceSize]))
  expect_true(grepl("rtnorm", updatedParams[dataClass == "potentialSettlements", disturbanceSize]))
})

test_that("calculateSize computes disturbanceSize correctly for line inputs", {
  disturbanceParameters[dataClass == "potentialSeismicLines", disturbanceSize := NA]
  updatedParams <- calculateSize(disturbanceParameters,
                                 disturbanceList,
                                 whichToUpdate = which(disturbanceParameters$dataClass == "potentialSeismicLines"))
  
  expect_true(!is.na(updatedParams[dataClass == "potentialSeismicLines", disturbanceSize]))
  expect_true(grepl("rtnorm", updatedParams[dataClass == "potentialSeismicLines", disturbanceSize]))
})

test_that("calculateSize calculates correct values for numeric stability", {
  tiny_poly <- st_polygon(list(rbind(c(0.01,0.01),c(0.01,0.02),c(0.02,0.02),c(0.02,0.01),c(0.01,0.01))))
  disturbanceList$testSector <- list(tinyDisturbance = to_sv(st_sfc(tiny_poly), "tinyDisturbance", "EPSG:3005"))
  disturbanceParameters <- rbind(disturbanceParameters, data.table(
    dataName="testSector", dataClass="tinyDisturbance", disturbanceType="Generating",
    disturbanceRate=NA_real_, disturbanceSize=NA_real_, disturbanceOrigin="tinyDisturbance",
    disturbanceEnd="", disturbanceInterval=1L, resolutionVector=list(NA)
  ))
  
  updatedParams <- calculateSize(disturbanceParameters,
                                 disturbanceList,
                                 whichToUpdate = which(disturbanceParameters$dataClass == "tinyDisturbance"))
  
  expect_true(!is.na(updatedParams[dataClass == "tinyDisturbance", disturbanceSize]))
  expect_true(grepl("rtnorm", updatedParams[dataClass == "tinyDisturbance", disturbanceSize]))
})

test_that("calculateSize handles zero-variance inputs with positive sigma", {
  skip_on_cran()
  
  # One single polygon -> sd(area) would be NA/0 without our fallback
  single_poly <- sf::st_polygon(list(rbind(c(0,0), c(0,1), c(1,1), c(1,0), c(0,0))))
  sfc <- sf::st_sfc(single_poly, crs = 3005)
  disturbanceList$testSector$singlePolyDisturbance <- to_sv(sfc, "singlePolyDisturbance", "EPSG:3005")
  
  disturbanceParameters <- rbind(
    disturbanceParameters,
    data.table::data.table(
      dataName="testSector", dataClass="singlePolyDisturbance", disturbanceType="Generating",
      disturbanceRate=NA_real_, disturbanceSize=NA_character_, disturbanceOrigin="singlePolyDisturbance",
      disturbanceEnd="", disturbanceInterval=1L, resolutionVector=list(NA)
    ),
    fill = TRUE
  )
  
  idx <- which(disturbanceParameters$dataClass == "singlePolyDisturbance")
  updatedParams <- calculateSize(disturbanceParameters, disturbanceList, whichToUpdate = idx)
  
  s <- updatedParams[dataClass == "singlePolyDisturbance", disturbanceSize]
  expect_true(!is.na(s))
  
  # Expect format: rtnorm(1, <mu>, <sigma>, lower = 0) with sigma > 0
  m <- regexec("^rtnorm\\(1,\\s*([0-9.]+),\\s*([0-9.]+),\\s*lower\\s*=\\s*0\\)$", s)
  caps <- regmatches(s, m)[[1]]
  expect_gt(length(caps), 2)                               # captured mu and sigma
  mu    <- as.numeric(caps[2])
  sigma <- as.numeric(caps[3])
  
  expect_true(is.finite(mu)    && mu > 0)
  expect_true(is.finite(sigma) && sigma > 0) 
  
  # Also ensure it's parsable/evaluable and non-negative
  val <- eval(parse(text = s))
  expect_true(is.numeric(val) && length(val) == 1L && is.finite(val) && val >= 0)
})

