# tests/testthat/test-saveDisturbances.R
# Set up temporary output directory and clear previous disturbance files
output_dir <- tempdir()
Paths <<- list(outputPath = output_dir) # Use <<- to assign in the global environment
# Clean up any existing disturbance files
existing <- list.files(output_dir, pattern = "^disturbances_", full.names = TRUE)
if (length(existing)) unlink(existing, recursive = TRUE)

# Core functionality tests

test_that("single SpatVector sector is saved as shapefile", {
  pts <- vect(sf::st_as_sf(data.frame(x = 1:2, y = 3:4), coords = c("x", "y")))
  disturbances <- list(sectorA = pts)
  expect_message(
    saveDisturbances(disturbances,
                     currentTime = "T1",
                     overwrite = TRUE,
                     runName = "runX"),
    "Layer: sectorA was likely not generated"
  )
  base <- file.path(output_dir, "disturbances_sectorA_T1_runX")
  expect_true(file.exists(paste0(base, ".shp")))
  expect_true(file.exists(paste0(base, ".shx")))
  expect_true(file.exists(paste0(base, ".dbf")))
})

test_that("single RasterLayer sector is saved as tiff", {
  r <- raster::raster(nrows = 5, ncols = 5)
  raster::values(r) <- 1
  disturbances <- list(sectorB = r)
  expect_message(
    saveDisturbances(disturbances,
                     currentTime = "IC",
                     overwrite = TRUE,
                     runName = "runY"),
    "Layer: sectorB was likely not generated"
  )
  tif <- file.path(output_dir, "disturbances_sectorB_IC_runY.tif")
  expect_true(file.exists(tif))
})

test_that("multi-layer sector saves only non-potential layers", {
  vec <- vect(sf::st_as_sf(data.frame(x = 10, y = 20), coords = c("x", "y")))
  rastL <- rast(nrows = 3, ncols = 3)
  values(rastL) <- 2
  disturbances <- list(sectorC = list(
    potentialLayer = vec,
    actual = vec,
    rastL = rastL
  ))
  expect_message(
    saveDisturbances(disturbances,
                     currentTime = 2025,
                     overwrite = FALSE,
                     runName = "testRun"),
    "Saving layer: sectorC -- actual"
  )
  shp_actual <- file.path(output_dir, "disturbances_sectorC_actual_2025_testRun.shp")
  tif_rastL <- file.path(output_dir, "disturbances_sectorC_rastL_2025_testRun.tif")
  expect_true(file.exists(shp_actual))
  expect_true(file.exists(tif_rastL))
  # potential layer should not be saved
  no_shp <- file.path(output_dir, "disturbances_sectorC_potentialLayer_2025_testRun.shp")
  expect_false(file.exists(no_shp))
})

test_that("NULL layers issue warnings and are not saved", {
  disturbances <- list(sectorD = list(layer1 = NULL))
  expect_warning(
    saveDisturbances(disturbances,
                     currentTime = "T2",
                     overwrite = TRUE,
                     runName = "r1"),
    "The layer for sectorD -- layer1 is NULL. Not saving."
  )
  file_path <- file.path(output_dir, "disturbances_sectorD_layer1_T2_r1.shp")
  expect_false(file.exists(file_path))
})

test_that("unsupported classes throw errors", {
  disturbances <- list(sectorE = list(layer1 = data.frame(a = 1)))
  expect_error(
    saveDisturbances(disturbances,
                     currentTime = "T3",
                     overwrite = TRUE,
                     runName = "r2"),
    "Objects of class"
  )
})

# Edge case tests

# 1. sf polygon input

test_that("sf polygon is saved as shapefile", {
  poly <- st_as_sf(
    data.frame(id = 1),
    geometry = st_sfc(st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,0)))))
  )
  disturbances <- list(sectorF = poly)
  expect_message(
    saveDisturbances(disturbances,
                     currentTime = "T6",
                     overwrite = TRUE,
                     runName = "runZ"),
    "Layer: sectorF was likely not generated"
  )
  shp <- file.path(output_dir, "disturbances_sectorF_T6_runZ.shp")
  expect_true(file.exists(shp))
})

# 2. SpatRaster input

test_that("SpatRaster sector is saved as tiff", {
  sr <- rast(nrows = 2, ncols = 2)
  values(sr) <- 5
  disturbances <- list(sectorG = sr)
  expect_message(
    saveDisturbances(disturbances,
                     currentTime = 2025,
                     overwrite = TRUE,
                     runName = "rg"),
    "Layer: sectorG was likely not generated"
  )
  tif <- file.path(output_dir, "disturbances_sectorG_2025_rg.tif")
  expect_true(file.exists(tif))
})

# 3. sp SpatialPointsDataFrame input

test_that("sp SpatialPointsDataFrame is saved as shapefile", {
  coords <- matrix(c(1,2,3,4), ncol = 2, byrow = TRUE)
  sp_pts <- SpatialPointsDataFrame(coords,
                                   data = data.frame(val = c(1,2)))
  disturbances <- list(sectorH = sp_pts)
  expect_message(
    saveDisturbances(disturbances,
                     currentTime = "IC",
                     overwrite = TRUE,
                     runName = "runSP"),
    "Layer: sectorH was likely not generated"
  )
  shp <- file.path(output_dir, "disturbances_sectorH_IC_runSP.shp")
  expect_true(file.exists(shp))
})

# 4. overwrite = FALSE on existing file triggers error

test_that("overwrite = FALSE with existing file errors", {
  pts <- vect(sf::st_as_sf(data.frame(x = 5, y = 6), coords = c("x", "y")))
  disturbances <- list(sectorI = pts)
  # first write succeeds
  saveDisturbances(disturbances, currentTime = "T4", overwrite = TRUE, runName = "runI")
  # second write should error when overwrite = FALSE
  expect_error(
    saveDisturbances(disturbances, currentTime = "T4", overwrite = FALSE, runName = "runI")
  )
})

# 5. existing Class column is preserved

test_that("existing Class column is preserved", {
  poly2 <- st_as_sf(
    data.frame(Class = c("A", "B"), id = 1:2),
    geometry = st_sfc(st_point(c(0,0)), st_point(c(1,1)))
  )
  disturbances <- list(sectorJ = poly2)
  saveDisturbances(disturbances,
                   currentTime = "T5",
                   overwrite = TRUE,
                   runName = "runJ")
  out <- vect(file.path(output_dir, "disturbances_sectorJ_T5_runJ.shp"))
  expect_equal(as.character(out$Class), c("A", "B"))
})

# 6. sector with only potential layers creates no files

test_that("sector with only potential layers creates no files", {
  vec2 <- vect(sf::st_as_sf(data.frame(x = 2, y = 2), coords = c("x", "y")))
  disturbances <- list(sectorK = list(potential1 = vec2))
  msg <- capture_messages(
    saveDisturbances(disturbances,
                     currentTime = "T6",
                     overwrite = TRUE,
                     runName = "runK")
  )
  expect_true(any(grepl("All disturbances saved for T6", msg)))
  files <- list.files(output_dir, pattern = "sectorK", full.names = TRUE)
  expect_length(files, 0)
})

# 7. case-sensitive potential filtering: uppercase 'PotentialX' is saved

test_that("uppercase 'PotentialX' layer is saved", {
  vec3 <- vect(sf::st_as_sf(data.frame(x = 3, y = 3), coords = c("x", "y")))
  disturbances <- list(sectorL = list(PotentialX = vec3))
  expect_message(
    saveDisturbances(disturbances,
                     currentTime = 2025,
                     overwrite = TRUE,
                     runName = "runL"),
    "Saving layer: sectorL -- PotentialX"
  )
  shp <- file.path(output_dir, "disturbances_sectorL_PotentialX_2025_runL.shp")
  expect_true(file.exists(shp))
})

# 8. empty disturbance list prints final message and writes nothing

test_that("empty disturbance list does nothing but confirms save", {
  msg <- capture_messages(
    saveDisturbances(list(),
                     currentTime = "T7",
                     overwrite = TRUE,
                     runName = "runEmpty")
  )
  expect_true(any(grepl("All disturbances saved for T7", msg)))
  files <- list.files(output_dir, pattern = "runEmpty", full.names = TRUE)
  expect_length(files, 0)
})

# 9. unusual runName and sector names handled correctly

test_that("unusual names with spaces and punctuation are saved", {
  pts <- vect(sf::st_as_sf(data.frame(x = 7, y = 8), coords = c("x", "y")))
  disturbances <- list("sec tor!" = pts)
  run_name <- "run:Special*"
  expect_message(
    saveDisturbances(disturbances,
                     currentTime = "T8",
                     overwrite = TRUE,
                     runName = run_name),
    "Layer: sec tor! was likely not generated"
  )
  # file name will include the raw sector and runName
  pattern <- paste0("disturbances_sec tor!_T8_", run_name)
  files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  expect_true(length(files) >= 1)
})
