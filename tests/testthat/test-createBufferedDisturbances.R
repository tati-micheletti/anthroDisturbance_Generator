library(testthat)
library(terra)
library(sf)
library(raster)

distList <- createDisturbanceList()
studyArea_sv <- terra::vect(st_as_sfc(st_bbox(c(xmin=0,ymin=0,xmax=10,ymax=10), crs=st_crs("EPSG:3005"))))
r_template  <- rast(xmin=0, xmax=10, ymin=0, ymax=10, res=1, crs="EPSG:3005")

# 2. RasterLayer coercion
test_that("coerces RasterLayer to SpatRaster or RasterLayer", {
  r_old <- raster(r_template)
  out <- createBufferedDisturbances(distList, bufferSize=1,
                                    rasterToMatch=r_old,
                                    studyArea=studyArea_sv,
                                    convertToRaster=TRUE)
  expect_true(inherits(out, "SpatRaster"))
})

# 3. Skip non-spatial with warning
test_that("skips non-spatial layers and warns once", {
  badList <- distList; badList$foo <- list(bar = 123)
  expect_warning(
    createBufferedDisturbances(badList, bufferSize=1,
                               rasterToMatch=r_template,
                               studyArea=studyArea_sv,
                               convertToRaster=TRUE),
    "skipping"
  )
})

# 4. Potential-only layers (raster)
test_that("only potential layers return zero-filled raster and skip app()", {
  onlyPot <- list(x = list(potentialA = distList$settlements$potentialSettlements))
  expect_message(
    outR <- createBufferedDisturbances(onlyPot, bufferSize=1,
                                       rasterToMatch=r_template,
                                       studyArea=studyArea_sv,
                                       convertToRaster=TRUE),
    "all-zero raster"
  )
  # if app() is skipped, rasterToMatch[] was used directly
  expect_true(inherits(outR, "SpatRaster"))
  vals <- terra::values(outR)[!is.na(terra::values(outR))]
  expect_true(all(vals == 0))
})

# 5. Potential-only layers (vector)
test_that("only potential layers error early for vector", {
  onlyPot <- list(x = list(potentialA = distList$settlements$potentialSettlements))
  expect_error(
    createBufferedDisturbances(onlyPot, bufferSize=1,
                               rasterToMatch=r_template,
                               studyArea=studyArea_sv,
                               convertToRaster=FALSE),
    "No non-potential disturbance layers"
  )
})

# 6. Core raster output
test_that("raster output is binary mask at key locations", {
  r_template <- terra::rast(r_template)               # ensure SpatRaster
  r_template <- terra::mask(terra::setValues(r_template, 0), studyArea_sv)
  outR <- createBufferedDisturbances(distList, bufferSize=1,
                                     rasterToMatch=r_template,
                                     studyArea=studyArea_sv,
                                     convertToRaster=TRUE)
  expect_true(inherits(outR, "SpatRaster"))
  vals <- terra::values(outR)
  expect_setequal(unique(vals[!is.na(vals)]), c(0,1))
  idx_in  <- cellFromXY(r_template, matrix(c(2,2),1,2))
  idx_out <- cellFromXY(r_template, matrix(c(9,9),1,2))
  expect_equal(vals[idx_in], 1)
  expect_equal(vals[idx_out], 0)
})

# 7. Core vector output
test_that("vector output returns dissolved SpatVector", {
  outV <- createBufferedDisturbances(distList, bufferSize=1,
                                     rasterToMatch=r_template,
                                     studyArea=studyArea_sv,
                                     convertToRaster=FALSE)
  expect_true(inherits(outV, "SpatVector"))
  expect_equal(nrow(outV), 1)
})

# 8. Empty layers
test_that("empty disturbance layers handled", {
  # raster mode
  emptyList <- list(a=list(empty=terra::vect()))
  outR <- createBufferedDisturbances(emptyList, bufferSize=1,
                                     rasterToMatch=r_template,
                                     studyArea=studyArea_sv,
                                     convertToRaster=TRUE)
  vals <- terra::values(outR)[!is.na(terra::values(outR))]
  expect_true(all(vals == 0))
  # vector mode errors
  expect_error(
    createBufferedDisturbances(emptyList, bufferSize=1,
                               rasterToMatch=r_template,
                               studyArea=studyArea_sv,
                               convertToRaster=FALSE)
  )
})

# 9. Zero buffer size behaves as union without re-vect error
test_that("bufferSize=0 does not error on mixed geometry types", {
  # Use expect_no_error instead of expect_silent
  # since the function produces informational messages
  expect_no_error({
    outV <- suppressMessages(
      createBufferedDisturbances(distList, bufferSize=0,
                                 rasterToMatch=r_template,
                                 studyArea=studyArea_sv,
                                 convertToRaster=FALSE)
    )
    # Additional validation
    expect_true(inherits(outV, "SpatVector"))
    expect_true(nrow(outV) > 0)
  })
})

# 11. Overlapping polygons
test_that("overlapping polygons dissolve into one feature", {
  p1 <- st_polygon(list(rbind(c(2,2),c(2,4),c(4,4),c(4,2),c(2,2))))
  p2 <- st_polygon(list(rbind(c(3,3),c(3,5),c(5,5),c(5,3),c(3,3))))
  sv1 <- terra::vect(st_sfc(p1, crs = st_crs("EPSG:3005")))
  sv2 <- terra::vect(st_sfc(p2, crs = st_crs("EPSG:3005")))
  outV <- createBufferedDisturbances(list(x=list(a=sv1,b=sv2)),
                                     bufferSize=0,
                                     rasterToMatch=r_template,
                                     studyArea=studyArea_sv,
                                     convertToRaster=FALSE)
  expect_equal(nrow(outV), 1)
})

# 12. Mixed raster+vector
test_that("mixed raster and vector layers processed together", {
  sv <- terra::vect(st_sfc(st_polygon(list(rbind(c(1,1),c(1,2),c(2,2),c(2,1),c(1,1)))), crs = st_crs("EPSG:3005")))
  r  <- r_template; values(r) <- 0; values(r)[5] <- 1
  dl <- list(m = list(vect=sv, rast=r))
  outV <- createBufferedDisturbances(dl, bufferSize=1,
                                     rasterToMatch=r_template,
                                     studyArea=studyArea_sv,
                                     convertToRaster=FALSE)
  expect_true(inherits(outV, "SpatVector"))
  expect_true(nrow(outV) >= 1)
})

# 13. NA values in raster input are handled by using na.rm in app()
test_that("NA values in raster input produce valid mask", {
  rNA <- r_template; terra::values(rNA) <- NA; terra::values(rNA)[1] <- 1
  outR <- createBufferedDisturbances(list(a=list(layer1=rNA)),
                                     bufferSize=1,
                                     rasterToMatch=r_template,
                                     studyArea=studyArea_sv,
                                     convertToRaster=TRUE)
  expect_true(inherits(outR, "SpatRaster"))
  vals <- terra::values(outR)
  expect_true(all(is.na(vals) | vals %in% c(0,1)))
})

# 15. Study area with holes
test_that("study area with holes yields valid raster mask", {
  outer <- rbind(c(0,0),c(0,10),c(10,10),c(10,0),c(0,0))
  hole  <- rbind(c(2,2),c(2,8),c(8,8),c(8,2),c(2,2))
  sfc_hole <- st_sfc(st_polygon(list(outer, hole)), crs = st_crs("EPSG:3005"))
  sv_hole  <- terra::vect(sfc_hole)
  outR <- createBufferedDisturbances(distList, bufferSize=1,
                                     rasterToMatch=r_template,
                                     studyArea=sv_hole,
                                     convertToRaster=TRUE)
  expect_true(inherits(outR, "SpatRaster"))
  vals <- terra::values(outR)
  expect_true(all(is.na(vals) | vals %in% c(0,1)))
})
