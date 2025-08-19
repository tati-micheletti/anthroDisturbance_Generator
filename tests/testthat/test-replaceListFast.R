library(testthat)
library(terra)

# build inputs
distList   <- createDisturbanceList()
distParam  <- createDisturbanceParameters(distList)
updAll     <- createUpdatedLayersAll()

# run function under test
out <- replaceListFast(
  disturbanceList       = distList,
  updatedLayersAll      = updAll,
  currentTime           = 42,
  disturbanceParameters = distParam
)

test_that("no updates returns NULL for pipelines", {
  expect_null(out$pipelines)
})

test_that("createdInSimulationTime is stamped on new SpatVector for settlements", {
  stamped <- out$settlements[["settlements"]]
  expect_s4_class(stamped, "SpatVector")
  expect_true("createdInSimulationTime" %in% names(stamped))
  expect_equal(unique(stamped$createdInSimulationTime), 42)
})

test_that("Enlarging type replaces old layer entirely for wind", {
  # old was empty, so we should see exactly the new feature
  new_extent <- ext(updAll$individuaLayers$wind$windTurbines)
  result     <- out$wind[["windTurbines"]]
  expect_s4_class(result, "SpatVector")
  expect_equal(terra::xmin(ext(result)), terra::xmin(new_extent))
  expect_equal(terra::xmax(ext(result)), terra::xmax(new_extent))
  expect_equal(terra::ymin(ext(result)), terra::ymin(new_extent))
  expect_equal(terra::ymax(ext(result)), terra::ymax(new_extent))
  # and no historical geometry inside:
  expect_equal(nrow(result), 1)
})

test_that("Merging raster & vector yields a SpatVector for mining", {
  merged <- out$mining[["mining"]]
  expect_s4_class(merged, "SpatVector")
  # should have more than the single original polygon
  expect_gt(nrow(merged), 1)
})

test_that("Geometry mismatch triggers buffering and merges both features for forestry", {
  merged <- out$forestry[["cutblocks"]]
  expect_s4_class(merged, "SpatVector")
  # original was 1 polygon, new was 1 line → after buffer+merge should be 2+
  expect_gt(nrow(merged), 1)
})

test_that("Special seismicLinesFirstYear branch merges correctly for oilGas seismicLines", {
  merged <- out$oilGas[["seismicLines"]]
  # should combine first_year + current → 2 features
  expect_s4_class(merged, "SpatVector")
  expect_equal(nrow(merged), 2)
})

test_that("Potential layers are preserved alongside updates", {
  # settlements
  names_set <- names(out$settlements)
  expect_true(all(c("potentialSettlements","settlements") %in% names_set))
  # oilGas
  names_oil <- names(out$oilGas)
  expect_true(all(c("potentialSeismicLines","seismicLines") %in% names_oil))
})

test_that("existing createdInSimulationTime in currDist is not overwritten", {
  # Prepare an updated layer with its own timestamp
  curr <- updAll$individuaLayers$settlements$settlements
  curr$createdInSimulationTime <- 99
  upd2 <- updAll
  upd2$individuaLayers$settlements$settlements <- curr
  
  out2 <- replaceListFast(
    disturbanceList       = distList,
    updatedLayersAll      = upd2,
    currentTime           = 42,
    disturbanceParameters = distParam
  )
  
  # Because 'settlements' is Enlarging, replaceListFast() returns curr unchanged:
  res <- out2$settlements[["settlements"]]
  expect_equal(unique(res$createdInSimulationTime), 99)
})

test_that("raster + raster merging yields SpatVector", {
  # Force both past and curr mining layers to be SpatRaster
  distList2 <- distList
  distList2$mining$mining <- updAll$individuaLayers$mining$mining
  
  out3 <- replaceListFast(
    disturbanceList       = distList2,
    updatedLayersAll      = updAll,
    currentTime           = 42,
    disturbanceParameters = distParam
  )
  merged <- out3$mining[["mining"]]
  
  expect_s4_class(merged, "SpatVector")
  expect_gt(nrow(merged), 1)
})

test_that("Enlarging with non-empty disturbanceEnd merges rather than replaces (pipelines$roads)", {
  # Create a small new roads layer
  new_roads <- vect(
    st_sfc(st_linestring(rbind(c(1,5), c(2,5))), crs = crs(r))
  )
  upd3 <- updAll
  upd3$individuaLayers$pipelines <- list(roads = new_roads)
  
  # Modify parameters so 'roads' is marked Enlarging *but* disturbanceEnd != ""
  distParam3 <- copy(distParam)
  distParam3[
    dataName == "pipelines" & disturbanceOrigin == "roads",
    `:=`(disturbanceType="Enlarging", disturbanceEnd="roads")
  ]
  
  out4 <- replaceListFast(
    disturbanceList       = distList,
    updatedLayersAll      = upd3,
    currentTime           = 42,
    disturbanceParameters = distParam3
  )
  merged_roads <- out4$pipelines[["roads"]]
  
  expect_s4_class(merged_roads, "SpatVector")
  # original had 1 feature, new 1 feature → should see 2
  expect_equal(nrow(merged_roads), 2)
})

test_that("merged mining vector has a 'Class' column set to the layer name", {
  merged <- out$mining[["mining"]]
  expect_true("Class" %in% names(merged))
  # All entries should have Class == "mining"
  expect_true(all(merged$Class == "mining"))
})

test_that("pastDist is stamped when no seismicLinesFirstYear provided", {
  upd4 <- updAll
  upd4$seismicLinesFirstYear <- NULL
  
  out5 <- replaceListFast(
    disturbanceList       = distList,
    updatedLayersAll      = upd4,
    currentTime           = 42,
    disturbanceParameters = distParam
  )
  merged <- out5$oilGas[["seismicLines"]]
  
  expect_true("createdInSimulationTime" %in% names(merged))
  # one feature from past (42-10 = 32), one from current (42)
  expect_setequal(merged$createdInSimulationTime, c(32, 42))
})

test_that("potential layers always come first in each sector", {
  nm <- names(out$settlements)
  expect_equal(nm[1], "potentialSettlements")
  expect_equal(nm[2], "settlements")
})

test_that("missing updates returns NULL for forestry when no updatedLayers provided", {
  upd2 <- updAll
  upd2$individuaLayers$forestry <- list()  # simulate no new forestry updates
  
  out2 <- replaceListFast(
    disturbanceList       = distList,
    updatedLayersAll      = upd2,
    currentTime           = 42,
    disturbanceParameters = distParam
  )
  
  expect_null(out2$forestry)
})

test_that("same-class, same-geometry merge for mining polygons (no buffering)", {
  # override mining update to a polygon so both past & curr are polygons
  new_poly2 <- vect(
    st_sfc(st_polygon(list(rbind(
      c(5.5,1), c(5.5,2), c(6,2), c(6,1), c(5.5,1)
    ))), crs = crs(r))
  )
  upd3 <- updAll
  upd3$individuaLayers$mining$mining <- new_poly2
  
  # ensure mining remains “Generating” so it merges rather than replaces
  distParam3 <- copy(distParam)
  distParam3[
    dataName=="mining" & disturbanceOrigin=="mining",
    disturbanceType := "Generating"
  ]
  
  out3 <- replaceListFast(
    disturbanceList       = distList,
    updatedLayersAll      = upd3,
    currentTime           = 42,
    disturbanceParameters = distParam3
  )
  
  merged <- out3$mining[["mining"]]
  expect_s4_class(merged, "SpatVector")
  # exactly two polygons (past + curr), no buffering has occurred
  expect_equal(nrow(merged), 2)
})

test_that("unsupported class combination throws a clear error for a non-Enlarging sector", {
  upd4 <- updAll
  # pipelines is GENERATING by default
  upd4$individuaLayers$pipelines <- list(
    roads = list(foo = "bar")
  )
  
  expect_error(
    replaceListFast(
      disturbanceList       = distList,
      updatedLayersAll      = upd4,
      currentTime           = 42,
      disturbanceParameters = distParam
    ),
    "currDist for 'pipelines::roads.foo' is invalid class: character"
  )
})

test_that("multi-layer updates with mixed types in pipelines sector", {
  # prepare new pipelines and roads updates
  new_pipe <- vect(st_sfc(st_point(c(2,2))))
  crs(new_pipe) <- crs(r)
  new_roads <- vect(st_sfc(st_linestring(rbind(c(1,5), c(2,5)))))
  crs(new_roads) <- crs(r)
  
  upd2 <- updAll
  upd2$individuaLayers$pipelines <- list(
    pipelines = new_pipe,
    roads     = new_roads
  )
  
  distParam2 <- copy(distParam)
  # make roads Enlarging (replace) and leave pipelines Generating (merge)
  distParam2[
    dataName == "pipelines" & disturbanceOrigin == "roads",
    `:=`(disturbanceType = "Enlarging", disturbanceEnd = "")
  ]
  
  out2 <- replaceListFast(distList, upd2, currentTime = 42, disturbanceParameters = distParam2)
  
  # both updated layers should appear
  expect_true(all(c("pipelines", "roads") %in% names(out2$pipelines)))
  
  # pipelines (Generating) merges: empty past + new = 1 feature
  expect_s4_class(out2$pipelines$pipelines, "SpatVector")
  expect_equal(nrow(out2$pipelines$pipelines), 1)
  
  # roads (Enlarging & disturbanceEnd == "") replaces: only new feature
  expect_s4_class(out2$pipelines$roads, "SpatVector")
  expect_equal(nrow(out2$pipelines$roads), 1)
})


test_that("Enlarging with disturbanceEnd != '' merges rather than replaces for settlements", {
  new_set <- vect(st_sfc(st_polygon(list(rbind(
    c(2,2), c(2,4), c(4,4), c(4,2), c(2,2)
  )))))
  crs(new_set) <- crs(r)
  
  upd2 <- updAll
  upd2$individuaLayers$settlements <- list(settlements = new_set)
  
  distParam2 <- copy(distParam)
  # make settlements Enlarging but with non-empty end -> merge
  distParam2[
    dataName == "settlements" & disturbanceOrigin == "settlements",
    `:=`(disturbanceType = "Enlarging", disturbanceEnd = "settlements")
  ]
  
  out2 <- replaceListFast(distList, upd2, currentTime = 42, disturbanceParameters = distParam2)
  merged <- out2$settlements$settlements
  
  expect_s4_class(merged, "SpatVector")
  # should combine past + current => 2 features
  expect_equal(nrow(merged), 2)
})


test_that("NULL pastDist with non-null currDist returns currDist for forestry", {
  distList2 <- distList
  distList2$forestry$cutblocks <- NULL
  
  new_fb <- vect(st_sfc(st_polygon(list(rbind(
    c(8,4), c(8,6), c(9,6), c(9,4), c(8,4)
  )))))
  crs(new_fb) <- crs(r)
  
  upd2 <- updAll
  upd2$individuaLayers$forestry <- list(cutblocks = new_fb)
  
  out2 <- replaceListFast(distList2, upd2, currentTime = 42, disturbanceParameters = distParam)
  result <- out2$forestry$cutblocks
  
  expect_s4_class(result, "SpatVector")
  expect_true("Class" %in% names(result))
  expect_equal(nrow(result), 1)
})

test_that("when pastDist is NULL and currDist is a raster, returns raster and no 'Class'", {
  skip_on_cran()

  r <- rast(nrows = 2, ncols = 2, xmin = 0, xmax = 10, ymin = 0, ymax = 10, crs = "EPSG:3857")
  values(r) <- 1
  
  # No past disturbance; only a new raster layer for 'mining::mining'
  disturbanceList <- list(mining = list(mining = NULL))
  updatedLayersAll <- list(individuaLayers = list(mining = list(mining = r)))
  disturbanceParameters <- data.table(
    dataName = "mining", disturbanceOrigin = "mining",
    disturbanceEnd = "", disturbanceType = "Generating"
  )
  
  out <- replaceListFast(
    disturbanceList = disturbanceList,
    updatedLayersAll = updatedLayersAll,
    currentTime = 2015L,
    disturbanceParameters = disturbanceParameters
  )
  
  expect_true(inherits(out$mining[["mining"]], "SpatRaster"))
  expect_false("Class" %in% names(out$mining[["mining"]]))  # no attempt to add an attribute to a raster
})


test_that("merged vector overwrites a numeric 'Class' with the layer name as character", {
  skip_on_cran()

  crs_str <- "EPSG:3857"
  # past vector (with numeric Class)
  past <- vect(matrix(c(0,0, 4,0, 4,4, 0,4, 0,0), ncol = 2, byrow = TRUE), type = "polygons", crs = crs_str)
  past$Class <- 1L  # numeric -> this caused NA when assigning a character later
  
  # current vector (no Class)
  curr <- vect(matrix(c(6,0, 10,0, 10,4, 6,4, 6,0), ncol = 2, byrow = TRUE), type = "polygons", crs = crs_str)
  
  disturbanceList   <- list(mining = list(mining = past))
  updatedLayersAll  <- list(individuaLayers = list(mining = list(mining = curr)))
  disturbanceParameters <- data.table(
    dataName = "mining", disturbanceOrigin = "mining",
    disturbanceEnd = "", disturbanceType = "Generating"
  )
  
  out <- replaceListFast(
    disturbanceList = disturbanceList,
    updatedLayersAll = updatedLayersAll,
    currentTime = 2015L,
    disturbanceParameters = disturbanceParameters
  )
  
  merged <- out$mining[["mining"]]
  
  expect_s4_class(merged, "SpatVector")
  expect_true("Class" %in% names(merged))
  expect_true(is.character(merged$Class))
  expect_true(all(merged$Class == "mining"))
})

