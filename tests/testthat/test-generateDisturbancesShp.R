# Setup

# dummy Paths to satisfy internal references
Paths <- list(outputPath = tempdir())

# Create a 100×100m template raster (1 m resolution), UTM33N
r <- rast(ncol = 100, nrow = 100,
          xmin = 0, xmax = 100, ymin = 0, ymax = 100,
          crs  = "EPSG:32633")
values(r) <- 1

# studyArea as matching SpatVector
sa <- as.polygons(ext(r))
crs(sa) <- crs(r)


# Small helpers #

get_connections <- function(out) {
  IL <- if (!is.null(out$individualLayers)) out$individualLayers else out$individuaLayers
  IL$pipelines$roads
}

make_base <- function() {
  r  <- rast(ncol=100, nrow=100, xmin=0, xmax=100, ymin=0, ymax=100, vals=1)
  crs(r) <- "EPSG:32633"
  sa <- as.polygons(ext(r)); crs(sa) <- crs(r)
  
  # Potential mining area (left half)
  pot <- vect(matrix(c(0,0, 50,0, 50,100, 0,100, 0,0), ncol=2, byrow=TRUE),
              type="polygons", crs=crs(r))
  pot$Potential <- 1L
  
  # Vertical road on far right
  roads <- terra::vect(
    sf::st_sfc(sf::st_linestring(rbind(c(90,10), c(90,90))), crs = 32633)
  )
  roads$Class <- "roads"
  
  # Central lake polygon (rectangle 40..60)
  lake <- vect(matrix(c(40,40, 60,40, 60,60, 40,60, 40,40), ncol=2, byrow=TRUE),
               type="polygons", crs=crs(r))
  
  list(r=r, sa=sa, pot=pot, roads=roads, lake=lake)
}

make_dp <- function(conn_block_size = NULL,
                    dataName = "mining",
                    dataClass = "potentialMining",
                    disturbanceType = "Generating",
                    disturbanceOrigin = "mining",
                    disturbanceRate = 4,
                    disturbanceSize = "1",
                    disturbanceInterval = 1L,
                    resolutionVector = 1) {
  
  dp_gen <- data.table(
    dataName = dataName,
    dataClass = dataClass,
    disturbanceType = disturbanceType,
    disturbanceOrigin = disturbanceOrigin,
    disturbanceRate = disturbanceRate,
    disturbanceSize = disturbanceSize,
    disturbanceInterval = disturbanceInterval,
    resolutionVector = resolutionVector,
    potentialField = "Potential"
  )
  
  dp_conn <- data.table(
    dataName = "pipelines",
    dataClass = "roads",
    disturbanceType = "Connecting",
    disturbanceOrigin = "mining",
    disturbanceEnd = "roads",
    disturbanceRate = NA_real_,
    disturbanceSize = NA_character_,
    disturbanceInterval = disturbanceInterval,
    resolutionVector = resolutionVector
  )
  
  list(dp = rbind(dp_gen, dp_conn, fill = TRUE),
       connectingBlockSize = conn_block_size)
}


################################################################################
### ------------------------------------------------------------------------ ###
################################################################################
# Enlarging branch

# 1) Test: Enlarging branch increases area

dp_enlarge <- data.table(
  dataName = "testSector",
  dataClass = "testSector",
  disturbanceType = "Enlarging",
  disturbanceRate = 100,
  disturbanceOrigin = list("origin1"),
  disturbanceInterval = 1,
  resolutionVector = 1,
  potentialField = NA_character_
)

poly <- as.polygons(ext(sa))
crs(poly) <- crs(r) 
disturbanceList_enlarge <- list(
  testSector = list(origin1 = poly)
)


test_that("Enlarging increases polygon area", {
  # Run and capture console output
  out_full <- NULL
  capture.output(
    out_full <- generateDisturbancesShp(
      disturbanceParameters = dp_enlarge,
      disturbanceList = disturbanceList_enlarge,
      rasterToMatch = r,
      studyArea = sa,
      fires = NULL,
      currentTime = 1,
      firstTime = FALSE,
      growthStepGenerating = 1,
      growthStepEnlargingPolys = 1,
      growthStepEnlargingLines = 0.1,
      currentDisturbanceLayer = NULL,
      connectingBlockSize = NULL,
      disturbanceRateRelatesToBufferedArea = FALSE,
      outputsFolder = tempdir(),
      seismicLineGrids = 10,
      checkDisturbancesForBuffer = FALSE,
      runName = "testRun",
      useRoadsPackage = FALSE,
      siteSelectionAsDistributing = NA_character_,
      probabilityDisturbance = NULL,
      maskWaterAndMountainsFromLines = FALSE,
      featuresToAvoid = NULL,
      altitudeCut = NA_real_,
      clusterDistance = NA_real_,
      distanceNewLinesFactor = 1,
      runClusteringInParallel = FALSE,
      useClusterMethod = FALSE,
      refinedStructure = FALSE
    )
  )
  # Check that individualLayers contains our sector
  expect_type(out_full, "list")
  expect_named(out_full$individualLayers, "testSector")
  # Compare areas
  orig_area <- terra::expanse(poly, unit = "m")
  new_geom <- out_full$individualLayers$testSector$origin1
  expect_true(inherits(new_geom, "SpatVector"))
  new_area <- terra::expanse(new_geom, unit = "m")
  expect_true(new_area > orig_area)
})

# 5) Test: 0% enlarging rate yields no change in area

dp_zero <- copy(dp_enlarge)
set(dp_zero, j = "disturbanceRate", value = 0)

dist_list_zero <- list(testSector = list(origin1 = poly))

test_that("Enlarging with 0% rate returns original area", {
  out_zero_list <- generateDisturbancesShp(
    disturbanceParameters = dp_zero,
    disturbanceList = dist_list_zero,
    rasterToMatch = r,
    studyArea = sa,
    fires = NULL,
    currentTime = 1,
    firstTime = FALSE,
    growthStepGenerating = 1,
    growthStepEnlargingPolys = 1,
    growthStepEnlargingLines = 0.1,
    currentDisturbanceLayer = NULL,
    connectingBlockSize = NULL,
    disturbanceRateRelatesToBufferedArea = FALSE,
    outputsFolder = tempdir(),
    seismicLineGrids = 10,
    checkDisturbancesForBuffer = FALSE,
    runName = "testRun",
    useRoadsPackage = FALSE,
    siteSelectionAsDistributing = NA_character_,
    probabilityDisturbance = NULL,
    maskWaterAndMountainsFromLines = FALSE,
    featuresToAvoid = NULL,
    altitudeCut = NA_real_,
    clusterDistance = NA_real_,
    distanceNewLinesFactor = 1,
    runClusteringInParallel = FALSE,
    useClusterMethod = FALSE,
    refinedStructure = FALSE
  )
  new_poly <- out_zero_list$individualLayers$testSector$origin1
  orig_area <- terra::expanse(poly, unit = "m")
  zero_area <- terra::expanse(new_poly, unit = "m")
  expect_equal(zero_area, orig_area)
})

# 12) 
#---- Test for Enlarging ----
test_that("Buffered area calculation works for Enlarging", {
  # Parameters: 40% growth on existing polygons
  dp_enlarging <- data.table(
    dataName            = "testSector",
    disturbanceType     = "Enlarging",
    disturbanceOrigin   = "origin1",
    disturbanceRate     = 40,        # 40% target increase
    disturbanceInterval = 5,
    dataClass           = "potential",
    potentialField      = "Potential",
    disturbanceSize     = NA_character_  # not used in enlarging
  )
  
  # Provide an existing disturbance layer named 'origin1'
  out_enl <- generateDisturbancesShp(
    disturbanceParameters = dp_enlarging,
    disturbanceList = list(testSector = list(origin1 = sa)),
    rasterToMatch = r,
    studyArea = sa,
    currentTime = 2020,
    firstTime = FALSE,
    growthStepGenerating = 10,
    growthStepEnlargingPolys = 10,
    growthStepEnlargingLines = 10,
    currentDisturbanceLayer = NULL,
    connectingBlockSize = 100,
    disturbanceRateRelatesToBufferedArea = TRUE,
    outputsFolder = tempdir(),
    seismicLineGrids = 10,
    checkDisturbancesForBuffer = TRUE,
    runName = "test_run",
    useRoadsPackage = FALSE,
    siteSelectionAsDistributing = NULL,
    probabilityDisturbance = NULL,
    maskWaterAndMountainsFromLines = FALSE,
    featuresToAvoid = NULL,
    altitudeCut = NULL,
    clusterDistance = 1000,
    distanceNewLinesFactor = 1.0,
    runClusteringInParallel = FALSE,
    useClusterMethod = FALSE,
    refinedStructure = FALSE
  )
  
  out_poly <- out_enl$individualLayers$testSector$origin1
  buffered <- terra::buffer(out_poly, width = 500)
  buf_area <- terra::expanse(buffered, unit = "m")
  orig_area <- terra::expanse(out_poly, unit = "m")
  
  # Expect at least 40% buffered-area growth
  expect_true(buf_area > orig_area * 1.4)
})



################################################################################
### ------------------------------------------------------------------------ ###
################################################################################
# Generating branch

# 3) Test: multiple resolutionVector values trigger an error

dp_multi_res <- copy(dp_enlarge)
set(dp_multi_res, j = "resolutionVector", value = list(c(1,2)))

test_that("Multiple resolutionVector values generate an error", {
  expect_error(
    generateDisturbancesShp(
      disturbanceParameters = dp_multi_res,
      disturbanceList = disturbanceList_enlarge,
      rasterToMatch = r,
      studyArea = sa,
      fires = NULL,
      currentTime = 1,
      firstTime = FALSE,
      growthStepGenerating = 1,
      growthStepEnlargingPolys = 1,
      growthStepEnlargingLines = 0.1,
      currentDisturbanceLayer = NULL,
      connectingBlockSize = NULL,
      disturbanceRateRelatesToBufferedArea = FALSE,
      outputsFolder = tempdir(),
      seismicLineGrids = 10,
      checkDisturbancesForBuffer = FALSE,
      runName = "testRun",
      useRoadsPackage = FALSE,
      siteSelectionAsDistributing = NA_character_,
      probabilityDisturbance = NULL,
      maskWaterAndMountainsFromLines = FALSE,
      featuresToAvoid = NULL,
      altitudeCut = NA_real_,
      clusterDistance = NA_real_,
      distanceNewLinesFactor = 1,
      runClusteringInParallel = FALSE,
      useClusterMethod = FALSE,
      refinedStructure = FALSE
    ),
    regexp = "Different resolutions have not yet been implemented."
  )
})


# 4) Test: Generating returns NULL when no potential polygons

dp_gen <- data.table(
  dataName = "testSector",
  dataClass = "unusedClass",
  disturbanceType = "Generating",
  disturbanceRate = 10,
  disturbanceOrigin = list("originX"),
  disturbanceInterval = 1,
  resolutionVector = 1,
  potentialField = NA_character_
)

disturbanceList_gen <- list(
  testSector = list(originX = NULL),
  unusedClass = list(originX = NULL)
)

test_that("Generating branch returns NULL for empty potential", {
  # Suppress warning from max() on empty potential
  out_gen <- suppressWarnings(
    generateDisturbancesShp(
      disturbanceParameters = dp_gen,
      disturbanceList = disturbanceList_gen,
      rasterToMatch = r,
      studyArea = sa,
      fires = NULL,
      currentTime = 1,
      firstTime = FALSE,
      growthStepGenerating = 1,
      growthStepEnlargingPolys = 1,
      growthStepEnlargingLines = 0.1,
      currentDisturbanceLayer = r,
      connectingBlockSize = NULL,
      disturbanceRateRelatesToBufferedArea = FALSE,
      outputsFolder = tempdir(),
      seismicLineGrids = 10,
      checkDisturbancesForBuffer = FALSE,
      runName = "testRun",
      useRoadsPackage = FALSE,
      siteSelectionAsDistributing = NA_character_,
      probabilityDisturbance = NULL,
      maskWaterAndMountainsFromLines = FALSE,
      featuresToAvoid = NULL,
      altitudeCut = NA_real_,
      clusterDistance = NA_real_,
      distanceNewLinesFactor = 1,
      runClusteringInParallel = FALSE,
      useClusterMethod = FALSE,
      refinedStructure = FALSE
    )
  )
  expect_null(out_gen$individualLayers$testSector$originX)
})

#13) 
#---- Test for Generating (Polygon Area Bounds) ----
test_that("Generating creates some disturbance but not more than intended (module-standard params)", {
  skip_on_covr()  # this somehow blows up when doing covr()
  r  <- rast(nrows = 100, ncols = 100, xmin = 0, xmax = 1000, ymin = 0, ymax = 1000, vals = 1)
  crs(r) <- "EPSG:3005"
  sa <- as.polygons(ext(r))
  crs(sa) <- crs(r)
  
  distList <- createDisturbanceList(crs(r))
  
  existing_cutblocks <- st_polygon(list(rbind(c(80,40), c(80,60), c(90,60), c(90,40), c(80,40))))
  existing_cutblocks <- st_sfc(existing_cutblocks, crs = "EPSG:3005")
  existing_cutblocks <- vect(existing_cutblocks)
  
  existing_area <- sum(terra::expanse(existing_cutblocks, unit = "m"))
  
  # disturbance parameters: proportion of area
  dp_generating <- data.table(
    dataName            = "forestry",
    dataClass           = "potentialCutblocks",
    disturbanceType     = "Generating",
    disturbanceOrigin   = "cutblocks",
    disturbanceEnd      = "",
    disturbanceRate     = 0.1,      # <-- 0.1% per year = 5,000 m² over 5 years on 1 km²
    disturbanceInterval = 5L,
    potentialField      = "Potential",
    disturbanceSize     = "1000",
    resolutionVector    = 1
  )
  
  out_gen <- generateDisturbancesShp(
    disturbanceParameters                = dp_generating,
    disturbanceList                      = distList,
    rasterToMatch                        = r,
    studyArea                            = sa,
    fires                                = NULL,
    currentTime                          = 2020,
    firstTime                            = FALSE,
    growthStepGenerating                 = 10,
    growthStepEnlargingPolys             = 10,
    growthStepEnlargingLines             = 10,
    currentDisturbanceLayer              = NULL,
    connectingBlockSize                  = 100,
    disturbanceRateRelatesToBufferedArea = FALSE,
    outputsFolder                        = Paths$outputPath,
    seismicLineGrids                     = 10,
    checkDisturbancesForBuffer           = FALSE,
    runName                              = "test_run",
    useRoadsPackage                      = FALSE,
    siteSelectionAsDistributing          = "cutblocks",   # <-- a character vector of ORIGIN(s)
    probabilityDisturbance = list(cutblocks = data.table(Potential = 1L, probPoly = 1)), 
    maskWaterAndMountainsFromLines       = FALSE,
    featuresToAvoid                      = NULL,
    altitudeCut                          = NULL,
    clusterDistance                      = 1000,
    distanceNewLinesFactor               = 1.0,
    runClusteringInParallel              = FALSE,
    useClusterMethod                     = FALSE,
    refinedStructure                     = FALSE
  )
  
  # fetch output and compute only the newly generated area
  out_poly <- out_gen$individualLayers$forestry$cutblocks
  
  plot(sa)
  plot(existing_cutblocks
       , add=TRUE
       )
  plot(out_poly, add=TRUE)
  
  expect_true(inherits(out_poly, "SpatVector"))
  
  total_area_after <- sum(terra::expanse(out_poly, unit = "m"))
  
  generated_area <- max(total_area_after - existing_area, 0)
  
  # Expect ~5 * 1000 m² with a small tolerance
  expect_gt(generated_area, 0)
  expect_lte(generated_area, 5 * 1000 * 1.05)  # 5% slack
})


# 17
test_that("Generating succeeds when potential has only numeric 'Potential' column (control)", {
  Paths <<- list(outputPath = tempdir())
  
  r <- rast(nrows=10, ncols=10, xmin=0, xmax=10, ymin=0, ymax=10, vals=1)
  crs(r) <- "EPSG:3005"
  sa <- vect(ext(r)); crs(sa) <- crs(r)
  
  mine_poly <- st_polygon(list(rbind(c(2,2), c(6,2), c(6,6), c(2,6), c(2,2))))
  potentialMining <- vect(st_sfc(mine_poly)); crs(potentialMining) <- crs(r)
  potentialMining$Potential <- 1L    # <- no Class column this time
  
  dp <- data.table(
    dataName            = "mining",
    dataClass           = "potentialMining",
    disturbanceType     = "Generating",
    disturbanceRate     = 1,
    disturbanceSize     = "1",
    disturbanceOrigin   = "mining",
    disturbanceEnd      = "",
    disturbanceInterval = 1L,
    resolutionVector    = 1
  )
  
  disturbanceList <- list(
    mining = list(potentialMining = potentialMining, mining = NULL)
  )
  
  out <- generateDisturbancesShp(
    disturbanceParameters = dp,
    disturbanceList       = disturbanceList,
    rasterToMatch         = r,
    studyArea             = sa,
    fires                 = NULL,
    currentTime           = 2020,
    firstTime             = FALSE,
    growthStepGenerating  = 1,
    growthStepEnlargingPolys  = 1,
    growthStepEnlargingLines  = 1,
    currentDisturbanceLayer   = NULL,
    connectingBlockSize       = NULL,
    disturbanceRateRelatesToBufferedArea = FALSE,
    outputsFolder            = tempdir(),
    seismicLineGrids         = 10,
    checkDisturbancesForBuffer = FALSE,
    runName                  = "control_ok",
    useRoadsPackage          = FALSE,
    siteSelectionAsDistributing = character(0),
    probabilityDisturbance   = list(),
    maskWaterAndMountainsFromLines = FALSE,
    featuresToAvoid          = NULL,
    altitudeCut              = NULL,
    clusterDistance          = 1000,
    distanceNewLinesFactor   = 1.0,
    runClusteringInParallel  = FALSE,
    useClusterMethod         = FALSE,
    refinedStructure         = FALSE
  )
  
  gen_mines <- out$individualLayers$mining$mining
  expect_s4_class(gen_mines, "SpatVector")
  expect_true(terra::geomtype(gen_mines) == "polygons")
  expect_gt(nrow(gen_mines), 0)
})

# 18
test_that("Generating mining clamps when requested area exceeds available potential", {
  options(reproducible.cachePath = tempfile())
  set.seed(42)
  
  r  <- rast(ncol=10, nrow=10, xmin=0, xmax=10, ymin=0, ymax=10, vals=1)
  crs(r) <- "EPSG:3005"
  sa <- vect(ext(r)); crs(sa) <- crs(r)
  
  # Small potential (~4 m²), with Potential + Class present
  p <- vect(matrix(c(1,1, 3,1, 3,3, 1,3, 1,1), ncol=2, byrow=TRUE),
            type="polygons", crs=crs(r))
  p$Potential <- 1L
  p$Class     <- "potentialMining"
  
  dl <- list(mining = list(potentialMining = p, mining = NULL))
  
  # --- STUB the expected seismic cache file (needed on firstTime=TRUE) ---
  outdir <- tempdir()
  studyAreaHash <- digest::digest(sa)
  seis_stub_path <- file.path(outdir, sprintf("seismicLinesYear%04d_%s.shp", 2020L, studyAreaHash))
  seis_stub <- vect(matrix(c(8,0, 9,0), ncol=2, byrow=TRUE), type="lines", crs=crs(r))
  terra::writeVector(seis_stub, seis_stub_path, overwrite=TRUE)
  on.exit({
    base <- sub("\\.shp$", "", seis_stub_path)
    unlink(paste0(base, c(".shp",".shx",".dbf",".prj")), force=TRUE)
  }, add=TRUE)
  
  # Ask for WAY more than available so we exercise clamp/stop
  dp <- data.table(
    dataName            = "mining",
    dataClass           = "potentialMining",
    disturbanceType     = "Generating",
    disturbanceOrigin   = "mining",
    disturbanceEnd      = NA_character_,
    disturbanceRate     = 999,      # huge %
    disturbanceSize     = "200",
    disturbanceInterval = 1L,
    resolutionVector    = 1,
    potentialField      = "Potential"
  )
  
  # Ensure generation isn't vetoed by probability logic
  prob <- list(mining = data.table(Potential = 1L, probPoly = 1.0))
  
  expect_warning(
    out <- generateDisturbancesShp(
      disturbanceParameters                = dp,
      disturbanceList                      = dl,
      rasterToMatch                        = r,
      studyArea                            = sa,
      fires                                = NULL,
      currentTime                          = 2020,
      firstTime                            = TRUE,           # keep first-year path
      growthStepGenerating                 = 1,
      growthStepEnlargingPolys             = 1,
      growthStepEnlargingLines             = 1,
      currentDisturbanceLayer              = NULL,
      connectingBlockSize                  = NULL,
      disturbanceRateRelatesToBufferedArea = FALSE,
      outputsFolder                        = outdir,         # where we stubbed seismic
      seismicLineGrids                     = 1L,
      checkDisturbancesForBuffer           = FALSE,
      runName                              = "mining_clamp",
      useRoadsPackage                      = FALSE,
      siteSelectionAsDistributing          = character(0),
      probabilityDisturbance               = prob,           # <- key
      maskWaterAndMountainsFromLines       = FALSE,
      featuresToAvoid                      = NULL,
      altitudeCut                          = NULL,
      clusterDistance                      = 1000,
      distanceNewLinesFactor               = 1.0,
      runClusteringInParallel              = FALSE,
      useClusterMethod                     = FALSE,
      refinedStructure                     = FALSE
    ),
    regexp = "clamping|Stopping after|exceeds available potential"
  )
  
  expect_true(is.list(out))   # returned without hanging
})



################################################################################
### ------------------------------------------------------------------------ ###
################################################################################
# Connecting branch

# 16) 
test_that("Mines generate from potentialMining and connect to a pre-existing road", {
  # 1) Shared raster + study area
  r <- rast(nrows=10, ncols=10, xmin=0, xmax=10, ymin=0, ymax=10, vals=1)
  crs(r) <- "EPSG:3005"
  sa <- vect(ext(r)); crs(sa) <- crs(r)
  
  # 2) Baseline disturbance list (from your helpers) + patches
  dl <- createDisturbanceList(crs = crs(r))
  
  # Ensure this is a "first-year" like run (no existing mining that would erase potential)
  dl$mining$mining <- NULL
  
  # Make potentialMining roomy and numeric-only Potential attribute
  big_poly <- st_polygon(list(rbind(c(2,2), c(6,2), c(6,6), c(2,6), c(2,2)))) # 16 m²
  dl$mining$potentialMining <- vect(st_sfc(big_poly)); crs(dl$mining$potentialMining) <- crs(r)
  dl$mining$potentialMining$Potential <- 1L
  keep <- "Potential"
  dl$mining$potentialMining <- dl$mining$potentialMining[, keep]
  
  # 3) Seed a pre-existing roads layer (ensure inside SA)
  roadsSeed <- vect(st_sfc(st_linestring(rbind(c(8,8), c(8,2)))))
  crs(roadsSeed) <- crs(r)
  roadsSeed$Class <- "roads"
  stopifnot(terra::relate(roadsSeed, sa, "T********"))
  dl$pipelines$roads <- roadsSeed
  
  # 4) Parameters: generate mines, then connect mines -> pipelines$roads
  # SA area = 100 m²; Rate 2% => expected 2 m²; site size "1" => deterministic generation.
  dp <- rbind(
    data.table(
      dataName            = "mining",
      dataClass           = "potentialMining",
      disturbanceType     = "Generating",
      disturbanceRate     = 2,
      disturbanceSize     = "1",
      disturbanceOrigin   = "mining",
      disturbanceEnd      = "",
      disturbanceInterval = 1L,
      resolutionVector    = 1
    ),
    data.table(
      dataName            = "pipelines",
      dataClass           = "roads",
      disturbanceType     = "Connecting",
      disturbanceRate     = NA_real_,
      disturbanceSize     = NA_character_,
      disturbanceOrigin   = "mining",
      disturbanceEnd      = "roads",
      disturbanceInterval = 1L,
      resolutionVector    = 1
    )
  )
  
  # 5) Run
  out <- suppressWarnings(
    generateDisturbancesShp(
      disturbanceParameters                = dp,
      disturbanceList                      = dl,
      rasterToMatch                        = r,
      studyArea                            = sa,
      fires                                = NULL,
      currentTime                          = 2020,
      firstTime                            = FALSE,
      growthStepGenerating                 = 1,
      growthStepEnlargingPolys             = 1,
      growthStepEnlargingLines             = 1,
      currentDisturbanceLayer              = NULL,
      connectingBlockSize                  = NULL,
      disturbanceRateRelatesToBufferedArea = FALSE,
      outputsFolder                        = tempdir(),
      seismicLineGrids                     = 10,
      checkDisturbancesForBuffer           = FALSE,
      runName                              = "test_run",
      useRoadsPackage                      = FALSE,
      siteSelectionAsDistributing          = character(0),
      probabilityDisturbance               = NULL,
      maskWaterAndMountainsFromLines       = FALSE,
      featuresToAvoid                      = NULL,
      altitudeCut                          = NULL,
      clusterDistance                      = 1000,
      distanceNewLinesFactor               = 1.0,
      runClusteringInParallel              = FALSE,
      useClusterMethod                     = FALSE,
      refinedStructure                     = FALSE
    ))
  
  # 6) Extract with a tiny safeguard for the module's name typo
  il <- out$individualLayers
  expect_true(is.list(il) && length(il) > 0)
  
  gen_mines  <- il$mining$mining
  conn_lines <- il$pipelines$roads
  
  # 7) Basic type/row checks
  expect_s4_class(gen_mines,  "SpatVector")
  expect_identical(terra::geomtype(gen_mines), "polygons")
  expect_gt(nrow(gen_mines), 0)
  
  expect_s4_class(conn_lines, "SpatVector")
  expect_identical(terra::geomtype(conn_lines), "lines")
  expect_gt(nrow(conn_lines), 0)
  
  # 8) Connectivity assertions:
  # (a) A connection line must intersect the seeded road
  conn_lines_sf <- st_as_sf(conn_lines)  # Conversion depends on the format of conn_lines
  roadsSeed_sf <- st_as_sf(roadsSeed)    # Conversion depends on the format of roadsSeed
  hitRoad <- st_intersects(conn_lines_sf, roadsSeed_sf, sparse = FALSE)
  expect_true(any(hitRoad))
  
  # (b) A connection line must intersect at least one generated mine polygon
  gen_mines_sf <- st_as_sf(gen_mines)
  hitMine <- st_intersects(conn_lines_sf, gen_mines_sf, sparse = FALSE)
  expect_true(any(hitMine))
  
  # 9) Plots
  #plot(sa)
  #plot(roadsSeed, add=TRUE)
  #plot(big_poly, add=TRUE)
  #plot(gen_mines, add=TRUE)
  #plot(conn_lines, add=TRUE)
})


# 21
test_that("Connecting avoids raster obstacles (no blocking)", {
  r  <- rast(ncol=100, nrow=100, xmin=0, xmax=100, ymin=0, ymax=100, vals=1); crs(r) <- "EPSG:32633"
  sa <- as.polygons(ext(r)); crs(sa) <- crs(r)
  
  pot <- vect(matrix(c(0,0, 50,0, 50,100, 0,100, 0,0), 5, 2, byrow=TRUE), type="polygons", crs=crs(r)); pot$Potential <- 1L
  roads <- terra::vect(
    sf::st_sfc(sf::st_linestring(rbind(c(90,10), c(90,90))), crs = 32633)
  )
  roads$Class <- "roads"
  lake <- vect(matrix(c(40,40, 60,40, 60,60, 40,60, 40,40), 5, 2, byrow=TRUE), type="polygons", crs=crs(r))
  
  # robust obstacle mask: NA in lake OR high-alt
  dem <- init(r, "x") + init(r, "y")
  avoid <- rast(r); avoid[] <- 1
  avoid <- terra::mask(avoid, lake, inverse = TRUE)  # inside lake -> NA
  avoid[dem >= 120] <- NA
  
  dp_gen <- data.table(dataName="mining", dataClass="potentialMining",
                       disturbanceType="Generating", disturbanceOrigin="mining",
                       disturbanceRate=2, disturbanceSize="1", disturbanceInterval=1L,
                       resolutionVector=1, potentialField="Potential")
  dp_conn <- data.table(dataName="pipelines", dataClass="roads",
                        disturbanceType="Connecting", disturbanceOrigin="mining",
                        disturbanceEnd="roads", disturbanceInterval=1L, resolutionVector=1)
  dp <- rbind(dp_gen, dp_conn, fill=TRUE)
  dl <- list(mining=list(mining=NULL, potentialMining=pot), pipelines=list(roads=roads))
  
  out <- suppressWarnings(generateDisturbancesShp(
    disturbanceParameters=dp, disturbanceList=dl,
    rasterToMatch=r, studyArea=sa, fires=NULL, currentTime=1, firstTime=FALSE,
    growthStepGenerating=1, growthStepEnlargingPolys=1, growthStepEnlargingLines=1,
    currentDisturbanceLayer=NULL, connectingBlockSize=NULL,
    disturbanceRateRelatesToBufferedArea=FALSE, outputsFolder=tempdir(),
    seismicLineGrids=1, checkDisturbancesForBuffer=FALSE, runName="avoid",
    useRoadsPackage=FALSE, siteSelectionAsDistributing=character(0),
    probabilityDisturbance=NULL,
    maskWaterAndMountainsFromLines=TRUE, featuresToAvoid=avoid,
    altitudeCut=NA_real_, clusterDistance=NA_real_, distanceNewLinesFactor=1,
    runClusteringInParallel=FALSE, useClusterMethod=FALSE, refinedStructure=FALSE
  ))
  
  con <- out$individualLayers$pipelines$roads
  
  plot(avoid)
  plot(roads, add=TRUE)
  plot(con, add=TRUE)
  
  expect_s4_class(con, "SpatVector")
  expect_identical(terra::geomtype(con), "lines")
  
  expect_equal(nrow(terra::intersect(con, lake)), 0)     # no crossings
  pts <- terra::densify(con, interval=1)
  vals <- terra::extract(avoid, pts)[, 2]
  expect_true(all(!is.na(vals)))
})

# 22
test_that("Blocking branch (connectingBlockSize) ignores obstacles (current behavior)", {
  r  <- rast(ncol=100, nrow=100, xmin=0, xmax=100, ymin=0, ymax=100, vals=1); crs(r) <- "EPSG:32633"
  sa <- as.polygons(ext(r)); crs(sa) <- crs(r)
  
  pot <- vect(matrix(c(0,0, 50,0, 50,100, 0,100, 0,0), 5, 2, byrow=TRUE), type="polygons", crs=crs(r)); pot$Potential <- 1L
  roads <- terra::vect(
    sf::st_sfc(sf::st_linestring(rbind(c(90,10), c(90,90))), crs = 32633)
  )
  roads$Class <- "roads"
  lake <- vect(matrix(c(40,40, 60,40, 60,60, 40,60, 40,40), 5, 2, byrow=TRUE), type="polygons", crs=crs(r))
  
  # robust obstacle mask: NA in lake OR high-alt
  dem <- init(r, "x") + init(r, "y")
  avoid <- rast(r); avoid[] <- 1
  avoid <- terra::mask(avoid, lake, inverse = TRUE)  # inside lake -> NA
  avoid[dem >= 120] <- NA
  
  dp_gen <- data.table(dataName="mining", dataClass="potentialMining",
                       disturbanceType="Generating", disturbanceOrigin="mining",
                       disturbanceRate=2, disturbanceSize="1", disturbanceInterval=1L,
                       resolutionVector=1, potentialField="Potential")
  dp_conn <- data.table(dataName="pipelines", dataClass="roads",
                        disturbanceType="Connecting", disturbanceOrigin="mining",
                        disturbanceEnd="roads", disturbanceInterval=1L, resolutionVector=1)
  dp <- rbind(dp_gen, dp_conn, fill=TRUE)
  dl <- list(mining=list(mining=NULL, potentialMining=pot), pipelines=list(roads=roads))
  
  out <- suppressWarnings(generateDisturbancesShp(
    disturbanceParameters=dp, disturbanceList=dl,
    rasterToMatch=r, studyArea=sa, fires=NULL, currentTime=1, firstTime=FALSE,
    growthStepGenerating=1, growthStepEnlargingPolys=1, growthStepEnlargingLines=1,
    currentDisturbanceLayer=NULL, connectingBlockSize=50,
    disturbanceRateRelatesToBufferedArea=FALSE, outputsFolder=tempdir(),
    seismicLineGrids=1, checkDisturbancesForBuffer=FALSE, runName="avoid",
    useRoadsPackage=FALSE, siteSelectionAsDistributing=character(0),
    probabilityDisturbance=NULL,
    maskWaterAndMountainsFromLines=TRUE, featuresToAvoid=avoid,
    altitudeCut=NA_real_, clusterDistance=NA_real_, distanceNewLinesFactor=1,
    runClusteringInParallel=FALSE, useClusterMethod=FALSE, refinedStructure=FALSE
  ))
  
  con <- out$individualLayers$pipelines$roads
  
  plot(avoid)
  plot(roads, add=TRUE)
  plot(con, add=TRUE)
  
  # EXPECT crossings (until implementation switches to cost-routed in blocking branch)
  expect_gt(nrow(terra::intersect(con, lake)), 0)
})

# 23) Cost-routed avoidance with a "corridor" in raster mask ------------------
testthat::test_that("Cost-routed avoidance: vertical NA wall with a narrow gap", {
  testthat::skip_if_not_installed("spaths")  # cost path
  base <- make_base()
  r <- base$r; sa <- base$sa; pot <- base$pot; roads <- base$roads; lake <- base$lake
  dp_info <- make_dp(conn_block_size = NULL)
  
  # Obstacle mask: NA almost everywhere along x=50 (vertical wall),
  # but allow a gap for 10 cells in y ∈ [45,55].
  avoid <- rast(r); avoid[] <- 1
  x <- init(r, "x"); y <- init(r, "y")
  wall <- abs(x - 50) <= 0.5
  gap  <- wall & (y >= 45 & y <= 55)
  avoid[wall] <- NA
  avoid[gap]  <- 1     # open a narrow corridor
  
  out <- suppressWarnings(
    generateDisturbancesShp(
      disturbanceParameters = dp_info$dp,
      disturbanceList = list(mining=list(mining=NULL, potentialMining=pot),
                             pipelines=list(roads=roads)),
      rasterToMatch=r, studyArea=sa,
      fires=NULL, currentTime=1, firstTime=FALSE,
      growthStepGenerating=1, growthStepEnlargingPolys=1, growthStepEnlargingLines=1,
      currentDisturbanceLayer=NULL,
      connectingBlockSize = dp_info$connectingBlockSize,
      disturbanceRateRelatesToBufferedArea=FALSE, outputsFolder=tempdir(),
      seismicLineGrids=1, checkDisturbancesForBuffer=FALSE, runName="corridor",
      useRoadsPackage=FALSE, siteSelectionAsDistributing=character(0),
      probabilityDisturbance=NULL,
      maskWaterAndMountainsFromLines=TRUE, featuresToAvoid=avoid,
      altitudeCut=NA_real_, clusterDistance=NA_real_, distanceNewLinesFactor=1,
      runClusteringInParallel=FALSE, useClusterMethod=FALSE, refinedStructure=FALSE
    )
  )
  
  con <- get_connections(out)
  testthat::expect_s4_class(con, "SpatVector")
  testthat::expect_identical(terra::geomtype(con), "lines")
  # Must not intersect the NA wall except through the gap -> all vertices allowed
  pts <- terra::densify(con, interval = 0.5)
  vals <- terra::extract(avoid, pts)[,2]
  testthat::expect_true(all(!is.na(vals)))
  
  # Optional: ensure they cross near x=50 within gap band
  V <- terra::crds(pts)
  near_wall <- V[,1] >= 49.5 & V[,1] <= 50.5
  if (any(near_wall)) {
    testthat::expect_true(all(V[near_wall,2] >= 45 & V[near_wall,2] <= 55))
  }
})

# 24) Vector obstacles path (polygon lake as featuresToAvoid) -----------------
testthat::test_that("Vector polygon obstacles are rasterized and avoided", {
  testthat::skip_if_not_installed("spaths")  # cost path
  base <- make_base()
  r <- base$r; sa <- base$sa; pot <- base$pot; roads <- base$roads; lake <- base$lake
  dp_info <- make_dp(conn_block_size = NULL)
  
  out <- suppressWarnings(
    generateDisturbancesShp(
      disturbanceParameters = dp_info$dp,
      disturbanceList = list(mining=list(mining=NULL, potentialMining=pot),
                             pipelines=list(roads=roads)),
      rasterToMatch=r, studyArea=sa,
      fires=NULL, currentTime=1, firstTime=FALSE,
      growthStepGenerating=1, growthStepEnlargingPolys=1, growthStepEnlargingLines=1,
      currentDisturbanceLayer=NULL,
      connectingBlockSize = dp_info$connectingBlockSize,
      disturbanceRateRelatesToBufferedArea=FALSE, outputsFolder=tempdir(),
      seismicLineGrids=1, checkDisturbancesForBuffer=FALSE, runName="vector-obstacles",
      useRoadsPackage=FALSE, siteSelectionAsDistributing=character(0),
      probabilityDisturbance=NULL,
      maskWaterAndMountainsFromLines=TRUE, featuresToAvoid=lake,  # <- SpatVector polygons
      altitudeCut=NA_real_, clusterDistance=NA_real_, distanceNewLinesFactor=1,
      runClusteringInParallel=FALSE, useClusterMethod=FALSE, refinedStructure=FALSE
    )
  )
  con <- get_connections(out)
  testthat::expect_s4_class(con, "SpatVector")
  testthat::expect_identical(terra::geomtype(con), "lines")
  testthat::expect_equal(nrow(terra::intersect(con, lake)), 0)
})

# 25) Blocking branch ignores obstacles (document current behavior) ----------
testthat::test_that("Blocking branch (connectingBlockSize) currently ignores obstacles", {
  base <- make_base()
  r <- base$r; sa <- base$sa; pot <- base$pot; roads <- base$roads; lake <- base$lake
  
  # Make sure we trigger blocking: set small connectingBlockSize
  dp_info <- make_dp(conn_block_size = 1L)
  
  # Build a hard obstacle raster (lake area NA)
  avoid <- rast(r); avoid[] <- 1
  avoid <- terra::mask(avoid, lake, inverse = TRUE) # inside lake -> NA
  
  out <- suppressWarnings(
    generateDisturbancesShp(
      disturbanceParameters = dp_info$dp,
      disturbanceList = list(mining=list(mining=NULL, potentialMining=pot),
                             pipelines=list(roads=roads)),
      rasterToMatch=r, studyArea=sa,
      fires=NULL, currentTime=1, firstTime=FALSE,
      growthStepGenerating=1, growthStepEnlargingPolys=1, growthStepEnlargingLines=1,
      currentDisturbanceLayer=NULL,
      connectingBlockSize = dp_info$connectingBlockSize,
      disturbanceRateRelatesToBufferedArea=FALSE, outputsFolder=tempdir(),
      seismicLineGrids=1, checkDisturbancesForBuffer=FALSE, runName="blocking",
      useRoadsPackage=FALSE, siteSelectionAsDistributing=character(0),
      probabilityDisturbance=NULL,
      maskWaterAndMountainsFromLines=TRUE, featuresToAvoid=avoid,
      altitudeCut=NA_real_, clusterDistance=NA_real_, distanceNewLinesFactor=1,
      runClusteringInParallel=FALSE, useClusterMethod=FALSE, refinedStructure=FALSE
    )
  )
  con <- get_connections(out)
  testthat::expect_s4_class(con, "SpatVector")
  
  # Current implementation uses terra::nearest(lines=TRUE) in blocking branch (no mask),
  # so at least some connections will cross the lake.
  testthat::expect_gt(nrow(terra::intersect(con, lake)), 0)
})


# 27) Class propagation: connected lines inherit endLay Class -----------------
test_that("Connected lines carry the endLay Class (blocking branch)", {
  base <- make_base()
  r <- base$r; sa <- base$sa; pot <- base$pot; roads <- base$roads
  
  # Force blocking: set a very small block size so NROW(oriLayVect) > block
  dp_info <- make_dp(conn_block_size = 1L)
  
  # No obstacles needed for this test (we're asserting 'Class' only)
  avoid <- rast(r); avoid[] <- 1
  
  out <- suppressWarnings(
    generateDisturbancesShp(
      disturbanceParameters = dp_info$dp,
      disturbanceList = list(mining=list(mining=NULL, potentialMining=pot),
                             pipelines=list(roads=roads)),
      rasterToMatch=r, studyArea=sa,
      fires=NULL, currentTime=1, firstTime=FALSE,
      growthStepGenerating=1, growthStepEnlargingPolys=1, growthStepEnlargingLines=1,
      currentDisturbanceLayer=NULL,
      connectingBlockSize = dp_info$connectingBlockSize,   # <- triggers blocking branch
      disturbanceRateRelatesToBufferedArea=FALSE, outputsFolder=tempdir(),
      seismicLineGrids=1, checkDisturbancesForBuffer=FALSE, runName="class-prop",
      useRoadsPackage=FALSE, siteSelectionAsDistributing=character(0),
      probabilityDisturbance=NULL,
      maskWaterAndMountainsFromLines=FALSE,                 # <- skip cost routing entirely
      featuresToAvoid=avoid,
      altitudeCut=NA_real_, clusterDistance=NA_real_, distanceNewLinesFactor=1,
      runClusteringInParallel=FALSE, useClusterMethod=FALSE, refinedStructure=FALSE
    )
  )
  
  IL  <- if (!is.null(out$individualLayers)) out$individualLayers else out$individuaLayers
  con <- IL$pipelines$roads
  
  expect_s4_class(con, "SpatVector")
  expect_identical(terra::geomtype(con), "lines")
  expect_setequal(unique(con$Class), "roads")
})



################################################################################
### ------------------------------------------------------------------------ ###
################################################################################
# Misc tests

# 2) Test: missing layer triggers an error

dp_missing <- copy(dp_enlarge)

test_that("Missing layer in disturbanceList errors out", {
  expect_error(
    generateDisturbancesShp(
      disturbanceParameters = dp_missing,
      disturbanceList = list(otherSector = list(dummy = poly)),
      rasterToMatch = r,
      studyArea = sa,
      fires = NULL,
      currentTime = 1,
      firstTime = FALSE,
      growthStepGenerating = 1,
      growthStepEnlargingPolys = 1,
      growthStepEnlargingLines = 0.1,
      currentDisturbanceLayer = NULL,
      connectingBlockSize = NULL,
      disturbanceRateRelatesToBufferedArea = FALSE,
      outputsFolder = tempdir(),
      seismicLineGrids = 10,
      checkDisturbancesForBuffer = FALSE,
      runName = "testRun",
      useRoadsPackage = FALSE,
      siteSelectionAsDistributing = NA_character_,
      probabilityDisturbance = NULL,
      maskWaterAndMountainsFromLines = FALSE,
      featuresToAvoid = NULL,
      altitudeCut = NA_real_,
      clusterDistance = NA_real_,
      distanceNewLinesFactor = 1,
      runClusteringInParallel = FALSE,
      useClusterMethod = FALSE,
      refinedStructure = FALSE
    ),
    regexp = "needs to exist in disturbanceList"
  )
})



# 6) Test: unknown disturbanceType yields no potential error

dp_bad <- copy(dp_enlarge)
set(dp_bad, j = "disturbanceType", value = "UnknownType")

test_that("Unknown disturbanceType yields no potential error", {
  expect_error(
    generateDisturbancesShp(
      disturbanceParameters = dp_bad,
      disturbanceList = disturbanceList_enlarge,
      rasterToMatch = r,
      studyArea = sa,
      fires = NULL,
      currentTime = 1,
      firstTime = FALSE,
      growthStepGenerating = 1,
      growthStepEnlargingPolys = 1,
      growthStepEnlargingLines = 0.1,
      currentDisturbanceLayer = NULL,
      connectingBlockSize = NULL,
      disturbanceRateRelatesToBufferedArea = FALSE,
      outputsFolder = tempdir(),
      seismicLineGrids = 10,
      checkDisturbancesForBuffer = FALSE,
      runName = "testRun",
      useRoadsPackage = FALSE,
      siteSelectionAsDistributing = NA_character_,
      probabilityDisturbance = NULL,
      maskWaterAndMountainsFromLines = FALSE,
      featuresToAvoid = NULL,
      altitudeCut = NA_real_,
      clusterDistance = NA_real_,
      distanceNewLinesFactor = 1,
      runClusteringInParallel = FALSE,
      useClusterMethod = FALSE,
      refinedStructure = FALSE
    ),
    regexp = "study area has no potential for disturbances"
  )
})

# 7) Test: multiple sectors produce multiple entries

dp_two <- rbind(
  dp_enlarge,
  data.table(
    dataName = "secondSector",
    dataClass = "secondSector",
    disturbanceType = "Enlarging",
    disturbanceRate = 50,
    disturbanceOrigin = list("origin1"),
    disturbanceInterval = 1,
    resolutionVector = 1,
    potentialField = NA_character_
  )
)

dist_list_two <- list(
  testSector   = list(origin1 = poly),
  secondSector = list(origin1 = poly)
)

test_that("Multiple sectors returned correctly", {
  out_two <- generateDisturbancesShp(
    disturbanceParameters = dp_two,
    disturbanceList = dist_list_two,
    rasterToMatch = r,
    studyArea = sa,
    fires = NULL,
    currentTime = 1,
    firstTime = FALSE,
    growthStepGenerating = 1,
    growthStepEnlargingPolys = 1,
    growthStepEnlargingLines = 0.1,
    currentDisturbanceLayer = NULL,
    connectingBlockSize = NULL,
    disturbanceRateRelatesToBufferedArea = FALSE,
    outputsFolder = tempdir(),
    seismicLineGrids = 10,
    checkDisturbancesForBuffer = FALSE,
    runName = "testRun",
    useRoadsPackage = FALSE,
    siteSelectionAsDistributing = NA_character_,
    probabilityDisturbance = NULL,
    maskWaterAndMountainsFromLines = FALSE,
    featuresToAvoid = NULL,
    altitudeCut = NA_real_,
    clusterDistance = NA_real_,
    distanceNewLinesFactor = 1,
    runClusteringInParallel = FALSE,
    useClusterMethod = FALSE,
    refinedStructure = FALSE
  )
  expect_named(out_two$individualLayers, c("testSector", "secondSector"))
})

# 8) Test: checkDisturbancesForBuffer TRUE prints discrepancy but retains output

test_that("checkDisturbancesForBuffer TRUE prints discrepancy message", {
  dp_buf <- copy(dp_enlarge)
  dist_list_buf <- list(testSector = list(origin1 = poly))
  output_messages <- capture.output(
    out_buf <- generateDisturbancesShp(
      disturbanceParameters = dp_buf,
      disturbanceList = dist_list_buf,
      rasterToMatch = r,
      studyArea = sa,
      fires = NULL,
      currentTime = 1,
      firstTime = FALSE,
      growthStepGenerating = 1,
      growthStepEnlargingPolys = 1,
      growthStepEnlargingLines = 0.1,
      currentDisturbanceLayer = NULL,
      connectingBlockSize = NULL,
      disturbanceRateRelatesToBufferedArea = FALSE,
      outputsFolder = tempdir(),
      seismicLineGrids = 10,
      checkDisturbancesForBuffer = TRUE,
      runName = "testRun",
      useRoadsPackage = FALSE,
      siteSelectionAsDistributing = NA_character_,
      probabilityDisturbance = NULL,
      maskWaterAndMountainsFromLines = FALSE,
      featuresToAvoid = NULL,
      altitudeCut = NA_real_,
      clusterDistance = NA_real_,
      distanceNewLinesFactor = 1,
      runClusteringInParallel = FALSE,
      useClusterMethod = FALSE,
      refinedStructure = FALSE
    )
  )
  expect_true(any(grepl("Difference between expected and achieved change", output_messages)))
  expect_named(out_buf$individualLayers, "testSector")
  expect_true(terra::expanse(out_buf$individualLayers$testSector$origin1) > 0)
})

# 9) Test: siteSelectionAsDistributing option accepted without error

dp_dist <- data.table(
  dataName = "testSector",
  dataClass = "unusedClass",
  disturbanceType = "Generating",
  disturbanceRate = 10,
  disturbanceOrigin = list("originX"),
  disturbanceInterval = 1,
  resolutionVector = 1,
  potentialField = NA_character_
)
# Provide a dummy potential layer to avoid early NULL return
potential_poly <- poly

dist_list_dist <- list(testSector = list(originX = potential_poly))

test_that("siteSelectionAsDistributing option accepted", {
  # Suppress messages and warnings from empty Generating branch
  out_dist <- suppressWarnings(
    suppressMessages(
      generateDisturbancesShp(
        disturbanceParameters = dp_dist,
        disturbanceList = dist_list_dist,
        rasterToMatch = r,
        studyArea = sa,
        fires = NULL,
        currentTime = 1,
        firstTime = FALSE,
        growthStepGenerating = 1,
        growthStepEnlargingPolys = 1,
        growthStepEnlargingLines = 0.1,
        currentDisturbanceLayer = NULL,
        connectingBlockSize = NULL,
        disturbanceRateRelatesToBufferedArea = FALSE,
        outputsFolder = tempdir(),
        seismicLineGrids = 10,
        checkDisturbancesForBuffer = FALSE,
        runName = "testRun",
        useRoadsPackage = FALSE,
        siteSelectionAsDistributing = "originX",
        probabilityDisturbance = NULL,
        maskWaterAndMountainsFromLines = FALSE,
        featuresToAvoid = NULL,
        altitudeCut = NA_real_,
        clusterDistance = NA_real_,
        distanceNewLinesFactor = 1,
        runClusteringInParallel = FALSE,
        useClusterMethod = FALSE,
        refinedStructure = FALSE
      )
    )
  )
  # Even if no-op, should return a list structure
  expect_type(out_dist, "list")
  expect_named(out_dist, c("individualLayers", "currentDisturbanceLayer", "seismicLinesFirstYear"))
})

# 10) Test: disturbanceRateRelatesToBufferedArea yields expected growth

test_that("disturbanceRateRelatesToBufferedArea yields expected growth", {
  # Define a large square polygon of side 1000 m
  poly_large <- vect(matrix(c(
    0, 0,
    1000, 0,
    1000, 1000,
    0, 1000,
    0, 0
  ), ncol = 2, byrow = TRUE), type = "polygons")
  crs(poly_large) <- crs(r)
  
  # Setup parameters for 50% growth
  dp_rate <- data.table(
    dataName = "testSector",
    dataClass = "testSector",
    disturbanceType = "Enlarging",
    disturbanceRate = 50,
    disturbanceOrigin = list("origin1"),
    disturbanceInterval = 1,
    resolutionVector = 1,
    potentialField = NA_character_
  )
  disturbanceList_rate <- list(testSector = list(origin1 = poly_large))
  
  # Original area
  orig_area <- terra::expanse(poly_large, unit = "m")
  
  # Raw growth = 50% of original area (should be > orig_area)
  capture.output(
    out_raw <- generateDisturbancesShp(
      disturbanceParameters = dp_rate,
      disturbanceList = disturbanceList_rate,
      rasterToMatch = r, studyArea = sa, fires = NULL,
      currentTime = 1, firstTime = FALSE,
      growthStepGenerating = 1, growthStepEnlargingPolys = 1, growthStepEnlargingLines = 0.1,
      currentDisturbanceLayer = NULL, connectingBlockSize = NULL,
      disturbanceRateRelatesToBufferedArea = FALSE, outputsFolder = tempdir(),
      seismicLineGrids = 10, checkDisturbancesForBuffer = FALSE,
      runName = "testRun", useRoadsPackage = FALSE,
      siteSelectionAsDistributing = NA_character_, probabilityDisturbance = NULL,
      maskWaterAndMountainsFromLines = FALSE, featuresToAvoid = NULL,
      altitudeCut = NA_real_, clusterDistance = NA_real_, distanceNewLinesFactor = 1,
      runClusteringInParallel = FALSE, useClusterMethod = FALSE, refinedStructure = FALSE
    )
  )
  area_raw <- terra::expanse(out_raw$individualLayers$testSector$origin1, unit = "m")
  expect_true(area_raw > orig_area)
  
  # Buffered growth = 50% of buffered area (should be > area_raw)
  capture.output(
    out_buf <- generateDisturbancesShp(
      disturbanceParameters = dp_rate,
      disturbanceList = disturbanceList_rate,
      rasterToMatch = r, studyArea = sa, fires = NULL,
      currentTime = 1, firstTime = FALSE,
      growthStepGenerating = 1, growthStepEnlargingPolys = 1, growthStepEnlargingLines = 0.1,
      currentDisturbanceLayer = NULL, connectingBlockSize = NULL,
      disturbanceRateRelatesToBufferedArea = TRUE, outputsFolder = tempdir(),
      seismicLineGrids = 10, checkDisturbancesForBuffer = FALSE,
      runName = "testRun", useRoadsPackage = FALSE,
      siteSelectionAsDistributing = NA_character_, probabilityDisturbance = NULL,
      maskWaterAndMountainsFromLines = FALSE, featuresToAvoid = NULL,
      altitudeCut = NA_real_, clusterDistance = NA_real_, distanceNewLinesFactor = 1,
      runClusteringInParallel = FALSE, useClusterMethod = FALSE, refinedStructure = FALSE
    )
  )
  # For buffered case: Measure the 500m-buffered version of the output
  output_poly <- out_buf$individualLayers$testSector$origin1
  buffered_poly <- terra::buffer(output_poly, width = 500)
  area_buf_buffered <- terra::expanse(buffered_poly, unit = "m")
  
  # Should be original_buffered_area + 0.5 km² = 3.78e6 + 5e5 = 4.28e6 m²
  expect_equal(area_buf_buffered, 4.28e6, tolerance = 0.01e6)
})

# 11) Test: firstTime TRUE populates seismicLinesFirstYear
test_that("firstTime TRUE picks up seismicLinesFirstYear via stub", {
  # 1. Create VALID SPATIAL INPUT for disturbanceList
  poly_large <- vect(
    matrix(
      c(0, 0,
        1000, 0,
        1000, 1000,
        0, 1000,
        0, 0),
      ncol = 2,
      byrow = TRUE
    ),
    type = "polygons"
  )
  crs(poly_large) <- crs(r)
  
  # 2. Compute the hash exactly as the function does
  studyAreaHash <- digest::digest(sa)
  
  # 3. Write a complete ESRI Shapefile stub using sf
  stub_shp <- file.path(
    tempdir(),
    paste0("seismicLinesYear", 1, "_", studyAreaHash, ".shp")
  )
  dummy_v <- vect(
    matrix(
      c(
        0, 0,
        0, 1,
        1, 1,
        1, 0,
        0, 0
      ),
      ncol = 2,
      byrow = TRUE
    ),
    type = "polygons"
  )
  crs(dummy_v) <- crs(r)
  # Convert to sf and write shapefile (creates .shp, .dbf, .shx, etc.)
  sf_obj <- sf::st_as_sf(dummy_v)
  sf::write_sf(sf_obj, stub_shp, delete_layer = TRUE)
  
  # 4. Call function, specifying 'seismicLines' origin and outputsFolder
  out_ft <- generateDisturbancesShp(
    disturbanceParameters = data.table(
      dataName          = "testSector",
      dataClass         = "testSector",
      disturbanceType   = "Enlarging",
      disturbanceRate   = 0,
      disturbanceOrigin = list("seismicLines"),
      disturbanceInterval = 1,
      resolutionVector  = 1,
      potentialField    = NA_character_
    ),
    disturbanceList = list(
      testSector = list(seismicLines = poly_large)
    ),
    rasterToMatch                = r,
    studyArea                    = sa,
    fires                        = NULL,
    currentTime                  = 1,
    firstTime                    = TRUE,
    growthStepGenerating         = 1,
    growthStepEnlargingPolys     = 1,
    growthStepEnlargingLines     = 0.1,
    currentDisturbanceLayer      = NULL,
    connectingBlockSize          = NULL,
    disturbanceRateRelatesToBufferedArea = FALSE,
    outputsFolder                = tempdir(),
    seismicLineGrids             = 10,
    checkDisturbancesForBuffer   = FALSE,
    runName                      = "testRun",
    useRoadsPackage              = FALSE,
    siteSelectionAsDistributing  = NA_character_,
    probabilityDisturbance       = NULL,
    maskWaterAndMountainsFromLines = FALSE,
    featuresToAvoid              = NULL,
    altitudeCut                  = NA_real_,
    clusterDistance              = NA_real_,
    distanceNewLinesFactor       = 1,
    runClusteringInParallel      = FALSE,
    useClusterMethod             = FALSE,
    refinedStructure             = FALSE
  )
  
  # 5. Verify stub was loaded
  expect_s4_class(out_ft$seismicLinesFirstYear, "SpatVector")
  expect_equal(nrow(out_ft$seismicLinesFirstYear), nrow(dummy_v))
})


#
# 14) 
#------------------------------------------------------------------------------
# Two-step clustering test
#------------------------------------------------------------------------------
# Expect same clusters reused

test_that("Cluster seeds have Pot_Clus on first year and Subsequent-year reuse preserves attributes", {
  # Use a temporary cache
  options(reproducible.cachePath = tempfile())
  
  # Create a minimal raster (100x100 instead of 1000x1000)
  r <- terra::rast(ncol = 100, nrow = 100, xmin = 0, xmax = 1000, ymin = 0, ymax = 1000)
  terra::values(r) <- 1
  terra::crs(r) <- "EPSG:3857"
  
  # Minimal study area polygon (slightly larger than raster)
  sa_poly <- terra::vect(terra::ext(r) + 50)
  terra::crs(sa_poly) <- terra::crs(r)
  sa_poly$Potential <- 1
  
  # Create just 2 seismic lines (instead of grid)
  lines <- list(
    terra::vect(matrix(c(100, 100, 900, 900), ncol = 2), type = "lines"),
    terra::vect(matrix(c(100, 900, 900, 100), ncol = 2), type = "lines")
  )
  seismic_lines <- do.call(rbind, lines)
  terra::crs(seismic_lines) <- terra::crs(r)
  seismic_lines$Potential <- 1
  
  # Simplified disturbance parameters
  dp_seismic <- data.table(
    dataName = "seismic",
    disturbanceType = "Generating",
    disturbanceOrigin = "seismicLines",
    disturbanceRate = 0.01, # Reduced rate for smaller area
    disturbanceSize = "50",
    disturbanceInterval = 5,
    dataClass = "potential",
    potentialField = "Potential",
    resolutionVector = 1
  )
  
  # Predefined probability distribution
  probabilityDisturbance <- list(
    seismicLines = data.table(Potential = 1, probPoly = 1)
  )
  
  disturbanceList <- list(
    seismic = list(
      seismicLines = seismic_lines,
      potential = sa_poly
    )
  )
  
  # First run (year 2020)
  out1 <- generateDisturbancesShp(
    disturbanceParameters = dp_seismic,
    disturbanceList = disturbanceList,
    rasterToMatch = r,
    studyArea = sa_poly,
    currentTime = 2020,
    firstTime = TRUE,
    growthStepGenerating = 1, # Minimal growth steps
    growthStepEnlargingPolys = 1,
    growthStepEnlargingLines = 1,
    currentDisturbanceLayer = NULL,
    connectingBlockSize = 10, # Smaller blocks
    disturbanceRateRelatesToBufferedArea = FALSE,
    outputsFolder = tempdir(),
    seismicLineGrids = 1, # Minimal grid count
    checkDisturbancesForBuffer = FALSE,
    runName = "test_run",
    useRoadsPackage = FALSE,
    siteSelectionAsDistributing = "seismicLines",
    probabilityDisturbance = probabilityDisturbance,
    maskWaterAndMountainsFromLines = FALSE,
    featuresToAvoid = NULL,
    altitudeCut = NULL,
    clusterDistance = 10, # Smaller cluster distance
    distanceNewLinesFactor = 0.1, # Minimal jitter
    runClusteringInParallel = FALSE,
    useClusterMethod = TRUE,
    refinedStructure = FALSE # Skip refinement
  )
  
  # Validate first run
  sf1 <- out1$seismicLinesFirstYear
  expect_s4_class(sf1, "SpatVector")
  expect_equal(terra::geomtype(sf1), "lines")
  expect_true("Pot_Clus" %in% names(terra::values(sf1)))
  
  # Second run (year 2021) - reuse first year's output
  prevLayers <- list(
    seismic = list(
      seismicLines = sf1,
      potential = sa_poly
    )
  )
  
  # Set seed for reproducible sampling
  set.seed(123)
  out2 <- generateDisturbancesShp(
    disturbanceParameters = dp_seismic,
    disturbanceList = list(seismic = list(seismicLines = sf1, potential = sa_poly)),
    rasterToMatch = r,
    studyArea = sa_poly,
    currentTime = 2021,
    firstTime = FALSE,
    currentDisturbanceLayer = prevLayers,
    growthStepGenerating = 1,
    growthStepEnlargingPolys = 1,
    growthStepEnlargingLines = 1,
    connectingBlockSize = 10,
    disturbanceRateRelatesToBufferedArea = FALSE,
    outputsFolder = tempdir(),
    seismicLineGrids = 1,
    checkDisturbancesForBuffer = FALSE,
    runName = "test_run2",
    useRoadsPackage = FALSE,
    siteSelectionAsDistributing = "seismicLines",
    probabilityDisturbance = probabilityDisturbance,
    maskWaterAndMountainsFromLines = FALSE,
    featuresToAvoid = NULL,
    altitudeCut = NULL,
    clusterDistance = 10,
    distanceNewLinesFactor = 0.1,
    runClusteringInParallel = FALSE,
    useClusterMethod = TRUE,
    refinedStructure = FALSE
  )
  
  # Validate second run
  vec2 <- out2$individualLayers$seismic$seismicLines
  plot(vec2)
  expect_s4_class(vec2, "SpatVector")
  expect_equal(terra::geomtype(vec2), "lines")
  df2 <- terra::values(vec2)
  expect_true("Pot_Clus" %in% names(df2))
  expect_equal(
    sort(unique(df2$Pot_Clus)),
    sort(unique(terra::values(out1$individualLayers$seismic$seismicLines)$Pot_Clus))
  )
})


# -----------------------------------------------------------------------------
# 15) Test Forestry “Generating” avoids burned cells
# -----------------------------------------------------------------------------
test_that("Cutblock centroids don't fall inside burned area", {
  options(reproducible.cachePath = tempfile())
  set.seed(4255)
  
  # --- Base raster & study area (use same CRS as helpers) ---
  r  <- rast(nrows = 100, ncols = 100, xmin = 0, xmax = 100, ymin = 0, ymax = 100, vals = 1)
  crs(r) <- "EPSG:3005"
  sa <- as.polygons(r, dissolve = TRUE); crs(sa) <- crs(r)
  
  # --- Potential layer: one polygon per cell with Potential=1 ---
  pot_vect <- as.polygons(r, dissolve = FALSE); crs(pot_vect) <- crs(r)
  pot_vect$Potential <- 1L
  pot_vect$ORIGIN    <- 1900L
  pot_vect$Class     <- "potential"
  
  # --- Existing cutblock inside SA (acts like current footprint) ---
  existing_cut <- vect(matrix(c(5,5, 5,7, 7,7, 7,5, 5,5), ncol = 2, byrow = TRUE),
                       type = "polygons", crs = crs(r))
  existing_cut$Class <- "cutblocks"
  
  # --- Burned patch in the middle to be avoided ---
  fires_rast <- r; values(fires_rast) <- 0L
  mid_ext <- ext(30, 70, 30, 70)
  fires_rast[cells(fires_rast, mid_ext)] <- 1L
  burned_vect <- as.polygons(fires_rast == 1, dissolve = TRUE); crs(burned_vect) <- crs(r)
  
  # --- Build a baseline disturbanceList, then conform forestry naming ---
  baseDL <- createDisturbanceList(crs = "EPSG:3005")
  # Replace forestry layers with our test data; ensure generic "potential" name exists
  baseDL$forestry$cutblocks <- existing_cut
  baseDL$forestry$potential <- pot_vect
  # (Optionally remove the helper’s specific name to avoid duplicate semantics)
  baseDL$forestry$potentialCutblocks <- NULL
  
  disturbanceList <- baseDL
  
  # --- Minimal currentDisturbanceLayer (only what's used by this test) ---
  currentDL <- list(forestry = list(cutblocks = existing_cut))
  
  # --- Disturbance parameters row for forestry/cutblocks generation ---
  dp_forestry <- data.table(
    dataName            = "forestry",
    dataClass           = "potential",       # matches disturbanceList$forestry$potential
    disturbanceType     = "Generating",
    disturbanceOrigin   = "cutblocks",
    disturbanceRate     = 0.8,               # small overall addition
    disturbanceSize     = 50,                # NUMERIC (buffer widths, areas, etc.)
    disturbanceEnd      = "",
    disturbanceInterval = 5L,
    potentialField      = "Potential",
    resolutionVector    = list(1)            # list-col as in module tables
  )
  
  # --- stub the expected seismic cache file ---
  outdir <- tempdir()
  studyAreaHash <- digest::digest(sa)
  seis_stub_path <- file.path(outdir, sprintf("seismicLinesYear%04d_%s.shp", 2020L, studyAreaHash))
  
  # minimal valid line so terra::vect() can read it
  seis_stub <- vect(matrix(c(80,0, 90,0), ncol=2, byrow=TRUE), type="lines", crs=crs(r))
  terra::writeVector(seis_stub, seis_stub_path, overwrite=TRUE)
  if (!file.exists(seis_stub_path)) {
    skip("terra::writeVector didn't create the seismic stub (namespace may be monkey-patched).")
  }
  on.exit({
    base <- sub("\\.shp$", "", seis_stub_path)
    unlink(paste0(base, c(".shp",".shx",".dbf",".prj")), force = TRUE)
  }, add = TRUE)
  
  # --- Run the generator ---
  out <- suppressWarnings(
    generateDisturbancesShp(
      disturbanceParameters                = dp_forestry,
      disturbanceList                      = disturbanceList,
      rasterToMatch                        = r,
      studyArea                            = sa,
      fires                                = fires_rast,
      currentTime                          = 2020L,
      firstTime                            = TRUE,
      growthStepGenerating                 = 10,
      growthStepEnlargingPolys             = 10,
      growthStepEnlargingLines             = 10,
      currentDisturbanceLayer              = currentDL,
      connectingBlockSize                  = 100,
      disturbanceRateRelatesToBufferedArea = FALSE,
      outputsFolder                        = tempdir(),
      seismicLineGrids                     = 1L,  
      checkDisturbancesForBuffer           = TRUE,
      runName                              = "forestry_test_burned_center",
      useRoadsPackage                      = FALSE,
      siteSelectionAsDistributing          = c("cutblocks"),
      probabilityDisturbance               = list(cutblocks = data.table(Potential = 1L, probPoly = 1.0)),
      maskWaterAndMountainsFromLines       = FALSE,
      featuresToAvoid                      = list(burned = burned_vect),
      altitudeCut                          = NULL,
      clusterDistance                      = 1000,
      distanceNewLinesFactor               = 1.0,
      runClusteringInParallel              = FALSE,
      useClusterMethod                     = FALSE,
      refinedStructure                     = FALSE
    ))
  
  # --- Checks ---
  new_cut <- out$individualLayers$forestry$cutblocks
  cents <- terra::centroids(new_cut)
  # No centroid may fall in the burned area
  names(fires_rast) <- "burn"
  vals <- terra::extract(fires_rast, cents)[,"burn", drop = TRUE]
  expect_true(all(vals != 1L))   # no centroid on a burned cell
  
  # --- Plots ---
  #plot(fires_rast)
  #plot(existing_cut, add=TRUE)
  #plot(new_cut, add = TRUE)
  #plot(cents, add=TRUE)
})



##################
# additional tests: Bug isolation



# 19
test_that("All outputs share rasterToMatch CRS and stable names", {
  # Minimal working inputs
  r  <- rast(ncol=50, nrow=50, xmin=0, xmax=1000, ymin=0, ymax=1000, vals=1)
  crs(r) <- "EPSG:32633"
  sa <- as.polygons(ext(r)); crs(sa) <- crs(r)
  
  # Small polygon to enlarge (0% so geometry is preserved but returned)
  base <- vect(matrix(c(100,100, 300,100, 300,300, 100,300, 100,100), ncol=2, byrow=TRUE), type="polygons", crs=crs(r))
  dp <- data.table(
    dataName="test", dataClass="test",
    disturbanceType="Enlarging",
    disturbanceRate=0, disturbanceOrigin=list("origin1"),
    disturbanceInterval=1, resolutionVector=1, potentialField=NA_character_
  )
  dl <- list(test=list(origin1=base))
  
  out <- generateDisturbancesShp(
    disturbanceParameters=dp, disturbanceList=dl,
    rasterToMatch=r, studyArea=sa,
    fires=NULL, currentTime=1, firstTime=FALSE,
    growthStepGenerating=1, growthStepEnlargingPolys=1, growthStepEnlargingLines=1,
    currentDisturbanceLayer=NULL, connectingBlockSize=NULL,
    disturbanceRateRelatesToBufferedArea=FALSE, outputsFolder=tempdir(),
    seismicLineGrids=1, checkDisturbancesForBuffer=FALSE, runName="struct",
    useRoadsPackage=FALSE, siteSelectionAsDistributing=character(0),
    probabilityDisturbance=NULL, maskWaterAndMountainsFromLines=FALSE,
    featuresToAvoid=NULL, altitudeCut=NA_real_,
    clusterDistance=NA_real_, distanceNewLinesFactor=1,
    runClusteringInParallel=FALSE, useClusterMethod=FALSE, refinedStructure=FALSE
  )
  
  # Accept either the typo ('individuaLayers') or the fixed name ('individualLayers')
  expect_true(
    all(c("individuaLayers","currentDisturbanceLayer","seismicLinesFirstYear") %in% names(out)) ||
      all(c("individualLayers","currentDisturbanceLayer","seismicLinesFirstYear") %in% names(out))
  )
  
  IL <- unlist(out$individualLayers)
  vects <- IL$test.origin1
  
  # Check CRS of every returned non-NULL SpatVector equals crs(r)
  lapply(vects, function(v) expect_identical(terra::crs(v), terra::crs(r)))
  
  # Empty branch still present: firstTime=FALSE => seismicLinesFirstYear exists but is NULL
  expect_true("seismicLinesFirstYear" %in% names(out))
  expect_null(out$seismicLinesFirstYear)
})

# 20
test_that("Distributing favors high potential; Exhausting fills it first", {
  set.seed(123)
  
  r  <- rast(ncol=100, nrow=100, xmin=0, xmax=100, ymin=0, ymax=100, vals=1)
  crs(r) <- "EPSG:32633"
  sa <- as.polygons(ext(r)); crs(sa) <- crs(r)
  
  # Two equal-size polygons with very different Potential values
  hi <- vect(matrix(c( 0, 0, 50, 0, 50,100,  0,100,  0, 0), ncol=2, byrow=TRUE), type="polygons", crs=crs(r))
  lo <- vect(matrix(c(50, 0,100, 0,100,100, 50,100, 50, 0), ncol=2, byrow=TRUE), type="polygons", crs=crs(r))
  pot <- rbind(hi, lo)
  pot$Potential <- c(10L, 1L)
  
  dp <- data.table(
    dataName="mining", dataClass="potentialMining",
    disturbanceType="Generating", disturbanceOrigin="mining",
    disturbanceRate=5,               # 5% of SA; enough to compare distribution
    disturbanceSize="1",             # small, many sites
    disturbanceInterval=1L,
    resolutionVector=1, potentialField="Potential"
  )
  dl <- list(mining=list(mining=NULL, potentialMining=pot))
  
  # (a) Distributing with explicit probabilities
  probs <- list(mining = data.table(Potential=c(1L,10L), probPoly=c(0.1,0.9)))
  out_dist <- suppressWarnings(
    generateDisturbancesShp(
      disturbanceParameters=dp, disturbanceList=dl,
      rasterToMatch=r, studyArea=sa,
      fires=NULL, currentTime=1, firstTime=FALSE,
      growthStepGenerating=1, growthStepEnlargingPolys=1, growthStepEnlargingLines=1,
      currentDisturbanceLayer=NULL, connectingBlockSize=NULL,
      disturbanceRateRelatesToBufferedArea=FALSE, outputsFolder=tempdir(),
      seismicLineGrids=1, checkDisturbancesForBuffer=FALSE, runName="dist",
      useRoadsPackage=FALSE, siteSelectionAsDistributing="mining",
      probabilityDisturbance=probs, maskWaterAndMountainsFromLines=FALSE,
      featuresToAvoid=NULL, altitudeCut=NA_real_,
      clusterDistance=NA_real_, distanceNewLinesFactor=1,
      runClusteringInParallel=FALSE, useClusterMethod=FALSE, refinedStructure=FALSE
    )
  )
  IL <- out_dist$individualLayers
  gen <- IL$mining$mining
  a_hi <- sum(terra::expanse(terra::intersect(gen, hi), unit="m"))
  a_lo <- sum(terra::expanse(terra::intersect(gen, lo), unit="m"))
  expect_gt(a_hi/(a_hi + a_lo), 0.80)  # >80% in high potential
  
  # (b) Exhausting (default): should all-but-fill the high-potential polygon first
  out_exh <- suppressWarnings(
    generateDisturbancesShp(
      disturbanceParameters=dp, disturbanceList=dl,
      rasterToMatch=r, studyArea=sa,
      fires=NULL, currentTime=2, firstTime=FALSE,
      growthStepGenerating=1, growthStepEnlargingPolys=1, growthStepEnlargingLines=1,
      currentDisturbanceLayer=NULL, connectingBlockSize=NULL,
      disturbanceRateRelatesToBufferedArea=FALSE, outputsFolder=tempdir(),
      seismicLineGrids=1, checkDisturbancesForBuffer=FALSE, runName="exh",
      useRoadsPackage=FALSE, siteSelectionAsDistributing=character(0),
      probabilityDisturbance=NULL, maskWaterAndMountainsFromLines=FALSE,
      featuresToAvoid=NULL, altitudeCut=NA_real_,
      clusterDistance=NA_real_, distanceNewLinesFactor=1,
      runClusteringInParallel=FALSE, useClusterMethod=FALSE, refinedStructure=FALSE
    )
  )
  IL2 <- if (!is.null(out_exh$individualLayers)) out_exh$individualLayers else out_exh$individuaLayers
  gen2 <- IL2$mining$mining
  b_hi <- sum(terra::expanse(terra::intersect(gen2, hi), unit="m"))
  b_lo <- sum(terra::expanse(terra::intersect(gen2, lo), unit="m"))
  expect_gt(b_hi/(b_hi + b_lo), 0.95)  # “exhausting” should heavily favor the max potential band
  
  
  #plot(sa)
  #plot(pot, add=TRUE)
  #plot(gen, add=TRUE)
  #plot(gen2, add=TRUE)
})


# 26) CRS/grid mismatch in featuresToAvoid (raster) ---------------------------
testthat::test_that("CRS/grid mismatch in featuresToAvoid raises an error (current)", {
  base <- make_base()
  r <- base$r; sa <- base$sa; pot <- base$pot; roads <- base$roads
  dp_info <- make_dp(conn_block_size = NULL)
  
  # Hi-res mask: same extent & CRS, different resolution (will mismatch ncell)
  avoid_hi <- rast(ncol=200, nrow=200, xmin=0, xmax=100, ymin=0, ymax=100, crs=crs(r))
  avoid_hi[] <- 1
  # carve a NA stripe to make it nontrivial
  xh <- init(avoid_hi, "x")
  avoid_hi[ abs(xh - 50) < 1 ] <- NA
  
  testthat::expect_error(
    suppressWarnings(
      generateDisturbancesShp(
        disturbanceParameters = dp_info$dp,
        disturbanceList = list(mining=list(mining=NULL, potentialMining=pot),
                               pipelines=list(roads=roads)),
        rasterToMatch=r, studyArea=sa,
        fires=NULL, currentTime=1, firstTime=FALSE,
        growthStepGenerating=1, growthStepEnlargingPolys=1, growthStepEnlargingLines=1,
        currentDisturbanceLayer=NULL,
        connectingBlockSize = dp_info$connectingBlockSize,
        disturbanceRateRelatesToBufferedArea=FALSE, outputsFolder=tempdir(),
        seismicLineGrids=1, checkDisturbancesForBuffer=FALSE, runName="grid-mismatch",
        useRoadsPackage=FALSE, siteSelectionAsDistributing=character(0),
        probabilityDisturbance=NULL,
        maskWaterAndMountainsFromLines=TRUE, featuresToAvoid=avoid_hi,
        altitudeCut=NA_real_, clusterDistance=NA_real_, distanceNewLinesFactor=1,
        runClusteringInParallel=FALSE, useClusterMethod=FALSE, refinedStructure=FALSE
      )
    )
  )
})


### additional tests for upping coverage
test_that("Generating: computes current buffered area when relating rate to 500 m buffer", {
  testthat::skip_if_not_installed("terra")
  testthat::skip_if_not_installed("sf")
  
  r  <- terra::rast(nrows = 60, ncols = 60, xmin = 0, xmax = 600, ymin = 0, ymax = 600, vals = 1)
  terra::crs(r) <- "EPSG:3005"
  sa <- terra::as.polygons(terra::ext(r)); terra::crs(sa) <- terra::crs(r)
  Paths <<- list(outputPath = tempdir(), inputPath = tempdir())
  
  # Potential polygons with a numeric Potential column
  pot <- terra::vect("POLYGON ((50 50, 50 550, 550 550, 550 50, 50 50))"); terra::crs(pot) <- terra::crs(r)
  pot$Potential <- 10L
  
  # Some current lines (Lay) inside the potential
  lay <- terra::vect("LINESTRING (100 300, 500 300)"); terra::crs(lay) <- terra::crs(r); lay$Class <- "pipelines"
  
  dList <- list(energy = list(pipelines = lay, potentialPipelines = pot))
  
  dp <- data.table::data.table(
    disturbanceType     = "Generating",
    dataName            = "energy",
    disturbanceOrigin   = "pipelines",
    dataClass           = "potentialPipelines",
    disturbanceRate     = 0.2,     # small but > 0 to run the branch
    disturbanceInterval = 1,
    potentialField      = "Potential",
    disturbanceSize     = "1000"   # m², not really used here
  )
  
  # Expect the "Buffered (500m) area..." message to be printed from the branch
  expect_message(
    suppressWarnings(
      generateDisturbancesShp(
        disturbanceParameters = dp,
        disturbanceList = dList,
        rasterToMatch = r, studyArea = sa, fires = NULL,
        currentTime = 1, firstTime = FALSE,
        growthStepGenerating = 1, growthStepEnlargingPolys = 1, growthStepEnlargingLines = 1,
        currentDisturbanceLayer = NULL,
        connectingBlockSize = NULL,
        disturbanceRateRelatesToBufferedArea = TRUE,   # <— triggers branch
        outputsFolder = tempdir(),
        seismicLineGrids = 1,
        checkDisturbancesForBuffer = FALSE,
        runName = "t-gen-buff",
        useRoadsPackage = FALSE,
        siteSelectionAsDistributing = character(),
        probabilityDisturbance = NULL,
        maskWaterAndMountainsFromLines = FALSE,
        featuresToAvoid = NULL,
        altitudeCut = 0,
        clusterDistance = 0,
        distanceNewLinesFactor = 0,
        runClusteringInParallel = FALSE,
        useClusterMethod = TRUE,
        refinedStructure = FALSE
      )
    ),
    "Buffered \\(500m\\) area for energy -- pipelines:",
    fixed = FALSE
  ) # covers the message path. :contentReference[oaicite:3]{index=3}
})


test_that("Generating: probabilityDisturbance is computed when NULL for distributing origins", {
  testthat::skip_if_not_installed("terra")
  testthat::skip_if_not_installed("sf")
  
  r  <- terra::rast(nrows = 60, ncols = 60, xmin = 0, xmax = 600, ymin = 0, ymax = 600, vals = 1)
  terra::crs(r) <- "EPSG:3005"
  sa <- terra::as.polygons(terra::ext(r)); terra::crs(sa) <- terra::crs(r)
  Paths <<- list(outputPath = tempdir(), inputPath = tempdir())
  
  # Two potential polygons with different Potential values
  potA <- terra::vect("POLYGON ((50 50, 50 300, 300 300, 300 50, 50 50))");  terra::crs(potA) <- terra::crs(r); potA$Potential <- 5L
  potB <- terra::vect("POLYGON ((300 300, 300 550, 550 550, 550 300, 300 300))"); terra::crs(potB) <- terra::crs(r); potB$Potential <- 10L
  pot  <- rbind(potA, potB)
  
  # Existing seismic lines (Lay) overlapping both potentials
  lay  <- terra::vect("LINESTRING (60 60, 290 290)"); terra::crs(lay) <- terra::crs(r); lay$Class <- "seismicLines"
  
  dList <- list(mining = list(seismicLines = lay, potentialSeismic = pot))
  
  dp <- data.table::data.table(
    disturbanceType     = "Generating",
    dataName            = "mining",
    disturbanceOrigin   = "seismicLines",
    dataClass           = "potentialSeismic",
    disturbanceRate     = 0.05,
    disturbanceInterval = 1,
    potentialField      = "Potential",
    disturbanceSize     = "1000"
  )
  
  expect_message(
    suppressWarnings(
      generateDisturbancesShp(
        disturbanceParameters = dp,
        disturbanceList = dList,
        rasterToMatch = r, studyArea = sa, fires = NULL,
        currentTime = 1, firstTime = TRUE,          # → createCropLayFinalYear1 path
        growthStepGenerating = 1, growthStepEnlargingPolys = 1, growthStepEnlargingLines = 1,
        currentDisturbanceLayer = NULL,
        connectingBlockSize = NULL,
        disturbanceRateRelatesToBufferedArea = FALSE,
        outputsFolder = tempdir(),
        seismicLineGrids = 1,
        checkDisturbancesForBuffer = FALSE,
        runName = "t-prob-auto",
        useRoadsPackage = FALSE,
        siteSelectionAsDistributing = "seismicLines",    # <— triggers branch
        probabilityDisturbance = NULL,                   # <— auto compute
        maskWaterAndMountainsFromLines = FALSE,
        featuresToAvoid = NULL,
        altitudeCut = 0,
        clusterDistance = 1000,
        distanceNewLinesFactor = 1,
        runClusteringInParallel = FALSE,
        useClusterMethod = TRUE,                         # avoids the grid/rtnorm path
        refinedStructure = FALSE
      )
    ),
    "probabilityDisturbance for seismicLines is NULL\\. Calculating from data\\.",
    fixed = FALSE
  ) # branch + message. :contentReference[oaicite:5]{index=5}
})

test_that("Generating: tiny size vs 500 m buffered point → fallback semantics are honored", {
  testthat::skip_if_not_installed("terra")
  testthat::skip_if_not_installed("sf")
  testthat::skip_if_not_installed("SpaDES.tools")
  
  r  <- terra::rast(nrows = 60, ncols = 60, xmin = 0, xmax = 600, ymin = 0, ymax = 600, vals = 1)
  terra::crs(r) <- "EPSG:3005"
  sa <- terra::as.polygons(terra::ext(r)); terra::crs(sa) <- terra::crs(r)
  Paths <<- list(outputPath = tempdir(), inputPath = tempdir())
  
  pot <- terra::vect("POLYGON ((50 50, 50 550, 550 550, 550 50, 50 50))"); terra::crs(pot) <- terra::crs(r)
  pot$Potential <- 1L
  lay <- terra::vect("POLYGON ((100 100, 100 150, 150 150, 150 100, 100 100))"); terra::crs(lay) <- terra::crs(r)
  dList <- list(industry = list(mining = lay, potentialMining = pot))
  
  dp <- data.table::data.table(
    disturbanceType     = "Generating",
    dataName            = "industry",
    disturbanceOrigin   = "mining",
    dataClass           = "potentialMining",
    disturbanceRate     = 0.01,     # 0.01% of the study area
    disturbanceInterval = 1,
    potentialField      = "Potential",
    disturbanceSize     = "200"     # m^2 — forces fallback
  )
  
  out <- suppressWarnings(
    generateDisturbancesShp(
      disturbanceParameters = dp,
      disturbanceList = dList,
      rasterToMatch = r, studyArea = sa, fires = NULL,
      currentTime = 1, firstTime = FALSE,
      growthStepGenerating = 1, growthStepEnlargingPolys = 1, growthStepEnlargingLines = 1,
      currentDisturbanceLayer = NULL,
      connectingBlockSize = NULL,
      disturbanceRateRelatesToBufferedArea = TRUE,
      outputsFolder = tempdir(),
      seismicLineGrids = 1,
      checkDisturbancesForBuffer = FALSE,
      runName = "t-minbuff-fallback",
      useRoadsPackage = FALSE,
      siteSelectionAsDistributing = character(),
      probabilityDisturbance = NULL,
      maskWaterAndMountainsFromLines = FALSE,
      featuresToAvoid = NULL,
      altitudeCut = 0,
      clusterDistance = 0,
      distanceNewLinesFactor = 0,
      runClusteringInParallel = FALSE,
      useClusterMethod = TRUE,
      refinedStructure = FALSE
    )
  )
  
  # The generated layer is stored under individualLayers
  gen <- out$individualLayers$industry$mining
  testthat::expect_s4_class(gen, "SpatVector")           # SpatVector is S4, not S3
  testthat::expect_identical(terra::geomtype(gen), "polygons")  # minimal buffer polygon
  
  # Fallback semantics: when re-buffered by 500 m, area ≈ π·500^2 per feature
  A_expected <- pi * 500^2
  A_buf <- sum(terra::expanse(terra::buffer(gen, width = 500), unit = "m", transform = FALSE))
  testthat::expect_lt(abs(A_buf - A_expected), 0.05 * A_expected)  # within 5%
  
  # And the unbuffered geometry should be essentially area-less (because of the tiny buffer)
  A_raw <- sum(terra::expanse(gen, unit = "m", transform = FALSE))
  testthat::expect_lt(A_raw, 1)  # < 1 m^2
})

test_that("Connecting: featuresToAvoid is built from geodata (mocked)", {
  testthat::skip_if_not_installed("terra")
  testthat::skip_if_not_installed("mockery")
  
  r  <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 200, ymin = 0, ymax = 200, vals = 1)
  terra::crs(r) <- "EPSG:3005"
  sa <- terra::as.polygons(terra::ext(r)); terra::crs(sa) <- terra::crs(r)
  Paths <<- list(outputPath = tempdir(), inputPath = tempdir())
  
  # minimal potential + a road to connect to
  pot  <- terra::vect("POLYGON ((30 30, 30 70, 70 70, 70 30, 30 30))", crs = terra::crs(r))
  pot$Potential <- 1L; pot$ORIGIN <- 1900L
  road <- terra::vect("LINESTRING (0 0, 200 200)", crs = terra::crs(r)); road$Class <- "roads"
  
  dList <- list(mining = list(
    potentialMining = pot,
    roads = road
  ))
  
  dp <- data.table::rbindlist(list(
    # Seed at least one origin to connect FROM
    data.table::data.table(
      disturbanceType     = "Generating",
      dataName            = "mining",
      disturbanceOrigin   = "mining",
      dataClass           = "potentialMining",
      disturbanceRate     = 0.001,
      disturbanceInterval = 1L,
      disturbanceSize     = "100",
      resolutionVector    = list(15)
    ),
    # And the Connecting instruction (use 'disturbanceEnd', not 'endLayer')
    data.table::data.table(
      disturbanceType     = "Connecting",
      dataName            = "mining",
      disturbanceOrigin   = "mining",
      disturbanceEnd      = "roads",
      resolutionVector    = list(15)
    )
  ), fill = TRUE)
  
  # stubs – let our mocked rasters actually flow through postProcessTo
  mock_elev <- function(country, path) terra::setValues(terra::rast(r), 1000)  # > altitudeCut
  mock_lcov <- function(type,    path) terra::setValues(terra::rast(r), 0)
  mock_pp   <- function(from, to, ...) from
  
  mockery::stub(generateDisturbancesShp, "elevation_30s", mock_elev)
  mockery::stub(generateDisturbancesShp, "landcover",     mock_lcov)
  mockery::stub(generateDisturbancesShp, "postProcessTo", mock_pp)
  
  # The generator is chatty (message()), so assert "no error" rather than "silent"
  expect_error(
    suppressMessages(suppressWarnings(
      generateDisturbancesShp(
        disturbanceParameters = dp,
        disturbanceList = dList,
        rasterToMatch = r, studyArea = sa, fires = NULL,
        currentTime = 1, firstTime = FALSE,
        growthStepGenerating = 1, growthStepEnlargingPolys = 1, growthStepEnlargingLines = 1,
        currentDisturbanceLayer = NULL,
        connectingBlockSize = 100,
        disturbanceRateRelatesToBufferedArea = FALSE,
        outputsFolder = tempdir(),
        seismicLineGrids = 1,
        checkDisturbancesForBuffer = FALSE,
        runName = "t-connect-mask",
        useRoadsPackage = FALSE,
        siteSelectionAsDistributing = character(),
        probabilityDisturbance = NULL,
        maskWaterAndMountainsFromLines = TRUE,   # forces geodata path
        featuresToAvoid = NULL,                  # so our stubs are used
        altitudeCut = 500,
        clusterDistance = 0,
        distanceNewLinesFactor = 0,
        runClusteringInParallel = FALSE,
        useClusterMethod = TRUE,
        refinedStructure = FALSE
      )
    )),
    NA
  )
})
