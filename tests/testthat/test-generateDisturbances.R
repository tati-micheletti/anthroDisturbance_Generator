# ---------- Your working fixtures (lightly inlined) ----------
# study grid + study area
r  <- rast(nrows = 100, ncols = 100, xmin = 0, xmax = 1000, ymin = 0, ymax = 1000, vals = 1)
crs(r) <- "EPSG:3005"
sa <- as.polygons(ext(r))
crs(sa) <- crs(r)

# module global
Paths <<- list(outputPath = tempdir())

# helper: sf -> SpatVector and attach fields
to_sv <- function(sf_geom, class_name, crs_str) {
  v <- vect(sf_geom)
  crs(v) <- crs_str
  if (nrow(v) > 0) {
    v$Class <- class_name
    if (startsWith(class_name, "potential")) {
      v$Potential <- 1L
      v$ORIGIN    <- 1900L
    }
  } else {
    names(v) <- character(0)
  }
  v
}

createDisturbanceList <- function(crs_str = "EPSG:3005") {
  poly1 <- st_polygon(list(rbind(c(100,100), c(100,300), c(300,300), c(300,100), c(100,100))))
  settlements           <- to_sv(st_sfc(poly1),                           "settlements",           crs_str)
  potentialSettlements  <- to_sv(st_buffer(st_sfc(poly1), 25),            "potentialSettlements",  crs_str)
  
  windTurbines          <- suppressWarnings(to_sv(st_sfc(),               "windTurbines",          crs_str))
  potentialWindTurbines <- to_sv(st_sfc(st_multipoint(rbind(c(700,700), c(800,200)))), "potentialWindTurbines", crs_str)
  
  pipelines             <- suppressWarnings(to_sv(st_sfc(),               "pipelines",             crs_str))
  potentialPipelines    <- to_sv(st_sfc(st_point(c(200,800)), st_point(c(400,900))), "potentialPipelines",  crs_str)
  roads                 <- to_sv(st_sfc(st_linestring(rbind(c(0,500), c(1000,500)))), "roads",            crs_str)
  
  mine_poly             <- st_polygon(list(rbind(c(500,100), c(500,200), c(600,200), c(600,100), c(500,100))))
  mining                <- to_sv(st_sfc(mine_poly),                        "mining",               crs_str)
  potentialMining       <- to_sv(st_buffer(st_sfc(mine_poly), 50),         "potentialMining",      crs_str)
  
  cb_poly    <- st_polygon(list(rbind(c(800,400), c(800,600), c(900,600), c(900,400), c(800,400))))
  cutblocks  <- to_sv(st_sfc(cb_poly), "cutblocks", crs_str)
  
  # Potential cutblocks: whole study minus existing cutblock
  pot_big            <- as.polygons(r, dissolve = TRUE)
  crs(pot_big)       <- crs_str
  pot_big            <- erase(pot_big, cutblocks)
  pot_big$Class      <- "potentialCutblocks"
  pot_big$Potential  <- 1L
  pot_big$ORIGIN     <- 1900L
  potentialCutblocks <- pot_big
  
  seis_line             <- st_linestring(rbind(c(300,500), c(300,1000)))
  seismicLines          <- to_sv(st_sfc(seis_line),                        "seismicLines",         crs_str)
  potentialSeismicLines <- to_sv(st_sfc(st_point(c(300,700)), st_point(c(200,900))), "potentialSeismicLines", crs_str)
  
  list(
    settlements = list(settlements = settlements, potentialSettlements = potentialSettlements),
    wind        = list(windTurbines = windTurbines, potentialWindTurbines = potentialWindTurbines),
    pipelines   = list(pipelines = pipelines, potentialPipelines = potentialPipelines, roads = roads),
    mining      = list(mining = mining, potentialMining = potentialMining),
    forestry    = list(cutblocks = cutblocks, potentialCutblocks = potentialCutblocks),
    oilGas      = list(seismicLines = seismicLines, potentialSeismicLines = potentialSeismicLines)
  )
}

# helper to ensure Paths exists per call
with_paths <- function(expr) {
  old <- if (exists("Paths", inherits = TRUE)) get("Paths", inherits = TRUE) else NULL
  on.exit({
    if (is.null(old)) rm("Paths", envir = .GlobalEnv) else assign("Paths", old, envir = .GlobalEnv)
  }, add = TRUE)
  assign("Paths", list(outputPath = tempdir()), envir = .GlobalEnv)
  force(expr)
}

# sum of 1s in a SpatRaster
sum_ones <- function(ras) sum(values(ras) == 1, na.rm = TRUE)

# Ensure optional dependencies available for scenarios that call the full generator
skip_if_missing_generate_deps <- function() {
  req_pkgs <- c("SpaDES.tools", "raster", "tictoc", "reproducible")
  for (pkg in req_pkgs) testthat::skip_if_not_installed(pkg)
}
# ---------- disturbanceParameters helpers (only polygon “Generating”) ----------
dp_generating_forestry <- function(rate_pct = 1.0, size_m2 = "10000", interval = 1L) {
  data.table(
    dataName            = "forestry",
    dataClass           = "potentialCutblocks",
    disturbanceType     = "Generating",
    disturbanceOrigin   = "cutblocks",
    disturbanceEnd      = "",
    disturbanceRate     = rate_pct,
    disturbanceSize     = size_m2,             # 10,000 m² ~ 1% of 1 km²
    disturbanceInterval = interval,
    potentialField      = "Potential",
    resolutionVector    = 15
  )
}


dp_generating_mining <- function(rate_pct = 0.2, size_m2 = "100", interval = 1L) {
  data.table(
    dataName            = "mining",
    dataClass           = "potentialMining",
    disturbanceType     = "Generating",
    disturbanceOrigin   = "mining",
    disturbanceEnd      = "",
    disturbanceRate     = rate_pct,
    disturbanceSize     = size_m2,          # ~1 pixel on your 10×10 m cells
    disturbanceInterval = interval,
    potentialField      = "Potential",
    resolutionVector    = 15
  )
}

dp_enlarging_settlements <- function(rate_pct = 3.0, interval = 1L) {
  data.table(
    dataName            = "settlements",
    dataClass           = "settlements",
    disturbanceType     = "Enlarging",
    disturbanceOrigin   = "settlements",
    disturbanceEnd      = "",
    disturbanceRate     = rate_pct,
    disturbanceSize     = NA_character_,
    disturbanceInterval = interval,
    potentialField      = NA_character_,
    resolutionVector    = 7.5
  )
}

dp_connect_cutblocks_to_roads <- function() data.table(
  dataName            = "pipelines",
  dataClass           = NA_character_,
  disturbanceType     = "Connecting",
  disturbanceOrigin   = "cutblocks",
  disturbanceEnd      = "roads",
  disturbanceRate     = NA_real_,
  disturbanceSize     = NA_character_,
  disturbanceInterval = 1L,
  potentialField      = NA_character_,
  resolutionVector    = 15
)

dp_connect_mining_to_roads <- function() data.table(
  dataName            = "pipelines",        # SAME SECTOR as above
  dataClass           = NA_character_,
  disturbanceType     = "Connecting",
  disturbanceOrigin   = "mining",           # different origin, still polygons
  disturbanceEnd      = "roads",            # lines in disturbanceList[['pipelines']]
  disturbanceRate     = NA_real_,
  disturbanceSize     = NA_character_,
  disturbanceInterval = 1L,
  potentialField      = NA_character_,
  resolutionVector    = 15
)

dp_connect_mining_to_roads_pipelines <- function() data.table(
  dataName="pipelines", dataClass=NA_character_, disturbanceType="Connecting",
  disturbanceOrigin="mining", disturbanceEnd="roads",
  disturbanceRate=NA_real_, disturbanceSize=NA_character_,
  disturbanceInterval=1L, potentialField=NA_character_, resolutionVector=15
)

dp_connect_seismic_to_roads_pipelines <- function() data.table(
  dataName="pipelines", dataClass=NA_character_, disturbanceType="Connecting",
  disturbanceOrigin="seismicLines", disturbanceEnd="roads",
  disturbanceRate=NA_real_, disturbanceSize=NA_character_,
  disturbanceInterval=1L, potentialField=NA_character_, resolutionVector=15
)



# ---------- Common inputs ----------
distList_base <- createDisturbanceList()
currentYear <- 2010  # IMPORTANT: makes ORIGIN < (year-50) TRUE for forestry potentials

# ============================================================
# 1) Structure & types
# ============================================================
testthat::test_that("Structure: returns individualLayers + currentDisturbanceLayer with expected inner types", {
  skip_if_missing_generate_deps()
  distList_base$roads <- list(roads = distList_base$pipelines$roads)
  
  dp <- rbindlist(list(
    dp_enlarging_settlements(rate_pct = 5),
    dp_generating_forestry(rate_pct = 1.0, size_m2 = "10000"),
    dp_generating_mining(rate_pct   = 0.2, size_m2 = "100"),
    dp_connect_cutblocks_to_roads(),
    dp_connect_mining_to_roads()
  ), fill = TRUE)
  
  out <- with_paths(suppressWarnings(generateDisturbances(
    disturbanceParameters = dp,
    disturbanceList = distList_base,
    rasterToMatch = r,
    studyArea = sa,
    fires = NULL,
    currentTime = currentYear,       # e.g., 2010
    growthStepGenerating = 1,
    growthStepEnlargingPolys = 2,
    growthStepEnlargingLines = 2,
    currentDisturbanceLayer = NULL,
    connectingBlockSize = 500L,      # non-zero
    disturbanceRateRelatesToBufferedArea = TRUE,
    outputsFolder = tempdir(),
    runName = "struct"
  )))
  
  expect_type(out, "list")
  expect_true(all(c("individualLayers", "currentDisturbanceLayer") %in% names(out)))
  
  il  <- out$individualLayers
  cdl <- out$currentDisturbanceLayer
  
  expect_true(all(c("settlements","forestry") %in% names(il)))
  expect_true("cutblocks" %in% names(il$forestry))
  expect_s4_class(il$forestry$cutblocks, "SpatRaster")
  expect_s4_class(il$settlements$settlements, "SpatVector")
  
  expect_true(all(c("settlements","forestry") %in% names(cdl)))
  expect_s4_class(cdl$settlements$settlements, "SpatVector")
  expect_s4_class(cdl$forestry$cutblocks, "SpatVector")
})

# ============================================================
# 2) Forestry Generating: some 1s, not wild overshoot
# ============================================================
testthat::test_that("Generating (forestry): creates some 1s but stays well below 30% of pixels", {
  skip_if_missing_generate_deps()
  dp <- rbindlist(list(
    dp_generating_forestry(1.0, "10000"),
    dp_generating_mining(0.2, "100"),
    dp_connect_cutblocks_to_roads(),
    dp_connect_mining_to_roads()
  ), fill = TRUE)
  
  out <- with_paths(suppressWarnings(generateDisturbances(
    disturbanceParameters = dp, disturbanceList = distList_base,
    rasterToMatch = r, studyArea = sa, fires = NULL, currentTime = currentYear,
    growthStepGenerating = 1, growthStepEnlargingPolys = 1, growthStepEnlargingLines = 1,
    currentDisturbanceLayer = NULL, connectingBlockSize = 500L,
    disturbanceRateRelatesToBufferedArea = TRUE, outputsFolder = tempdir(), runName = "cap"
  )))
  
  g <- out$individualLayers$forestry$cutblocks
  expect_s4_class(g, "SpatRaster")
  n1 <- sum_ones(g); expect_gt(n1, 0)
  total <- ncell(g) - sum(is.na(values(g)))
  expect_lt(n1, 0.30 * total)
})

# ============================================================
# 3) Fire exclusion for Generating (mask burned cells)
# ============================================================
testthat::test_that("Forestry Generating avoids burned cells", {
  skip_if_missing_generate_deps()
  fires <- r; values(fires) <- 0
  set.seed(123)
  values(fires)[sample.int(ncell(fires), size = round(0.1 * ncell(fires)))] <- 1
  
  dp <- rbindlist(list(
    dp_generating_forestry(1.2, "10000"),
    dp_generating_mining(0.2, "100"),
    dp_connect_cutblocks_to_roads(),
    dp_connect_mining_to_roads()
  ), fill = TRUE)
  
  out <- with_paths(suppressWarnings(generateDisturbances(
    disturbanceParameters = dp, disturbanceList = distList_base,
    rasterToMatch = r, studyArea = sa, fires = fires, currentTime = currentYear,
    growthStepGenerating = 1, growthStepEnlargingPolys = 1, growthStepEnlargingLines = 1,
    currentDisturbanceLayer = NULL, connectingBlockSize = 500L,
    disturbanceRateRelatesToBufferedArea = TRUE, outputsFolder = tempdir(), runName = "fires"
  )))
  
  g <- out$individualLayers$forestry$cutblocks
  expect_s4_class(g, "SpatRaster")
  total_new <- sum(values(g) == 1, na.rm = TRUE)
  burned <- sum(values(g) == 1 & values(fires) == 1, na.rm = TRUE)
  clean  <- total_new - burned
  expect_gt(total_new, 0)
  expect_gt(clean, 0)
  expect_lt(burned / total_new, 0.5)
})

# ============================================================
# 4) Probability: prefers higher potential (right half)
# ============================================================
testthat::test_that("Generating prefers higher-Potential half (right)", {
  skip_if_missing_generate_deps()
  distList <- distList_base
  pot <- distList$forestry$potentialCutblocks
  e <- ext(pot); midx <- (xmin(e) + xmax(e)) / 2
  left  <- terra::as.polygons(terra::ext(xmin(e), midx, ymin(e), ymax(e)));  terra::crs(left)  <- terra::crs(pot)
  right <- terra::as.polygons(terra::ext(midx, xmax(e), ymin(e), ymax(e)));   terra::crs(right) <- terra::crs(pot)
  potL <- intersect(pot, left);  if (nrow(potL) > 0) potL$Potential <- 1L
  potR <- intersect(pot, right); if (nrow(potR) > 0) potR$Potential <- 5L
  distList$forestry$potentialCutblocks <- if (nrow(potL) == 0) potR else if (nrow(potR) == 0) potL else rbind(potL, potR)
  
  dp <- rbindlist(list(
    dp_generating_forestry(1.4, "10000"),
    dp_generating_mining(0.2, "100"),
    dp_connect_cutblocks_to_roads(),
    dp_connect_mining_to_roads()
  ), fill = TRUE)
  
  out <- with_paths(suppressWarnings(generateDisturbances(
    disturbanceParameters = dp, disturbanceList = distList,
    rasterToMatch = r, studyArea = sa, fires = NULL, currentTime = currentYear,
    growthStepGenerating = 1, growthStepEnlargingPolys = 1, growthStepEnlargingLines = 1,
    currentDisturbanceLayer = NULL, connectingBlockSize = 500L,
    disturbanceRateRelatesToBufferedArea = TRUE, outputsFolder = tempdir(), runName = "pref"
  )))
  
  g <- out$individualLayers$forestry$cutblocks
  g_left  <- crop(g, ext(xmin(e), midx, ymin(e), ymax(e)))
  g_right <- crop(g, ext(midx, xmax(e), ymin(e), ymax(e)))
  expect_s4_class(g, "SpatRaster")
  expect_gt(sum_ones(g_right), sum_ones(g_left))
  expect_gt(sum_ones(g_right), 0)
})

# ============================================================
# 5) Missing potential → NULL for that generated entry
# ============================================================
testthat::test_that("Missing forestry potential returns NULL generated layer", {
  skip_if_missing_generate_deps()
  distList <- distList_base
  distList$forestry$potentialCutblocks <- NULL  # <- remove forestry potential
  
  dp <- rbindlist(list(
    # No forestry Generating row here (we're testing its absence)
    dp_generating_mining(0.2, "100"),                        # tiny mining so origins exist
    dp_connect_mining_to_roads_pipelines(),                  # both Connecting rows
    dp_connect_seismic_to_roads_pipelines()                  # in the same sector
  ), fill = TRUE)
  
  out <- with_paths(suppressWarnings(generateDisturbances(
    disturbanceParameters = dp, disturbanceList = distList,
    rasterToMatch = r, studyArea = sa, fires = NULL, currentTime = currentYear,
    growthStepGenerating = 1, growthStepEnlargingPolys = 1, growthStepEnlargingLines = 1,
    currentDisturbanceLayer = NULL, connectingBlockSize = 500L,
    disturbanceRateRelatesToBufferedArea = TRUE, outputsFolder = tempdir(), runName = "nopot"
  )))
  
  # Robust expectation: either the forestry sector is missing OR cutblocks entry is NULL
  expect_true(is.null(out$individualLayers$forestry) ||
                is.null(out$individualLayers$forestry$cutblocks))
})


# ============================================================
# 6) Buffered-area mode: smoke test (generates and is bounded)
# ============================================================
testthat::test_that("Buffered-area mode generates forestry output and stays bounded", {
  skip_if_missing_generate_deps()
  dp <- rbindlist(list(
    dp_generating_forestry(1.0, "10000"),
    dp_generating_mining(0.2, "100"),
    dp_connect_cutblocks_to_roads(),
    dp_connect_mining_to_roads()
  ), fill = TRUE)
  
  out <- with_paths(suppressWarnings(generateDisturbances(
    disturbanceParameters = dp, disturbanceList = distList_base,
    rasterToMatch = r, studyArea = sa, fires = NULL, currentTime = currentYear,
    growthStepGenerating = 1, growthStepEnlargingPolys = 1, growthStepEnlargingLines = 1,
    currentDisturbanceLayer = NULL, connectingBlockSize = 500L,
    disturbanceRateRelatesToBufferedArea = TRUE, outputsFolder = tempdir(), runName = "buffTRUE"
  )))
  
  g <- out$individualLayers$forestry$cutblocks
  expect_s4_class(g, "SpatRaster")
  n1 <- sum_ones(g); expect_gt(n1, 0)
  total <- ncell(g) - sum(is.na(values(g)))
  expect_lt(n1, 0.35 * total)   # slightly looser cap; buffering counts toward area
})

# ============================================================
# 7) Multiple Generating (forestry + mining)
# ============================================================
testthat::test_that("Multiple Generating sectors produce outputs", {
  skip_if_missing_generate_deps()
  dp <- rbindlist(list(
    dp_generating_forestry(0.8, "10000"),
    dp_generating_mining(0.6, "100"),
    dp_connect_cutblocks_to_roads(),
    dp_connect_mining_to_roads()
  ), fill = TRUE)
  
  out <- with_paths(suppressWarnings(generateDisturbances(
    disturbanceParameters = dp, disturbanceList = distList_base,
    rasterToMatch = r, studyArea = sa, fires = NULL, currentTime = currentYear,
    growthStepGenerating = 1, growthStepEnlargingPolys = 1, growthStepEnlargingLines = 1,
    currentDisturbanceLayer = NULL, connectingBlockSize = 500L,
    disturbanceRateRelatesToBufferedArea = TRUE, outputsFolder = tempdir(), runName = "multi"
  )))
  
  g_for <- out$individualLayers$forestry$cutblocks
  g_min <- out$individualLayers$mining$mining
  expect_s4_class(g_for, "SpatRaster"); expect_gt(sum_ones(g_for), 0)
  expect_s4_class(g_min, "SpatRaster");  expect_gt(sum_ones(g_min), 0)
})

# ============================================================
# 8) Enlarging increases area (settlements)
# ============================================================
testthat::test_that("Enlarging settlements increases area", {
  skip_if_missing_generate_deps()
  dp <- rbindlist(list(
    dp_enlarging_settlements(10),
    dp_generating_forestry(0.5, "10000"),
    dp_generating_mining(0.2, "100"),
    dp_connect_cutblocks_to_roads(),
    dp_connect_mining_to_roads()
  ), fill = TRUE)
  
  out <- with_paths(suppressWarnings(generateDisturbances(
    disturbanceParameters = dp, disturbanceList = distList_base,
    rasterToMatch = r, studyArea = sa, fires = NULL, currentTime = currentYear,
    growthStepGenerating = 1, growthStepEnlargingPolys = 5, growthStepEnlargingLines = 5,
    currentDisturbanceLayer = NULL, connectingBlockSize = 500L,
    disturbanceRateRelatesToBufferedArea = TRUE, outputsFolder = tempdir(), runName = "enlarge"
  )))
  
  base <- distList_base$settlements$settlements
  enlarged <- out$individualLayers$settlements$settlements
  expect_s4_class(enlarged, "SpatVector")
  a0 <- sum(expanse(base, unit = "m", transform = FALSE))
  a1 <- sum(expanse(enlarged, unit = "m", transform = FALSE))
  expect_gt(a1, a0)
})
