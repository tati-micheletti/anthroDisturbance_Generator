
module_name <- "anthroDisturbance_Generator"

find_modulePath <- function(mod = module_name) {
  # Start from the test dir and walk up to the module root
  module_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = TRUE)
  stopifnot(basename(module_root) == mod)
  # SpaDES expects modulePath to be the directory that CONTAINS the module folder
  dirname(module_root)
}

# ----- Visual diagnostics: inputs + buffer/erase previews ---------------------
plot_disturbance_diagnostics <- function(distList, sa, r, rstCurrentBurn,
                                         seismic_erase_width = 50,
                                         outdir = tempdir(), tag = "pre_run") {
  stopifnot(inherits(sa, "SpatVector") || inherits(sa, "SpatRaster"))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # helpers ----
  ac <- function(col, a = 0.35) grDevices::adjustcolor(col, alpha.f = a)
  hasv <- function(x) inherits(x, "SpatVector") && terra::nrow(x) > 0
  safe_valid <- function(v) try(suppressWarnings(terra::makeValid(v)), silent = TRUE)
  safe_erase <- function(pot, existing, buf = 0) {
    if (!hasv(pot)) return(pot)
    pot <- if (inherits(tmp <- safe_valid(pot), "try-error")) pot else tmp
    if (!hasv(existing)) return(pot)
    ex  <- if (buf > 0) terra::buffer(existing, width = buf) else existing
    ex  <- if (inherits(tmp <- safe_valid(ex), "try-error")) ex else tmp
    suppressWarnings(terra::erase(pot, ex))
  }
  plot_base <- function(main) {
    terra::plot(terra::as.polygons(terra::ext(r)), border = NA, col = "white")
    terra::plot(sa, col = ac("grey90", 1), border = "grey40", add = TRUE)
    title(main, cex.main = 1.0)
  }
  
  # derive availabilities (what the generator conceptually uses) ----
  settl      <- distList$settlements$settlements
  pot_settl  <- distList$settlements$potentialSettlements
  cut        <- distList$forestry$cutblocks
  pot_cut    <- distList$forestry$potentialCutblocks
  mine       <- distList$mining$mining
  pot_mine   <- distList$mining$potentialMining
  seis       <- distList$oilGas$seismicLines
  pot_seis   <- distList$oilGas$potentialSeismicLines
  roads      <- distList$roads$roads
  pipes      <- distList$pipelines$pipelines
  wind       <- distList$wind$windTurbines
  pot_wind   <- distList$wind$potentialWindTurbines
  
  avail_cut   <- safe_erase(pot_cut,  cut,  buf = 0)
  avail_mine  <- safe_erase(pot_mine, mine, buf = 0)
  seis_buf50  <- if (hasv(seis)) terra::buffer(seis, width = seismic_erase_width) else NULL
  avail_seis1 <- safe_erase(pot_seis, seis_buf50, buf = 0)  # first-year "availableArea" preview
  
  burned_v <- try({
    pol <- terra::as.polygons(rstCurrentBurn, values = TRUE, na.rm = TRUE)
    terra::subset(pol, pol[[1]] == 1)
  }, silent = TRUE)
  if (inherits(burned_v, "try-error")) burned_v <- NULL
  
  # file 1: overview of raw inputs --------------------------------------------
  f1 <- file.path(outdir, sprintf("diagnostics_%s_inputs.png", tag))
  grDevices::png(f1, width = 1600, height = 1200, res = 120)
  op <- par(mfrow = c(2, 3), mar = c(2, 2, 2, 1))
  
  plot_base("Study area + fire mask")
  if (!is.null(burned_v)) terra::plot(burned_v, col = ac("orange", .3), border = "orange", add = TRUE)
  
  plot_base("Settlements (existing + potential)")
  if (hasv(pot_settl)) terra::plot(pot_settl, col = ac("gold"), border = "gold3", add = TRUE)
  if (hasv(settl))     terra::plot(settl,     col = ac("red"),  border = "red3",  add = TRUE)
  
  plot_base("Forestry (cutblocks + potential)")
  if (hasv(pot_cut)) terra::plot(pot_cut, col = ac("darkseagreen2"), border = "darkseagreen4", add = TRUE)
  if (hasv(cut))     terra::plot(cut,     col = ac("darkgreen"),     border = "darkgreen",     add = TRUE)
  
  plot_base("Mining (existing + potential)")
  if (hasv(pot_mine)) terra::plot(pot_mine, col = ac("tan"),  border = "tan4", add = TRUE)
  if (hasv(mine))     terra::plot(mine,     col = ac("sienna"), border = "sienna4", add = TRUE)
  
  plot_base("Oil&Gas (seismic lines + potential)")
  if (hasv(pot_seis)) terra::plot(pot_seis, col = ac("lightblue"), border = "steelblue3", add = TRUE)
  if (hasv(seis))     terra::plot(seis,     col = "black", lwd = 1.5, add = TRUE)
  
  plot_base("Wind (turbines + potential)")
  if (hasv(pot_wind)) terra::plot(pot_wind, col = ac("olivedrab2"), border = "olivedrab4", add = TRUE)
  if (hasv(wind))     terra::points(wind,   pch = 16, cex = 0.6, col = "black")
  
  par(op); grDevices::dev.off()
  
  # file 2: buffer/erase previews ----------------------------------------------
  f2 <- file.path(outdir, sprintf("diagnostics_%s_ops.png", tag))
  grDevices::png(f2, width = 1600, height = 1200, res = 120)
  op <- par(mfrow = c(2, 3), mar = c(2, 2, 2, 1))
  
  plot_base("Forestry availability = potential - cutblocks")
  if (hasv(pot_cut))  terra::plot(pot_cut,  col = ac("darkseagreen2", .25), border = NA, add = TRUE)
  if (hasv(cut))      terra::plot(cut,      col = ac("darkgreen", .6),       border = "darkgreen", add = TRUE)
  if (hasv(avail_cut))terra::plot(avail_cut,col = ac("chartreuse3", .4),      border = "chartreuse4", add = TRUE)
  
  plot_base("Mining availability = potential - mines")
  if (hasv(pot_mine))   terra::plot(pot_mine,   col = ac("tan", .25), border = NA, add = TRUE)
  if (hasv(mine))       terra::plot(mine,       col = ac("sienna", .6), border = "sienna4", add = TRUE)
  if (hasv(avail_mine)) terra::plot(avail_mine, col = ac("navajowhite3", .5), border = "peru", add = TRUE)
  
  plot_base(sprintf("Seismic: potential – %dm buffer(lines)", seismic_erase_width))
  if (hasv(pot_seis))  terra::plot(pot_seis,  col = ac("lightblue", .25),  border = NA, add = TRUE)
  if (hasv(seis_buf50))terra::plot(seis_buf50,col = ac("purple", .35), border = "purple4", add = TRUE)
  if (hasv(avail_seis1))terra::plot(avail_seis1, col = ac("turquoise3", .5), border = "turquoise4", add = TRUE)
  
  plot_base("Roads & pipelines (context)")
  if (hasv(roads)) terra::plot(roads, col = "grey40", lwd = 1, add = TRUE)
  if (hasv(pipes)) terra::plot(pipes, col = "grey20", lwd = 1, add = TRUE)
  
  plot_base("Fire mask (raster polygons)")
  if (!is.null(burned_v)) terra::plot(burned_v, col = ac("orange", .4), border = "orange3", add = TRUE)
  
  plot_base("Settlements (no erase preview; enlarging in-sim)")
  if (hasv(pot_settl)) terra::plot(pot_settl, col = ac("gold", .3), border = "gold4", add = TRUE)
  if (hasv(settl))     terra::plot(settl,     col = ac("red", .5),  border = "red4",  add = TRUE)
  
  par(op); grDevices::dev.off()
  
  message("Wrote: ", f1, "\nWrote: ", f2)
  invisible(c(f1, f2))
}


testthat::test_that("Primary sim runs one year and updates outputs (vector mode)", {
  testthat::skip_if_not_installed("SpaDES.core")
  testthat::skip_if_not_installed("SpaDES.tools")
  testthat::skip_if_not_installed("reproducible")
  testthat::skip_if_not_installed("terra")
  testthat::skip_if_not_installed("data.table")
  testthat::skip_if_not_installed("sf")
  testthat::skip_if_not_installed("fasterize")
  testthat::skip_if_not_installed("raster")
  testthat::skip_if_not_installed("tictoc")
  
  # Attach only what the module calls unqualified
  attach_pkg_safely <- function(pkg) {
    if (paste0("package:", pkg) %in% search()) return(invisible(TRUE))
    ok <- try(suppressPackageStartupMessages(library(pkg, character.only = TRUE, quietly = TRUE)), silent = TRUE)
    if (!inherits(ok, "try-error")) return(invisible(TRUE))
    loadNamespace(pkg); attachNamespace(pkg); invisible(TRUE)
  }
  attach_pkg_safely("reproducible")
  attach_pkg_safely("SpaDES.tools")
  
  withr::local_options(
    spades.useRequire       = FALSE,
    spades.moduleCodeChecks = FALSE,
    Require.install         = FALSE,
    Require.unload          = FALSE
  )
  
  # ---- Fixtures from helpers -------------------------------------------------
  # helpers define: r (10 km x 10 km EPSG:3005), sa, Paths, .DEFAULT_RESOLUTION
  # build a rich disturbanceList, then prune to just sectors we want
  distList_full <- createDisturbanceList(crs(r))     # helpers file
  dp_full       <- createDisturbanceParameters(distList_full)  # helpers file
  
  # keep only: settlements (Enlarging) + forestry (Generating)
  keep_rows <- dp_full$dataName %in% c("settlements","forestry","oilGas")
  dp <- data.table::copy(dp_full[keep_rows])
  dp[dataName == "oilGas" & disturbanceType == "Generating",
     disturbanceRate := 0] 
  
  # force a simple, robust forestry size (avoid truncnorm paths)
  # ~20k m² blocks so 0.2% of 100 km² (~200k m²) is reachable in ~10 blocks
  dp[dataName == "forestry" & disturbanceType == "Generating",
     disturbanceSize := "max(stats::rnorm(1, mean = 20000, sd = 3000), 1)"]
  
  # sanity-check helper fixtures
  check_fixtures(distList_full, dp)                  # helpers file
  
  # prune the list to just what we need (settlements, forestry)
  disturbanceList <- list(
    settlements = distList_full$settlements[c("settlements","potentialSettlements")],
    forestry    = distList_full$forestry[c("cutblocks","potentialCutblocks")],
    oilGas      = distList_full$oilGas[c("seismicLines","potentialSeismicLines")]
  )
  
  # minimal DisturbanceRate object (unused if rates provided, but harmless)
  DisturbanceRate <- data.table::data.table(
    dataName          = c("settlements","forestry"),
    dataClass         = c("settlements","potentialCutblocks"),
    disturbanceType   = c("Enlarging", "Generating"),
    disturbanceOrigin = c("settlements","cutblocks"),
    disturbanceRate   = c(0.05, 0.05)
  )
  
  # seeds for reproducibility of module events (names per Appendix C events)
  paramsList <- list(
    anthroDisturbance_Generator = list(
      generatedDisturbanceAsRaster       = FALSE,  # vector mode
      saveInitialDisturbances            = FALSE,
      saveCurrentDisturbances            = FALSE,
      runInterval                        = 1,
      growthStepEnlargingPolys           = 5,
      growthStepEnlargingLines           = 5,
      growthStepGenerating               = 1,
      disturbanceRateRelatesToBufferedArea = FALSE,
      seismicLineGrids                   = 1,
      useClusterMethod                   = FALSE,
      runClusteringInParallel            = FALSE,
      useRoadsPackage                    = FALSE,
      maskWaterAndMountainsFromLines     = FALSE,
      siteSelectionAsDistributing        = NA_character_,
      connectingBlockSize                = 50L,
      outputsFolder                      = tempdir(),
      runName                            = "unit",
      featuresToAvoid                    = NULL,
      altitudeCut                        = Inf,
      clusterDistance                    = 10,
      distanceNewLinesFactor             = 1,
      refinedStructure                   = FALSE,
      disturbFirstYear                   = FALSE,
      .seed = list(anthroDisturbance_Generator = list(
        calculatingSize=1L, calculatingRate=2L, generatingDisturbances=3L, updatingDisturbanceList=4L
      ))
    )
  )
  
  DEM <- r; terra::values(DEM) <- 300
  rstCurrentBurn <- r; terra::values(rstCurrentBurn) <- 0
  
  modPath <- find_modulePath()
  
  sim <- SpaDES.core::simInit(
    times   = list(start = 2000, end = 2001),
    params  = paramsList,
    modules = "anthroDisturbance_Generator",
    paths   = list(modulePath = modPath, inputPath = tempdir(), outputPath = tempdir()),
    objects = list(
      studyArea             = sa,
      rasterToMatch         = r,
      DEM                   = DEM,
      rstCurrentBurn        = rstCurrentBurn,
      DisturbanceRate       = DisturbanceRate,
      disturbanceList       = disturbanceList,
      disturbanceParameters = dp,
      disturbanceDT         = data.table::data.table(dataName="forestry", dataClass="potentialCutblocks")
    ),
    loadPkgs = FALSE
  )
  
  simOut <- SpaDES.core::spades(sim, until = 2)
  
  # ---- Assertions: year layer exists & settlements grew ----------------------
  testthat::expect_true("Year2001" %in% names(simOut$currentDisturbanceLayer))
  lay <- simOut$currentDisturbanceLayer[["Year2001"]]
  testthat::expect_true(is.list(lay))  # vector mode returns a list of SpatVectors
  
  # settlements enlarging increased area
  a0 <- sum(terra::expanse(disturbanceList$settlements$settlements, unit = "m"))
  a1 <- sum(terra::expanse(simOut$disturbanceList$settlements$settlements, unit = "m"))
  testthat::expect_true(a1 >= a0)
  
  # forestry generating did not shrink total area (should add polygons)
  f0 <- sum(terra::expanse(disturbanceList$forestry$cutblocks, unit = "m"))
  f1 <- sum(terra::expanse(simOut$disturbanceList$forestry$cutblocks, unit = "m"))
  testthat::expect_true(f1 >= f0)
})
