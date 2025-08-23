# tests/testthat/helper-disturbance-fixtures.R
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(data.table)
  library(testthat)
})

# --- shared raster / study area ------------------------------------------------
r  <- rast(nrows = 100, ncols = 100, xmin = 0, xmax = 1000, ymin = 0, ymax = 1000, vals = 1)
crs(r) <- "EPSG:3005"
sa <- as.polygons(ext(r))
crs(sa) <- crs(r)

# Some module code expects this global
Paths <<- list(outputPath = tempdir())

# Resolution used by module for line/point buffering.
# Keep scalar, the generator expects a single unique numeric (15 m is a safe default).
.DEFAULT_RESOLUTION <- 15

# Wrap sf -> SpatVector and attach "Class"; for *potential* layers add defaults if missing
to_sv <- function(sf_geom, class_name, crs_str) {
  v <- vect(sf_geom); crs(v) <- crs_str
  if (nrow(v) > 0) {
    v$Class <- class_name
    if (startsWith(class_name, "potential")) {
      if (!"Potential" %in% names(v)) v$Potential <- 1L
      if (!"ORIGIN"    %in% names(v)) v$ORIGIN    <- 1900L
    }
  } else {
    names(v) <- character(0)
  }
  v
}

#— 1) create the *baseline* disturbanceList, including a dummy sectorX -------
createDisturbanceList <- function(crs_str = "EPSG:3005") {
  # settlements (tiny seed)
  poly1 <- st_polygon(list(rbind(c( 50,  50), c( 50, 150), c(150, 150), c(150,  50), c( 50,  50))))
  settlements          <- to_sv(st_sfc(poly1), "settlements",          crs_str)
  potentialSettlements <- to_sv(st_buffer(st_sfc(poly1), 150), "potentialSettlements", crs_str)
  
  # wind (no existing, a couple of potential points)
  windTurbines          <- suppressWarnings(to_sv(st_sfc(), "windTurbines", crs_str))
  potentialWindTurbines <- to_sv(st_sfc(st_multipoint(rbind(c(800,800), c(850,250)))), "potentialWindTurbines", crs_str)
  
  # pipelines (+ roads)
  pipelines          <- suppressWarnings(to_sv(st_sfc(), "pipelines", crs_str))
  potentialPipelines <- to_sv(st_sfc(st_point(c(200,850)), st_point(c(400,900))), "potentialPipelines", crs_str)
  roads              <- to_sv(st_sfc(st_linestring(rbind(c(0,500), c(1000,500)))), "roads", crs_str)
  
  # mining (one existing pad + buffered potential)
  mine_poly        <- st_polygon(list(rbind(c(500,100), c(500,200), c(600,200), c(600,100), c(500,100))))
  mining           <- to_sv(st_sfc(mine_poly), "mining", crs_str)
  potentialMining  <- to_sv(st_buffer(st_sfc(mine_poly), 150), "potentialMining", crs_str) # has Potential by to_sv()
  
  # forestry: small existing cutblock + two broad potential regions (Potential 1/2)
  cb_poly   <- st_polygon(list(rbind(c(800,400), c(800,600), c(900,600), c(900,400), c(800,400))))
  cutblocks <- to_sv(st_sfc(cb_poly), "cutblocks", crs_str)
  
  pot_f1 <- st_polygon(list(rbind(c(0,0), c(0,1000), c(500,1000), c(500,0), c(0,0))))
  pot_f2 <- st_polygon(list(rbind(c(500,0), c(500,1000), c(1000,1000), c(1000,0), c(500,0))))
  potentialCutblocks <- to_sv(st_sfc(pot_f1, pot_f2), "potentialCutblocks", crs_str)
  potentialCutblocks$Potential <- c(1L, 2L)
  # avoid overlap with existing cutblock
  potentialCutblocks <- erase(potentialCutblocks, cutblocks)
  
  # oil & gas: one seismic line + polygonal potential (Potential 1/2)
  seis_line <- st_linestring(rbind(c(300,500), c(300,1000)))
  seismicLines <- to_sv(st_sfc(seis_line), "seismicLines", crs_str)
  
  pot_seis_poly1 <- st_polygon(list(rbind(c(200,600), c(200,1000), c(400,1000), c(400,600), c(200,600))))
  pot_seis_poly2 <- st_polygon(list(rbind(c(100,450), c(100,650),  c(200,650),  c(200,450), c(100,450))))
  potentialSeismicLines <- to_sv(st_sfc(pot_seis_poly1, pot_seis_poly2), "potentialSeismicLines", crs_str)
  potentialSeismicLines$Potential <- c(1L, 2L)
  
  list(
    settlements = list(settlements = settlements, potentialSettlements = potentialSettlements),
    wind        = list(windTurbines = windTurbines, potentialWindTurbines = potentialWindTurbines),
    pipelines   = list(pipelines = pipelines, potentialPipelines = potentialPipelines, roads = roads),
    mining      = list(mining = mining, potentialMining = potentialMining),
    forestry    = list(cutblocks = cutblocks, potentialCutblocks = potentialCutblocks),
    oilGas      = list(seismicLines = seismicLines, potentialSeismicLines = potentialSeismicLines)
  )
}

#— 2) disturbanceParameters with realistic rates/sizes and correct wiring -----
# Notes:
# - For Generating rows: dataClass = potential layer, disturbanceOrigin = existing/current layer
# - For Enlarging rows:  dataClass unused; set = origin for clarity
# - Sizes chosen so they exceed minimal effective buffer widths and yield plausible lengths/areas.
createDisturbanceParameters <- function(distList, res = .DEFAULT_RESOLUTION) {
  row_gen <- function(sector, origin, potential, rate, size_expr, interval = 1L) {
    data.table(
      dataName            = sector,
      dataClass           = potential,      # potential layer used for siting
      disturbanceType     = "Generating",
      disturbanceRate     = rate,           # proportion of total STUDY AREA per time step (e.g., 0.02 = 2%)
      disturbanceSize     = size_expr,      # string expression, evaluated later
      disturbanceOrigin   = origin,         # existing/current layer (excluded when siting)
      disturbanceEnd      = "",
      disturbanceInterval = interval,
      resolutionVector    = res
    )
  }
  row_enl <- function(sector, origin, rate, interval = 1L) {
    data.table(
      dataName            = sector,
      dataClass           = origin,
      disturbanceType     = "Enlarging",
      disturbanceRate     = rate,
      disturbanceSize     = NA_character_,
      disturbanceOrigin   = origin,
      disturbanceEnd      = "",
      disturbanceInterval = interval,
      resolutionVector    = res
    )
  }
  
  # ---- Generating (plausible, buffer-safe sizes) ----
  # Forestry cutblocks: ~10–30 ha; mean 20 ha (200,000 m²), sd 6 ha
  # Mining pads: ~3–12 ha; mean 8 ha (80,000 m²), sd 2 ha
  # Seismic (as area of buffered line): choose mean 30,000 m², sd 9,000 m²
  # Wind turbines: keep module default pixel-sized pad (62,500 m²) when unknown
  gen_rows <- rbindlist(list(
    row_gen("forestry", "cutblocks",      "potentialCutblocks",     0.020, "rtnorm(1, 200000, 60000,  lower=0)"),
    row_gen("mining",   "mining",         "potentialMining",        0.005, "rtnorm(1,  80000, 20000,  lower=0)"),
    row_gen("oilGas",   "seismicLines",   "potentialSeismicLines",  0.010, "rtnorm(1,  30000,  9000,  lower=0)"),
    row_gen("wind",     "windTurbines",   "potentialWindTurbines",  0.001, "62500") # module default for turbines
  ))
  
  # ---- Enlarging ----
  # Settlements: gentle growth (0.5% of existing buffered area per step)
  enl_rows <- rbindlist(list(
    row_enl("settlements", "settlements", 0.005),
    row_enl("pipelines",   "roads",       0.000) 
  ))
  
  rbindlist(list(gen_rows, enl_rows), fill = TRUE)
}

#— 3) updatedLayersAll that actually hits multiple branches --------------------
createUpdatedLayersAll <- function() {
  # Enlarging vectors
  new_settlement_vect <- vect(st_sfc(st_polygon(list(rbind(
    c( 80,  80), c( 80, 200), c(200, 200), c(200,  80), c( 80,  80)
  ))))); crs(new_settlement_vect) <- crs(r)
  
  new_wind_vect <- vect(st_sfc(st_point(c(850, 850)))); crs(new_wind_vect) <- crs(r)
  
  new_dummy_vect <- vect(st_sfc(st_point(c(0, 0)))); crs(new_dummy_vect) <- crs(r)
  
  # A Raster for mining to force raster→vector conversion path
  new_mining_raster <- rast(nrows = 5, ncols = 5, xmin = 500, xmax = 600, ymin = 100, ymax = 200, vals = 1)
  crs(new_mining_raster) <- crs(r)
  
  # Different-geometry SpatVector for forestry (linestring)
  new_forestry_line <- vect(st_sfc(st_linestring(rbind(c(800, 500), c(900, 500)))))
  crs(new_forestry_line) <- crs(r)
  
  # Seismic lines + first year
  new_seis_vect <- vect(st_sfc(st_linestring(rbind(c(300, 600), c(300, 900)))))
  crs(new_seis_vect) <- crs(r)
  
  first_year_seis <- vect(st_sfc(st_linestring(rbind(c(250, 500), c(250, 800)))))
  crs(first_year_seis) <- crs(r)
  
  list(
    individualLayers  = list(
      settlements = list(settlements = new_settlement_vect),
      wind        = list(windTurbines = new_wind_vect),
      pipelines   = list(), # no updates → tests “no updates” path
      mining      = list(mining = new_mining_raster),
      forestry    = list(cutblocks = new_forestry_line),
      oilGas      = list(seismicLines = new_seis_vect),
      sectorX     = list(layerX = new_dummy_vect)
    ),
    seismicLinesFirstYear = first_year_seis
  )
}

#— 4) quick sanity checks to fail fast during tests ----------------------------
check_fixtures <- function(distList, dpar) {
  stopifnot(is.list(distList), inherits(dpar, "data.table"))
  # Generating rows must map to existing + potential layers
  gen <- dpar[disturbanceType == "Generating"]
  for (i in seq_len(nrow(gen))) {
    s   <- gen$dataName[i]
    pot <- gen$dataClass[i]
    org <- gen$disturbanceOrigin[i]
    if (is.null(distList[[s]][[pot]])) stop(sprintf("Missing potential layer %s/%s", s, pot))
    if (is.null(distList[[s]][[org]])) stop(sprintf("Missing origin layer %s/%s", s, org))
  }
  # Potential layers must carry 'Potential'; seismic potentials must be polygons
  for (s in names(distList)) for (lay in names(distList[[s]])) {
    v <- distList[[s]][[lay]]
    if (startsWith(lay, "potential")) {
      if (!"Potential" %in% names(v)) stop(sprintf("No 'Potential' in %s/%s", s, lay))
      if (lay == "potentialSeismicLines" && geomtype(v) != "polygons")
        stop("potentialSeismicLines must be polygons")
    }
  }
  # Resolution numeric singleton; rates present
  if (length(unique(dpar$resolutionVector)) != 1L || !is.numeric(unique(dpar$resolutionVector)))
    stop("resolutionVector must be a single numeric value")
  if (any(is.na(dpar[disturbanceType %in% c("Generating","Enlarging"), disturbanceRate])))
    stop("Provide non-NA disturbanceRate for Generating/Enlarging")
  invisible(TRUE)
}
