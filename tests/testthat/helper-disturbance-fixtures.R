# tests/testthat/helper-disturbance-fixtures.R
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(data.table)
  library(testthat)
})

# --- study raster / study area: 10 km x 10 km (100 km²) ----------------------
r  <- rast(nrows = 200, ncols = 200, xmin = 0, xmax = 10000, ymin = 0, ymax = 10000, vals = 1)
crs(r) <- "EPSG:3005"
sa <- as.polygons(ext(r)); crs(sa) <- crs(r)

# Some module code expects this global
Paths <<- list(outputPath = tempdir())

# Resolution used by module for line/point buffering (see Enlarging logic).
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

#— 1) create the baseline disturbanceList -------------------------------------
createDisturbanceList <- function(crs_str = "EPSG:3005") {
  # settlements (seed ~1200x1200 m)
  poly1 <- st_polygon(list(rbind(c( 800,  800), c( 800, 2000), c(2000, 2000), c(2000,  800), c( 800,  800))))
  settlements          <- to_sv(st_sfc(poly1), "settlements",           crs_str)
  potentialSettlements <- to_sv(st_buffer(st_sfc(poly1), 1500), "potentialSettlements", crs_str)
  
  # wind (no existing; a few potential points)
  windTurbines          <- suppressWarnings(to_sv(st_sfc(), "windTurbines", crs_str))
  potentialWindTurbines <- to_sv(st_sfc(
    st_point(c(8500,8500)), st_point(c(9000,2500)), st_point(c(2500,9000)), st_point(c(6000,6000))
  ), "potentialWindTurbines", crs_str)
  
  # pipelines (+ roads)
  pipelines          <- suppressWarnings(to_sv(st_sfc(), "pipelines", crs_str))
  potentialPipelines <- to_sv(st_sfc(st_point(c(2000,8500)), st_point(c(4000,9000))), "potentialPipelines", crs_str)
  roads              <- to_sv(st_sfc(
    st_linestring(rbind(c(0,5000), c(10000,5000))),
    st_linestring(rbind(c(1000,1000), c(9000,9000)))
  ), "roads", crs_str)
  
  # mining (existing pad ~80,000 m²) + buffered potential
  mine_poly      <- st_polygon(list(rbind(c(500,100), c(900,100), c(900,300), c(500,300), c(500,100)))) # 400x200 m
  mining         <- to_sv(st_sfc(mine_poly), "mining", crs_str)
  potentialMining <- to_sv(st_buffer(st_sfc(mine_poly), 1500), "potentialMining", crs_str)
  
  # forestry: existing ~20,000 m² + two broad potential regions (Potential 1/2)
  cb_poly   <- st_polygon(list(rbind(c(8000,4000), c(8000,6000), c(9000,6000), c(9000,4000), c(8000,4000))))
  cutblocks <- to_sv(st_sfc(cb_poly), "cutblocks", crs_str)
  
  pot_f1 <- st_polygon(list(rbind(c(   0,   0), c(   0,10000), c(5000,10000), c(5000,   0), c(   0,   0))))
  pot_f2 <- st_polygon(list(rbind(c(5000,   0), c(5000,10000), c(10000,10000), c(10000,  0), c(5000,   0))))
  potentialCutblocks <- to_sv(st_sfc(pot_f1, pot_f2), "potentialCutblocks", crs_str)
  potentialCutblocks$Potential <- c(1L, 2L)
  potentialCutblocks <- erase(potentialCutblocks, cutblocks) # avoid overlap
  
  # oil & gas: several existing seismic lines + polygonal potentials (1/2)
  seis_lines <- st_sfc(
    st_linestring(rbind(c(3000,2000), c(3000,9500))),
    st_linestring(rbind(c(7000,1500), c(7000,9000))),
    st_linestring(rbind(c(1500,7000), c(8500,7000))) # cross line to vary angles
  )
  seismicLines <- to_sv(seis_lines, "seismicLines", crs_str)
  
  pot_seis_poly1 <- st_polygon(list(rbind(c(2000,6000), c(2000,10000), c(4000,10000), c(4000,6000), c(2000,6000))))
  pot_seis_poly2 <- st_polygon(list(rbind(c(1000,4500), c(1000,6500),  c(2000,6500),  c(2000,4500), c(1000,4500))))
  potentialSeismicLines <- to_sv(st_sfc(pot_seis_poly1, pot_seis_poly2), "potentialSeismicLines", crs_str)
  potentialSeismicLines$Potential <- c(1L, 2L) # polygons, as required
  
  list(
    settlements = list(settlements = settlements, potentialSettlements = potentialSettlements),
    wind        = list(windTurbines = windTurbines, potentialWindTurbines = potentialWindTurbines),
    pipelines   = list(pipelines = pipelines, potentialPipelines = potentialPipelines, roads = roads),
    mining      = list(mining = mining, potentialMining = potentialMining),
    forestry    = list(cutblocks = cutblocks, potentialCutblocks = potentialCutblocks),
    oilGas      = list(seismicLines = seismicLines, potentialSeismicLines = potentialSeismicLines)
  )
}

#— 2) disturbanceParameters (rates in PERCENT; sizes mostly inferred) ----------
# NOTE: The module converts percent to fraction internally (Rate <- disturbanceRate/100),
# so '2' means 2%, not 0.02. Keep this consistent across tests.
createDisturbanceParameters <- function(distList, res = .DEFAULT_RESOLUTION) {
  row_gen <- function(sector, origin, potential, rate_percent, size_expr = NA_character_, interval = 1L) {
    data.table(
      dataName            = sector,
      dataClass           = potential,
      disturbanceType     = "Generating",
      disturbanceRate     = rate_percent,   # percent, e.g., 0.2 = 0.2%
      disturbanceSize     = size_expr,      # NA => let calculateSize() infer from 'origin'
      disturbanceOrigin   = origin,
      disturbanceEnd      = "",
      disturbanceInterval = interval,
      resolutionVector    = res
    )
  }
  row_enl <- function(sector, origin, rate_percent, interval = 1L) {
    data.table(
      dataName            = sector,
      dataClass           = origin,
      disturbanceType     = "Enlarging",
      disturbanceRate     = rate_percent,
      disturbanceSize     = NA_character_,
      disturbanceOrigin   = origin,
      disturbanceEnd      = "",
      disturbanceInterval = interval,
      resolutionVector    = res
    )
  }
  
  # ---- Generating ----
  gen_rows <- rbindlist(list(
    # Forestry: target ~0.2% of 100 km² = 0.2 km² = 200,000 m² per step → ~10 x 20k m² blocks
    row_gen("forestry", "cutblocks",      "potentialCutblocks",     0.2, NA_character_),
    # Mining: ~0.1% = 100,000 m² per step → ~1 pad of ~80k m²
    row_gen("mining",   "mining",         "potentialMining",        0.1, NA_character_),
    # Seismic (area ~ of buffered lines sample): aim ~0.1% = 100,000 m²
    row_gen("oilGas",   "seismicLines",   "potentialSeismicLines",  0.1, "rtnorm(1, 30000, 9000, lower=0)"),
    # Wind: each turbine ~62,500 m² ⇒ 0.0625% of 100 km²; set 0.10% to allow ~1–2 turbines/step
    row_gen("wind",     "windTurbines",   "potentialWindTurbines",  0.10, "62500")
  ))
  
  # ---- Enlarging (settlements only; removed odd 'pipelines enlarging roads') ----
  enl_rows <- rbindlist(list(
    row_enl("settlements", "settlements", 0.5)  # 0.5% growth on (buffered) existing area
  ))
  
  rbindlist(list(gen_rows, enl_rows), fill = TRUE)
}

#— 3) updatedLayersAll that hits different branches ----------------------------
createUpdatedLayersAll <- function() {
  new_settlement_vect <- vect(st_sfc(st_polygon(list(rbind(
    c(1000,1000), c(1000,2200), c(2200,2200), c(2200,1000), c(1000,1000)
  ))))); crs(new_settlement_vect) <- crs(r)
  
  new_wind_vect <- vect(st_sfc(st_point(c(9000, 9000)))); crs(new_wind_vect) <- crs(r)
  
  new_dummy_vect <- vect(st_sfc(st_point(c(0, 0)))); crs(new_dummy_vect) <- crs(r)
  
  # A Raster for mining to force raster→vector conversion path
  new_mining_raster <- rast(nrows = 10, ncols = 10, xmin = 500, xmax = 900, ymin = 100, ymax = 300, vals = 1)
  crs(new_mining_raster) <- crs(r)
  
  # Different-geometry SpatVector for forestry (linestring)
  new_forestry_line <- vect(st_sfc(st_linestring(rbind(c(8000, 5000), c(9000, 5000)))))
  crs(new_forestry_line) <- crs(r)
  
  # Seismic lines + first year
  new_seis_vect <- vect(st_sfc(st_linestring(rbind(c(3000, 6500), c(3000, 9000)))))
  crs(new_seis_vect) <- crs(r)
  
  first_year_seis <- vect(st_sfc(st_linestring(rbind(c(2500, 5000), c(2500, 8000)))))
  crs(first_year_seis) <- crs(r)
  
  list(
    individualLayers  = list(
      settlements = list(settlements = new_settlement_vect),
      wind        = list(windTurbines = new_wind_vect),
      pipelines   = list(), # no updates → “no updates” path
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
  
  # Rates should be expressed in PERCENT because module divides by 100
  if (any(dpar$disturbanceRate > 100 | dpar$disturbanceRate < 0))
    stop("disturbanceRate must be in [0, 100] percent")
  
  invisible(TRUE)
}
