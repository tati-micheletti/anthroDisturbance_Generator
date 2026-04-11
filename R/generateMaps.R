if (FALSE){
  library("Require")
  Require("raster")
  Require("reproducible")
  Require("terra")
  Require("data.table")
  Require("sf")
  Require("tictoc")
  Require("rgdal")
  
  grepMulti <- function(x, patterns, unwanted = NULL) {
    rescued <- sapply(x, function(fun) all(sapply(X = patterns, FUN = grepl, fun)))
    recovered <- x[rescued]
    if (!is.null(unwanted)){
      discard <- sapply(recovered, function(fun) all(sapply(X = unwanted, FUN = grepl, fun)))
      afterFiltering <- recovered[!discard]
      return(afterFiltering)
    } else {
      return(recovered)
    }
  }
  bufferCells <- function(ras, valToLookFor, newValue){
    rDT <- data.table(pixelID = 1:ncell(ras),
                      vals = as.numeric(values(ras)))
    whichCells <- rDT[vals == valToLookFor, pixelID]
    ADJ <- unique(na.omit(as.numeric(adjacent(ras, cells = whichCells, 
                                              directions = "queen", 
                                              include = TRUE))))
    ras[ADJ] <- newValue
    return(ras)
  }
  
  
  ## SHAPEFILES
  # POWERLINES
  allForRas <- grepMulti(x = list.files(path = file.path(getwd(), "outputs"), 
                                        full.names = TRUE), 
                         patterns = c("disturbances_Energy_powerLines", ".shp"), 
                         unwanted = "_IC")
  powerlines <- lapply(allForRas, shapefile)
  
  # plot(powerlines[[5]], col = "blue"); plot(powerlines[[4]], col = "red", add = TRUE); plot(powerlines[[3]], col = "orange", add = TRUE);plot(powerlines[[2]], col = "purple", add = TRUE);plot(powerlines[[1]], col = "navyblue", add = TRUE)
  
  # PIPELINE
  allForRas <- grepMulti(x = list.files(path = file.path(getwd(), "outputs"), 
                                        full.names = TRUE), 
                         patterns = c("disturbances_oilGas_pipeline", ".shp"), 
                         unwanted = "_IC")
  pipelines <- lapply(allForRas, shapefile)
  
  # plot(pipelines[[5]], col = "blue"); plot(pipelines[[4]], col = "red", add = TRUE); plot(pipelines[[3]], col = "orange", add = TRUE);plot(pipelines[[2]], col = "purple", add = TRUE); plot(pipelines[[1]], col = "navyblue", add = TRUE)
  
  # SEISMIC LINES
  allForRas <- grepMulti(x = list.files(path = file.path(getwd(), "outputs"), 
                                        full.names = TRUE), 
                         patterns = c("disturbances_oilGas_seismicLines", ".shp"), 
                         unwanted = "_IC")
  seismic <- lapply(allForRas, shapefile)
  
  # plot(seismic[[5]], col = "blue"); plot(seismic[[4]], col = "red", add = TRUE); plot(seismic[[3]], col = "orange", add = TRUE);plot(seismic[[2]], col = "purple", add = TRUE);plot(seismic[[1]], col = "navyblue", add = TRUE)
  
  
  # ROADS
  allForRas <- grepMulti(x = list.files(path = file.path(getwd(), "outputs"), 
                                        full.names = TRUE), 
                         patterns = c("disturbances_roads_roads", ".shp"), 
                         unwanted = "_IC")
  roads <- lapply(allForRas, shapefile)
  
  # plot(roads[[5]], col = "blue"); plot(roads[[4]], col = "red", add = TRUE); plot(roads[[3]], col = "orange", add = TRUE);plot(roads[[2]], col = "purple", add = TRUE);plot(roads[[1]], col = "navyblue", add = TRUE)
  
  # SETTLEMENTS
  allForRas <- grepMulti(x = list.files(path = file.path(getwd(), "outputs"), 
                                        full.names = TRUE), 
                         patterns = c("disturbances_settlements_settlements", ".shp"), 
                         unwanted = "IC")
  sett <- lapply(allForRas, shapefile)
  
  # Other lines: disturbances_settlements_otherLines_IC
  allForRas <- grepMulti(x = list.files(path = file.path(getwd(), "outputs"), 
                                        full.names = TRUE), 
                         patterns = c("disturbances_settlements_otherLines_IC", ".shp"))
  otherLines <- lapply(allForRas, shapefile)
  
  # Other polygons: disturbances_settlements_otherPolygons_IC
  allForRas <- grepMulti(x = list.files(path = file.path(getwd(), "outputs"), 
                                        full.names = TRUE), 
                         patterns = c("disturbances_settlements_otherPolygons_IC", ".shp"))
  otherPolys <- lapply(allForRas, shapefile)
  
  ## RASTERS
  # WINDTURBINE
  allForRas <- grepMulti(x = list.files(path = file.path(getwd(), "outputs"), 
                                        full.names = TRUE), 
                         patterns = c("disturbances_Energy_windTurbines", ".tif"), 
                         unwanted = ".xml")
  wind <- stack(lapply(allForRas, raster))
  # FORESTRY
  allForRas <- grepMulti(x = list.files(path = file.path(getwd(), "outputs"), 
                                        full.names = TRUE), 
                         patterns = c("disturbances_forestry", ".tif"), 
                         unwanted = ".xml")
  stk <- stack(lapply(allForRas, raster))
  forestry <- stack(lapply(allForRas, raster))
  # MINING
  allForRas <- grepMulti(x = list.files(path = file.path(getwd(), "outputs"), 
                                        full.names = TRUE), 
                         patterns = c("disturbances_mining_mining", ".tif"), 
                         unwanted = ".xml")
  mines <- stack(lapply(allForRas, raster))
  # OILGAS
  allForRas <- grepMulti(x = list.files(path = file.path(getwd(), "outputs"), 
                                        full.names = TRUE), 
                         patterns = c("disturbances_oilGas_oilGas", ".tif"), 
                         unwanted = ".xml")
  oil <- stack(lapply(allForRas, raster))
  
  # LOAD STUDY AREAS
  
  SA <- prepInputs(url = "https://drive.google.com/open?id=1fYvNPwovjNtTABoGcegrvdFGkNfCUsxf", 
                   destinationPath = file.path(getwd(), "outputs"))
  SA <- vect(SA)
  SArep <- project(SA, crs(oil))
  Ede <- buffer(SArep, width = 10^5)
  
  Boo <- prepInputs(url = "https://drive.google.com/file/d/1x_fQEKHW2nGbqo1JvCpDwmVuTPYAavl3/view?usp=sharing", 
                    destinationPath = file.path(getwd(), "outputs"))
  Boo <- vect(Boo)
  Boo <- project(Boo, crs(oil))
  
  Herds <- prepInputs(url = "https://drive.google.com/file/d/18xEvO5ktCSRLRZC3KrsymTBMh_4MBB7Z/view?usp=sharing", 
                      destinationPath = file.path(getwd(), "outputs"))
  Herds <- vect(Herds)
  Herds <- project(Herds, crs(oil))
  Herds <- crop(Herds, Boo)
  Herds <- subset(Herds, Herds$PROV_TERR %in% c("NWT"))
  
  ## PRELIMINARY RESULTS
  # Rasters: oil, mines, forestry, wind
  # Shapefiles: sett, roads, seismic, pipelines, powerlines
  # NOT DIRECTLY USED: otherPolys, otherLines
  
  rasNames <- c("oil", "mines", "forestry", "wind")
  shapesNames <- c("sett", "roads", "seismic", "pipelines", "powerlines")
  
  # Cute NWT Map
  par(bg = NA)
  raster::plot(Boo, 
               col = "#072C62",
               axes = FALSE, 
               box = FALSE)
  legend(x = 0.7*xmax(Boo),
         y = ymax(Boo),
         legend = "Northwest\nTerritories",
         fill = "#072C62",
         cex = 0.95,
         box.lwd = 0)
  dev.copy(png,'NT_1.png')
  dev.off()
  
  par(bg = NA)
  raster::plot(Boo, 
               col = viridisLite::viridis(NROW(Boo)),
               axes = FALSE, 
               box = FALSE)
  legend(x = xmax(Boo), y = ymax(Boo), 
         legend = Boo$REGION,
         fill = viridisLite::viridis(NROW(Boo)),
         cex = 0.95,
         box.lwd = 0)
  dev.copy(png,'NT_2.png')
  dev.off()
  
  par(bg = NA)
  
  Box <- aggregate(Boo)
  raster::plot(Box, 
               col = "white",
               axes = FALSE, 
               box = FALSE)
  raster::plot(Herds, 
               add = TRUE,
               col = viridisLite::viridis(NROW(Herds)),
               axes = FALSE, 
               box = FALSE)
  dev.copy(png,'NT_4.png')
  dev.off()
  
  par(bg = NA)
  raster::plot(Boo, 
               col = viridisLite::viridis(NROW(Boo)),
               axes = FALSE, 
               box = FALSE)
  legend(x = xmax(Boo), y = ymax(Boo),
         legend = Boo$REGION,
         fill = viridisLite::viridis(NROW(Boo)),
         cex = 0.95,
         box.lwd = 0)
  dev.copy(png,'NT_6.png')
  dev.off()
  
  raster::plot(Boo, 
               col = viridisLite::viridis(NROW(Boo)),
               axes = FALSE, 
               box = FALSE)
  dev.copy(png,'NT_6.png')
  dev.off()
  
  par(bg = "white")
  
  # 1. Current rate of change of individual sectors over the study area
  # Calculate how many 1's we have in each map
  # Calculate how many 0's we have in each map
  # Do 1 over 0 to calculate the rate of change for each of the sectors
  
  rateOfChange <- rbindlist(lapply(c(rasNames, shapesNames), function(TYPE){ 
    tic(paste0("Time elapsed for ", TYPE, ":"))
    dist <- get(TYPE)
    if (class(dist) == "RasterStack"){
      dist[] <- dist[]
      tb <- rbindlist(lapply(1:nlayers(dist), function(N){
        TB <- data.table(table(dist[[N]][]))
        TB2 <- data.table(Zero = TB$N[1],
                          One = TB$N[2])
        TB2[, Year := 2001+(N*10)]
        return(TB2)
      }))
      tb[, total := Zero + One]
      tb[, percArea := One/total]
    } else {
      tb <- rbindlist(lapply(1:length(dist), function(N){
        TB <- data.table(Zero = sum(expanse(Boo)), 
                         One = sum(expanse(vect(dist[[N]]))))
        TB[, Year := 2001+(N*10)]
        return(TB)
      }))
      tb[, percArea := One/Zero]
    }
    changePerc <- numeric(NROW(tb))
    for (i in 1:NROW(tb)){
      changePerc[i+1] <- (tb$percArea[(i+1)]-tb$percArea[(i)])/tb$percArea[(i)]
    }
    tb[, rateChange := na.omit(changePerc)]
    changePercFirst <- numeric(NROW(tb))
    for (i in 1:NROW(tb)){
      changePercFirst[i+1] <- (tb$percArea[(i+1)]-tb$percArea[(1)])/tb$percArea[(1)]
    }
    tb[, totalRateChange := na.omit(changePercFirst)]
    tb[, sector := TYPE]
    DT <- tb[, c("sector", "Year", "percArea", "rateChange", "totalRateChange")]
    toc()
    return(DT)
  }))
  options(scipen = 100)
  DT <- rateOfChange[Year %in% c(2011, 2051) & percArea > 0,]
  DT[, percArea := percArea*100]
  DT[, rateChange := rateChange*100]
  DT[, totalRateChange := totalRateChange*100]
  write.csv(DT, file.path(getwd(), "outputs", "rateOfChange.csv"))
  
  # 2. Rate of change of buffered disturbances for all sectors within caribou ranges
  # For year 2011 and 2051:
  shpRas <- lapply(shapesNames, function(TYPE){
    # Buffer all shapefiles
    # Convert shapefiles to rasters
    tic(paste0("Time elapsed for converting ", TYPE, " from shp to raster:"))
    shp <- get(TYPE)
    first <- shp[[1]]
    last <- shp[[length(shp)]]
    bFirst <- aggregate(buffer(first, width = 500))
    bLast <- aggregate(buffer(last, width = 500))
    bFirstSF <- st_as_sf(bFirst)
    bLastSF <- st_as_sf(bLast)
    bFirstRas <- fasterize::fasterize(sf = bFirstSF, raster = oil[[1]])
    names(bFirstRas) <- TYPE
    bLastRas <- fasterize::fasterize(sf = bLastSF, raster = oil[[1]])
    names(bLastRas) <- TYPE
    toc()
    return(list(Year2011 = bFirstRas, Year2051 = bLastRas))
  })
  # Organize the list!
  
  Year2011shp <- raster::stack(lapply(shpRas, `[[`, "Year2011"))
  names(Year2011shp) <- paste0(shapesNames, "_2011")
  writeRaster(Year2011shp, filename = file.path(getwd(), "outputs", "shpRas_2011"), 
              format = "GTiff")
  Year2051shp <- raster::stack(lapply(shpRas, `[[`, "Year2051"))
  names(Year2051shp) <- paste0(shapesNames, "_2051")
  writeRaster(Year2051shp, filename = file.path(getwd(), "outputs", "shpRas_2051"), 
              format = "GTiff")
  
  # Buffer all rasters
  ras <- lapply(rasNames, function(TYPE){
    tic(paste0("Time elapsed for buffering ", TYPE, ": "))
    r <- get(TYPE)
    first <- r[[1]]
    last <- r[[nlayers(r)]]
    first <- rast(first)
    last <- rast(last)
    first[] <- first[]
    last[] <- last[]
    # bFirst <- terra::buffer(first, width = 500) # ~~~~~~~~~~~~> SUPER SLOW!!!
    # bLast <- terra::buffer(last, width = 500) # ~~~~~~~~~~~~> SUPER SLOW!!!
    
    # Below works because we are interested in 500m (ECCC guideline) which is,
    # coincidentaly, double the resolution of the current rasters (250m)
    # This should be generalized for other resolutions before used so!
    # To generalize, we need to identify the res/2 and then iterate over 
    # the number of pixels we need --> we reuse the identified adjecent cells and 
    # identify their adjecent ones until the total buffer is reached. Might end up 
    # being too time consuming, though. Should see the 'buffer' mechanism for bigger 
    # buffers. 
    bFirst <- bufferCells(ras = first, 
                          valToLookFor = 1, 
                          newValue = 1)
    bLast <- bufferCells(ras = last, 
                         valToLookFor = 1, 
                         newValue = 1)
    
    # Buffer is taking a huge amount of time. 
    # Finding the adjecent cells in all directions of the cells that are == 1
    # is much quicker. Then convert those to 1
    toc()
    return(list(Year2011 = raster(bFirst), 
                Year2051 = raster(bLast)))
  })
  
  Year2011ras <- stack(lapply(ras, `[[`, "Year2011"))
  names(Year2011ras) <- paste0(rasNames, "_2011")
  writeRaster(Year2011ras, filename = file.path(getwd(), "outputs", "rasRas_2011"), 
              format = "GTiff")
  Year2051ras <- stack(lapply(ras, `[[`, "Year2051"))
  names(Year2051ras) <- paste0(rasNames, "_2051")
  writeRaster(Year2051ras, filename = file.path(getwd(), "outputs", "rasRas_2051"), 
              format = "GTiff")
  
  # Join all rasters
  Year2011 <- stack(Year2011shp, Year2011ras)
  Year2051 <- stack(c(Year2051shp, Year2051ras))
  
  # Convert all pixels that are disturbed to 1, and non-disturbed to 0
  tempRas <- wind[[1]]
  tempRas[] <- tempRas[]
  tempRas[!is.na(tempRas[])] <- 0
  names(tempRas) <- "templateRas"
  
  Year2011_all <- tempRas
  for (i in 1:nlayers(Year2011)){
    tic(paste0("Updating ", names(Year2011)[i], ": "))
    Year2011_all[which(Year2011[[i]][] == 1)] <- 1
    toc()
  }
  
  Year2051_all <- tempRas
  for (i in 1:nlayers(Year2051)){
    tic(paste0("Updating ", names(Year2051)[i], ": "))
    Year2051_all[which(Year2051[[i]][] == 1)] <- 1
    toc()
  }
  
  # Summarize by range polygon in terms of pixels & range plan
  # FOR 2011 ##########################################################
  Year2011 <- rast(Year2011_all)
  # BOO SHAPEFILE
  booDist2011 <- terra::extract(x = Year2011, y = Boo)
  booDist2011_DT <- data.table(booDist2011)
  names(booDist2011_DT) <- c("groupID", "boo_2011")
  booDist2011_DT <- merge(booDist2011_DT, data.table(groupID = 1:NROW(Boo),
                                                     Boo[["REGION"]]))
  names(booDist2011_DT) <- c("groupID", "boo_2011", "REGION")
  cols <- c("notDisturbed", "disturbed")
  booDist2011_DT[, (cols) := .(sum(boo_2011 == 0, na.rm = TRUE), 
                               sum(boo_2011 == 1, na.rm = TRUE)), 
                 by = "REGION"]
  Boo_2011 <- unique(booDist2011_DT[, c("REGION", "notDisturbed", "disturbed")])
  Boo_2011[, Year := 2011]
  Boo_2011[, proportionDisturbed := disturbed/(notDisturbed+disturbed)]
  Boo_2011[, totalAreaHa := expanse(Boo, unit="ha")]
  
  # HERD SHAPEFILE
  herdsDist2011 <- terra::extract(x = Year2011, y = Herds)
  herdsDist2011 <- data.table(herdsDist2011)
  names(herdsDist2011) <- c("groupID", "herd_2011")
  herdsDist2011 <- merge(herdsDist2011, data.table(groupID = 1:NROW(Herds),
                                                   Herds[["HERD"]]))
  names(herdsDist2011) <- c("groupID", "herd_2011", "REGION")
  cols <- c("notDisturbed", "disturbed")
  herdsDist2011[, (cols) := .(sum(herd_2011 == 0, na.rm = TRUE), 
                              sum(herd_2011 == 1, na.rm = TRUE)), 
                by = "REGION"]
  Herd_2011 <- unique(herdsDist2011[, c("REGION", "notDisturbed", "disturbed")])
  Herd_2011[, Year := 2011]
  Herd_2011[, proportionDisturbed := disturbed/(notDisturbed+disturbed)]
  Herd_2011[, totalAreaHa := expanse(Herds, unit="ha")]
  
  # FOR 2051 ##########################################################
  Year2051 <- rast(Year2051_all)
  # BOO SHAPEFILE
  booDist2051 <- terra::extract(x = Year2051, y = Boo)
  booDist2051_DT <- data.table(booDist2051)
  names(booDist2051_DT) <- c("groupID", "boo_2051")
  booDist2051_DT <- merge(booDist2051_DT, data.table(groupID = 1:NROW(Boo),
                                                     Boo[["REGION"]]))
  names(booDist2051_DT) <- c("groupID", "boo_2051", "REGION")
  cols <- c("notDisturbed", "disturbed")
  booDist2051_DT[, (cols) := .(sum(boo_2051 == 0, na.rm = TRUE), 
                               sum(boo_2051 == 1, na.rm = TRUE)), 
                 by = "REGION"]
  Boo_2051 <- unique(booDist2051_DT[, c("REGION", "notDisturbed", "disturbed")])
  Boo_2051[, Year := 2051]
  Boo_2051[, proportionDisturbed := disturbed/(notDisturbed+disturbed)]
  Boo_2051[, totalAreaHa := expanse(Boo, unit="ha")]
  
  # HERD SHAPEFILE
  herdsDist2051 <- terra::extract(x = Year2051, y = Herds)
  herdsDist2051 <- data.table(herdsDist2051)
  names(herdsDist2051) <- c("groupID", "herd_2051")
  herdsDist2051 <- merge(herdsDist2051, data.table(groupID = 1:NROW(Herds),
                                                   Herds[["HERD"]]))
  names(herdsDist2051) <- c("groupID", "herd_2051", "REGION")
  cols <- c("notDisturbed", "disturbed")
  herdsDist2051[, (cols) := .(sum(herd_2051 == 0, na.rm = TRUE), 
                              sum(herd_2051 == 1, na.rm = TRUE)), 
                by = "REGION"]
  Herd_2051 <- unique(herdsDist2051[, c("REGION", "notDisturbed", "disturbed")])
  Herd_2051[, Year := 2051]
  Herd_2051[, proportionDisturbed := disturbed/(notDisturbed+disturbed)]
  Herd_2051[, totalAreaHa := expanse(Herds, unit="ha")]
  
  # JOIN ALL DATA!
  DT <- rbind(Boo_2011, Boo_2051,
              Herd_2011, Herd_2051)
  write.csv(x = DT, file = file.path(getwd(), "outputs", "disturbedTableByRegionRaw.csv"))
  DT[, c("notDisturbed", "disturbed") := NULL]
  DTc <- dcast(DT, REGION + totalAreaHa ~ Year, value.var = "proportionDisturbed")
  DTc[, totalDisturbanceChange := ((`2051`-`2011`)/`2011`)]  # Total increase over 2011 in 40 years
  DTc[, disturbanceChangePerYear := ((`2051`-`2011`)/`2011`)/(2051-2011)] # Avergare yearly increase between 2011 and 2051 
  DTc[, disturbanceChangePerYearPerc := disturbanceChangePerYear*100]
  names(DTc)[names(DTc) == 2011] <- "disturbedPercent2011"
  names(DTc)[names(DTc) == 2051] <- "disturbedPercent2051"
  write.csv(x = DTc, file = file.path(getwd(), "outputs", "disturbedTableByRegion.csv"))
  
  # 3. Time series of total disturbance buffered
  # Get maps from previous calcs
  Require("lattice")
  Require("rasterVis")
  Require("viridis")
  Require("maptools")
  
  both <- Year2011+Year2051
  both <- crop(both, Boo)
  both <- mask(both, Boo)
  both[both == 2] <- 2051
  both[both == 1] <- 2011
  both[is.na(both)] <- 0
  bothRat <- as.factor(both)
  bothRat <- raster(bothRat)
  att <- "ID"
  
  Pal <- c("white", "cornflowerblue", "darkslateblue") 
  
  plot(bothRat, col = Pal)
  bothRat <- rast(bothRat)
  bothRat <- mask(bothRat, Boo)
  
  png(filename = file.path(getwd(), "outputs", "totalDisturbance.png"),
      width = 21, height = 29,
      units = "cm", res = 300)
  par(bg = NA)
  raster::plot(bothRat, 
               main = "Total disturbance buffered in 2011 and 2051",
               col = c("white", "orange", "darkslateblue"),
               axes = FALSE, 
               box = FALSE)
  plot(Boo, add = TRUE, lwd = 1, 
       border = "black")
  dev.off()
  
  # 4. Time series of non-buffered disturbance types (Edéhzhíe National Wildlife Area region)
  
  nonBufferedDist <- lapply(c(rasNames, shapesNames), function(TYPE){
    tic(paste0("Time elapsed for ", TYPE, ":"))
    dist <- get(TYPE)
    if (class(dist) == "RasterStack"){
      dist[] <- dist[]
      allYears <- lapply(rev(1:nlayers(dist)), function(INDEX){
        d <- dist[[INDEX]]
        d <- rast(d)
        # Crop and mask to Ede
        d <- crop(d, Ede)
        d <- mask(d, Ede)
        # Convert to SpatVect
        distVec <- as.polygons(d,
                               trunc = TRUE, 
                               dissolve = TRUE, 
                               values = FALSE,
                               na.rm = TRUE, 
                               na.all = TRUE, 
                               extent = FALSE)
        # Add data info
        Y <- 2001+(INDEX*10)
        distVec[["TYPE"]] <- TYPE
        distVec[["YEAR"]] <- Y 
        distVec[["Type_Year"]] <- paste0(TYPE, "_", Y)
        distVec[["sizeHa"]] <- expanse(distVec, unit = "ha")
        distVec[["Class"]] <- TYPE  
        distVec <- subset(x = distVec, subset = distVec$sizeHa == min(distVec$sizeHa), )
        return(distVec)
      })
    } else {
      allYears <- lapply(rev(1:length(dist)), function(INDEX){
        d <- dist[[INDEX]]
        if (class(d) != "SpatVector")
          d <- vect(d)
        # Crop and mask to Ede
        d <- crop(d, Ede)
        distVec <- mask(d, Ede)
        # Add data info
        Y <- 2001+(INDEX*10)
        distVec[["TYPE"]] <- TYPE
        distVec[["YEAR"]] <- Y 
        distVec[["Type_Year"]] <- paste0(TYPE, "_", Y)
        distVec[["sizeHa"]] <- expanse(distVec, unit = "ha")
        # If not polygon, convert
        if (geomtype(distVec) != "polygons")
          distVec <- buffer(distVec, width = 7.5) ## The 7.5 is the half of the 15m, resolution of original data
        return(distVec)
      })
    }
    allYearsShp <- do.call(rbind, allYears)
    toc()
    # PLOTTING
    whichPlot <- "YEAR"
    par(bg = NA)
    plot(x = allYearsShp,
         y = whichPlot,
         type = "classes",
         col = rainbow(length(unlist(allYearsShp[[whichPlot]]))),
         border = rainbow(length(unlist(allYearsShp[[whichPlot]]))),
         box = FALSE,
         axes = FALSE,
         main = TYPE
    )
    plot(Ede, add = TRUE, border = "black")
    dev.copy(png, file.path(getwd(), "outputs", paste0("disturbance_",TYPE,".png")))
    dev.off()
    return(allYearsShp)
  })
  nonBufferedDist <- do.call(rbind, nonBufferedDist)
  
  # Now subset year == 2011 and plot all types, and year == 20151 and plot all types
  allDist2011 <- subset(nonBufferedDist, nonBufferedDist[["YEAR"]] == 2011)
  whichPlot <- "TYPE"
  par(bg = NA)
  COLS <- c("red", "orange", "goldenrod1", "forestgreen", "blue", "navy", "purple")
  for (N in 1:length(unique(allDist2011$TYPE)[!unique(allDist2011$TYPE) %in% c("wind", 
                                                                               "powerlines")])){
    i <- unique(allDist2011$TYPE)[!unique(allDist2011$TYPE) %in% c("wind", 
                                                                   "powerlines")][N]
    print(paste0("Making the figure for ", i))
    plot.window(xlim = c(xmin(Ede), xmax(Ede)), 
                ylim = c(ymin(Ede), ymax(Ede)))
    plot(x = allDist2011[allDist2011$TYPE == i],
         y = whichPlot,
         type = "classes",
         col = COLS[N],
         border = COLS[N],
         box = FALSE,
         axes = FALSE,
         main = NA
    )
    plot(Ede, add = TRUE, border = "black")
    dev.copy(png,file.path(getwd(), "outputs", paste0("disturbance_",i,"2011.png")))
    dev.off()
  }
  writeOGR(obj = as(st_as_sf(allDist2011), "Spatial"), dsn = file.path(getwd(), "outputs"), 
           layer = "disturbances_all_2011", 
           driver = "ESRI Shapefile")
  
  allDist2051 <- subset(nonBufferedDist, nonBufferedDist[["YEAR"]] == 2051)
  for (N in 1:length(unique(allDist2051$TYPE)[!unique(allDist2051$TYPE) %in% c("wind", 
                                                                               "powerlines")])){
    i <- unique(allDist2051$TYPE)[!unique(allDist2051$TYPE) %in% c("wind", 
                                                                   "powerlines")][N]
    print(paste0("Making the figure for ", i))
    plot.window(xlim = c(xmin(Ede), xmax(Ede)), 
                ylim = c(ymin(Ede), ymax(Ede)))
    plot(x = allDist2051[allDist2051$TYPE == i],
         y = whichPlot,
         type = "classes",
         col = COLS[N],
         border = COLS[N],
         box = FALSE,
         axes = FALSE,
         main = NA
    )
    plot(Ede, add = TRUE, border = "black")
    dev.copy(png,file.path(getwd(), "outputs", paste0("disturbance_",i,"2051.png")))
    dev.off()
  }
  writeOGR(obj = as(st_as_sf(allDist2051), "Spatial"), dsn = file.path(getwd(), "outputs"), 
           layer = "disturbances_all_2051", 
           driver = "ESRI Shapefile")
  
  
  # 5. Upload all data
  # [ ] Raw maps
  Raw <- list.files(path = file.path(getwd(), "outputs"), 
                    pattern = "disturbances_", full.names = TRUE)
  # [ ] Results 1
  R1 <- file.path(getwd(), "outputs", "rateOfChange.csv") # STILL TO RUN!
  # [ ] Results 2
  R2 <- c(file.path(getwd(), "outputs", "disturbedTableByRegion.csv"), 
          file.path(getwd(), "outputs", "disturbedTableByRegionRaw.csv"))
  
  # [ ] Results 3
  R3 <- c(file.path(getwd(), "outputs", "totalDisturbance.png"),
          file.path(getwd(), "outputs", "rasRas_2051.tif"),
          file.path(getwd(), "outputs", "rasRas_2011.tif"),
          file.path(getwd(), "outputs", "shpRas_2051.tif"),
          file.path(getwd(), "outputs", "shpRas_2011.tif"))
  
  lapply(R1, drive_upload, path = as_id("108QahyJg35V0C7nq0CR6VTM1rA6Lvvcm"))
  lapply(R2, drive_upload, path = as_id("108QahyJg35V0C7nq0CR6VTM1rA6Lvvcm"))
  lapply(R3, drive_upload, path = as_id("108QahyJg35V0C7nq0CR6VTM1rA6Lvvcm"))
  lapply(Raw, drive_upload, path = as_id("108QahyJg35V0C7nq0CR6VTM1rA6Lvvcm"))
  
  
  lastFile <- file.path(getwd(), "outputs", "disturbances_forestry_cutblocks_2041.tif.aux.xml")
  
  Raw2 <- Raw[which(Raw == lastFile):length(Raw)]
  Require("googledrive")
  lapply(Raw2, drive_upload, path = as_id("108QahyJg35V0C7nq0CR6VTM1rA6Lvvcm"))
  
}
