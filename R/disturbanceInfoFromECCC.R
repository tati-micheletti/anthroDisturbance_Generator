disturbanceInfoFromECCC <- function(studyArea, 
                                    RTM, 
                                    classesAvailable,
                                    totalstudyAreaVAreaSqKm,
                                    disturbanceList,
                                    diffYears = "2010_2015",
                                    destinationPath,
                                    bufferedDisturbances = TRUE,
                                    maskOutLinesFromPolys = TRUE,
                                    aggregateSameDisturbances = FALSE,
                                    archiveNEW = "ECCC_2015_anthro_dist_corrected_to_NT1_2016_final.zip",
                                    targetFileNEW = c("BEADlines2015_NWT_corrected_to_NT1_2016.shp", 
                                                      "BEADpolys2015_NWT_corrected_to_NT1_2016.shp"),
                                    urlNEW = paste0("https://drive.google.com/file/d/1sxAa0wwwt7iwiHD7zB0DDnjfqyIQjKI2"), # James Hodsons' corrections,
                                    archiveOLD = paste0("Boreal-ecosystem-anthropogenic-disturbance-",
                                                        "vector-data-2008-2010.zip"),
                                    targetFileOLD = c("EC_borealdisturbance_linear_2008_2010_FINAL_ALBERS.shp", 
                                                      "EC_borealdisturbance_polygonal_2008_2010_FINAL_ALBERS.shp"),
                                    urlOLD = paste0("https://www.ec.gc.ca/data_donnees/STB-DGST/003/Boreal-ecosystem",
                                                    "-anthropogenic-disturbance-vector-data-2008-OLD.zip")){
  
  # If the table doesn't have disturbance rate, we can calculate it based on data!
  
  # NOTE [UPDATE 22-12-23]: If using a study area different (i.e., buffered) than just the caribou NT1 range, for example,
  # the 2015 data seems to be cropped already, which results in much smaller total disturbed area than 2015.
  # For now, we will generate the scenarios without getting information from data (i.e., using the % calculated 
  # over the known size of area and total disturbed area).
  
  # NOTE: According to ECCC report -- 2019 Species at Risk Act Conservation Agreement for the 
  # Conservation of the Boreal Caribou new human disturbances total to 0.2% of the total area 
  # (NT1) a year 
  
  # Calculate the difference in time between both layers 
  parts <- strsplit(diffYears, "_")[[1]]
  num1 <- suppressWarnings(as.numeric(parts[1]))
  num2 <- suppressWarnings(as.numeric(parts[2]))
  yearDistance <- num2 - num1
  
  # Original ECCC file:
  # urlNEW <- paste0("https://data-donnees.ec.gc.ca/data/species/developplans/NEW-",
  #                   "anthropogenic-disturbance-footprint-within-boreal-caribou-ranges",
  #                   "-across-canada-as-interpreted-from-NEW-landsat-satellite-imagery/",
  #                   "Anthro_Disturb_Perturb_30m_NEW.zip")
  # archive = "NorthwestTerritories_30m_Disturb_Perturb_Line.zip",
  # targetFile = "NorthwestTerritories_30m_Disturb_Perturb_Line.shp",
  AD_NEW_Lines <- prepInputs(url = urlNEW,
                              archive = archiveNEW,
                              alsoExtract = "similar",
                              studyArea = studyArea,
                              rasterToMatch = RTM,
                              fun = "terra::vect",
                              targetFile = targetFileNEW[1],
                              destinationPath = destinationPath)
  if (!is(AD_NEW_Lines, "SpatVector"))
    AD_NEW_Lines <- terra::vect(AD_NEW_Lines)
  
  AD_NEW_Polys <- prepInputs(url = urlNEW,
                              alsoExtract = "similar",
                              archive = archiveNEW,
                              studyArea = studyArea,
                              rasterToMatch = RTM,
                              fun = "terra::vect",
                              targetFile = targetFileNEW[2],
                              destinationPath = destinationPath)
  if (!is(AD_NEW_Polys, "SpatVector"))
    AD_NEW_Polys <- terra::vect(AD_NEW_Polys)
  
  AD_OLD_Lines <- prepInputs(url = urlOLD,
                              archive = archiveOLD,
                              alsoExtract = "similar",
                              studyArea = studyArea,
                              rasterToMatch = RTM,
                              fun = "terra::vect",
                              targetFile = targetFileOLD[1],
                              destinationPath = destinationPath)
  if (!is(AD_OLD_Lines, "SpatVector"))
    AD_OLD_Lines <- terra::vect(AD_OLD_Lines)
  
  AD_OLD_Polys <- prepInputs(url = urlOLD,
                              archive = archiveOLD,
                              alsoExtract = "similar",
                              studyArea = studyArea,
                              rasterToMatch = RTM,
                              fun = "terra::vect",
                              targetFile = targetFileOLD[2],
                              destinationPath = destinationPath)
  if (!is(AD_OLD_Polys, "SpatVector"))
    AD_OLD_Polys <- terra::vect(AD_OLD_Polys)
  
  bufferSize <- if (bufferedDisturbances) 500 else 30

  # First we need to buffer all lines to 30 (bufferedDisturbances = FALSE) as they are not measured 
  # otherwise, or all disturbances (bufferedDisturbances = TRUE) to 500m. The 30m come from the 
  # resolution of the original data.
  
  AD_NEW_Lines <- terra::buffer(x = AD_NEW_Lines, width = bufferSize)
  if (bufferedDisturbances)
    AD_NEW_Polys <- terra::buffer(x = AD_NEW_Polys, width = bufferSize)
  
  AD_OLD_Lines <- terra::buffer(x = AD_OLD_Lines, width = bufferSize)
  if (bufferedDisturbances)
    AD_OLD_Polys <- terra::buffer(x = AD_OLD_Polys, width = bufferSize)
  
  if (maskOutLinesFromPolys){
    # Exclude buffered polygons from lines to not double count them
    AD_NEW_Lines <- terra::mask(AD_NEW_Lines, mask = AD_NEW_Polys, inverse = TRUE)
    AD_OLD_Lines <- terra::mask(AD_OLD_Lines, mask = AD_OLD_Polys, inverse = TRUE)
  }
  
  # Bind both data so I can extract the area 
  AD_NEW_all <- rbind(AD_NEW_Lines, AD_NEW_Polys)
  AD_OLD_all <- rbind(AD_OLD_Lines, AD_OLD_Polys)
  
  # Class conversion needs to happen before aggregation
  # Cleanup the data. Some classes are not in both datasets.
  # NEW
  # TODO Implement here a cleanup method: whatever is not in one list, needs to be deleted from the other!
  # The whole NT1 needs only the "NotDisturbance" cleaned up
  # toChange <- which(AD_NEW_all$Class == "Well site")
  # AD_NEW_all$Class[toChange] <- "Oil/Gas"
  # toChange <- which(AD_NEW_all$Class == "Airstrip")
  # AD_NEW_all$Class[toChange] <- "Road"
  AD_NEW_all <- subset(AD_NEW_all, AD_NEW_all$Class != "NotDisturbance", select = "Class")
  # OLD
  # toChange <- which(AD_OLD_all$Class == "Well site")
  # AD_OLD_all$Class[toChange] <- "Oil/Gas"
  # toChange <- which(AD_OLD_all$Class == "Reservoir")
  # AD_OLD_all$Class[toChange] <- "Oil/Gas"  
  # toChange <- which(AD_OLD_all$Class == "Powerline")
  # AD_OLD_all$Class[toChange] <- "Pipeline"

  # Check if all disturbances that are currently in the area (which are available
  # in the layer disturbanceList) are being aggregated to the NEW layer
  # These might have been missed by ECCC's product due to resolution, for example
  # We assume they are newer rather than older.
  allClassesAvailable <- unique(AD_NEW_all$Class)
  nonPotLay <- extractNonPotentialLayers(disturbanceList)
  laysToADD <- lapply(1:nrow(nonPotLay), function(INDEX){
    layIndex <- disturbanceList[[nonPotLay[INDEX, Sector]]][[nonPotLay[INDEX, dataClass]]]
    if (geomtype(layIndex) != geomtype(AD_NEW_all)){
      layIndex <- terra::buffer(x = layIndex, width = bufferSize)
    }
    return(layIndex)
    })
  laysToADD <- do.call(rbind, laysToADD)
  AD_NEW_all <- rbind(AD_NEW_all, laysToADD)
  
  if (aggregateSameDisturbances){
    AD_OLD_all <- terra::aggregate(AD_OLD_all, by = "Class", dissolve = TRUE)
    AD_NEW_all <- terra::aggregate(AD_NEW_all, by = "Class", dissolve = TRUE)
  }
  
  AD_NEW_all$Area_sqKm <- terra::expanse(AD_NEW_all, transform = FALSE, unit = "km")
  AD_OLD_all$Area_sqKm <- terra::expanse(AD_OLD_all, transform = FALSE, unit = "km")

    # 3.2. Bring all to a data.table and summarize
  AD_NEW_DT <- as.data.table(as.data.frame(AD_NEW_all[, c("Class", "Area_sqKm")]))
  AD_OLD_DT <- as.data.table(as.data.frame(AD_OLD_all[, c("Class", "Area_sqKm")]))
  
  # Test if after data cleanup we have the same classes in both datasets 
  if (!all(unique(AD_NEW_DT$Class) %in% unique(AD_OLD_DT$Class))){
    warning(paste0("Not all classes of OLD are in the NEW data. \nThis might",
                   " be due to the lack of the disturbances at that point."), immediate. = TRUE)
    missingClassOLD <- setdiff(AD_NEW_DT$Class, AD_OLD_DT$Class)
    missingClassNEW <- setdiff(AD_OLD_DT$Class, AD_NEW_DT$Class)
    if (length(missingClassOLD) > 0){
      message(paste0("Adding the following disturbances as 0 to OLD: "))
      paste(missingClassOLD, sep = ", ")
      AD_OLD_DT <- rbind(AD_OLD_DT, data.table(Class = missingClassOLD,
                                                     Area_sqKm = 0))
    }
    if (length(missingClassNEW) > 0){
      message(paste0("Adding the following disturbances as 0 to NEW: "))
      paste(missingClassNEW, sep = ", ")
      AD_NEW_DT <- rbind(AD_NEW_DT, data.table(Class = missingClassNEW,
                                                 Area_sqKm = 0))
    }
  }

  AD_NEW_DT_summ <- AD_NEW_DT[, totalArea_sqKm := sum(Area_sqKm), by = "Class"]
  AD_NEW_DT_summ <- unique(AD_NEW_DT_summ[, Area_sqKm := NULL]) # For when aggregateSameDisturbances = FALSE
  AD_OLD_DT_summ <- AD_OLD_DT[, totalArea_sqKm := sum(Area_sqKm), by = "Class"]
  AD_OLD_DT_summ <- unique(AD_OLD_DT_summ[, Area_sqKm := NULL]) # For when aggregateSameDisturbances = FALSE
  
  # 3.3. Now I can see how much the disturbances changed from OLD to NEW in sq Km
  AD_NEW_DT_summ[, Year := "yearNEW"]
  AD_OLD_DT_summ[, Year := "yearOLD"]
  
  AD_change <- rbind(AD_NEW_DT_summ, AD_OLD_DT_summ)
  
  AD_changed <- dcast(data = AD_change, formula = Class ~ Year, value.var = "totalArea_sqKm")
  AD_changed[is.na(yearOLD), yearOLD := 0]
  
  AD_changed[, disturbProportionInAreaOLD := yearOLD/totalstudyAreaVAreaSqKm]
  AD_changed[, disturbProportionInAreaNEW := yearNEW/totalstudyAreaVAreaSqKm]

  AD_changed[, totalProportionAreaDisturbedOLD := sum(disturbProportionInAreaOLD)]
  AD_changed[, totalProportionAreaDisturbedNEW := sum(disturbProportionInAreaNEW)]

  AD_changed[, totalPercentAreaDisturbedOLD := totalProportionAreaDisturbedOLD*100]
  AD_changed[, totalPercentAreaDisturbedNEW := totalProportionAreaDisturbedNEW*100]

  AD_changed[, relativeChangePerClassPerYear := ((yearNEW-yearOLD)/yearOLD)/yearDistance]
  AD_changed[, relativeChangePerClassPerYearPerc := relativeChangePerClassPerYear*100]

  # [ NOTE: Using percent is not a good practice because the really high percentages have a really small total area.
  # And the largest change (seismic lines in smaller study area, reduction has a HUGE area but comparatively really small %)

  message(paste0("Percentage of total human disturbance ", 
                 ifelse(bufferedDisturbances, " buffered by 500m ", " unbuffered "),
                 ifelse(aggregateSameDisturbances, " aggregated by class ", " unaggregated "),
                 ifelse(maskOutLinesFromPolys, " with masked out lines from polygons ", 
                        " without masking lines out of polygons "),
                 " in the study area based on OLD (", num1,") and NEW (", num2,") data are: ",
                 paste0(round(unique(AD_changed[["totalPercentAreaDisturbedOLD"]]), 3),
                        "% and ",
                        round(unique(AD_changed[["totalPercentAreaDisturbedNEW"]]), 3),
                        "%, respectively.")))

    # Here we used all disturbances to calculate the proportion that belongs to roads and other lines and polys
    # This helps building the connection between totalDisturbanceRate and each disturbance class
    toFill <- AD_changed[Class %in% classesAvailable[["classToSearch"]]]
    toUse <- merge(toFill, classesAvailable[, c("classToSearch", "dataClass")], all.x = TRUE, 
                    by.x = "Class", by.y = "classToSearch")
    toUse <- toUse[, c("dataClass", "disturbProportionInAreaOLD", "disturbProportionInAreaNEW")]
    toUse[, proportionAreaSqKmChangedPerYear := (disturbProportionInAreaNEW-disturbProportionInAreaOLD)/yearDistance]
    toUse <- toUse[, c("dataClass", "proportionAreaSqKmChangedPerYear")]
    
    proportionTable <- dcast(toUse, dataClass ~ ., fun.agg = sum, 
                     value.var = "proportionAreaSqKmChangedPerYear")
    names(proportionTable) <- c("dataClass", "proportionAreaSqKmChangedPerYear")
    
    # If we have anything reducing from one year to the next, at this time we will remove from here
    if (any(proportionTable[["proportionAreaSqKmChangedPerYear"]] < 0)){
      whichDistReducing <- proportionTable[["dataClass"]][proportionTable[["proportionAreaSqKmChangedPerYear"]] < 0]
      # Percentage of total disturbance belonging to each sector
      warning(paste0("The disturbance(s) '", paste0(whichDistReducing, collapse = ", "), "' are reducing ",
                     "between the years OLD and NEW in the study area according to ECCC data. These ",
                     "will be excluded from the calculation of totalDisturbance"), immediate. = TRUE)
      proportionTable <- proportionTable[proportionAreaSqKmChangedPerYear >= 0, ]
    }
    # ttDistubanceRate: total disturbance proportion across the studyArea per year
    ttDistubanceRate <- sum(proportionTable[["proportionAreaSqKmChangedPerYear"]])
    proportionTable[, proportionOfTotalDisturbance := proportionAreaSqKmChangedPerYear/ttDistubanceRate]
    write.csv(x = AD_changed, file.path(destinationPath, 
                                        paste0("anthropogenicDisturbance_ECCC_", diffYears,"_", digest(studyArea),".csv")))
    write.csv(x = proportionTable, file.path(destinationPath, 
                                             paste0("proportionTable_ECCC_", diffYears,"_", digest(studyArea),".csv")))
    return(list(AD_changed = AD_changed,
                proportionTable = proportionTable))
}

