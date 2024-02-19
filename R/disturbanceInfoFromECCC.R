disturbanceInfoFromECCC <- function(studyArea, 
                                    RTM, 
                                    classesAvailable,
                                    totalstudyAreaVAreaSqKm,
                                    destinationPath,
                                    bufferedDisturbances = TRUE,
                                    maskOutLinesFromPolys = TRUE,
                                    aggregateSameDisturbances = FALSE){
  # If the table doesn't have disturbance rate, we can calculate it based on data!
  
  # NOTE [UPDATE 22-12-23]: If using a study area different (i.e., buffered) than just the caribou NT1 range, for example,
  # the 2015 data seems to be cropped already, which results in much smaller total disturbed area than 2015.
  # For now, we will generate the scenarios without getting information from data (i.e., using the % calculated 
  # over the known size of area and total disturbed area).
  
  # NOTE: According to ECCC report -- 2019 Species at Risk Act Conservation Agreement for the 
  # Conservation of the Boreal Caribou new human disturbances total to 0.2% of the total area 
  # (NT1) a year 
  
  # Original ECCC file:
  # url2015 <- paste0("https://data-donnees.ec.gc.ca/data/species/developplans/2015-",
  #                   "anthropogenic-disturbance-footprint-within-boreal-caribou-ranges",
  #                   "-across-canada-as-interpreted-from-2015-landsat-satellite-imagery/",
  #                   "Anthro_Disturb_Perturb_30m_2015.zip")
  # archive = "NorthwestTerritories_30m_Disturb_Perturb_Line.zip",
  # targetFile = "NorthwestTerritories_30m_Disturb_Perturb_Line.shp",
  url2015 <- paste0("https://drive.google.com/file/d/1sxAa0wwwt7iwiHD7zB0DDnjfqyIQjKI2") # James Hodsons' corrections
  AD_2015_Lines <- prepInputs(url = url2015,
                              archive = "ECCC_2015_anthro_dist_corrected_to_NT1_2016_final.zip",
                              alsoExtract = "similar",
                              studyArea = studyArea,
                              rasterToMatch = RTM,
                              fun = "terra::vect",
                              targetFile = "BEADlines2015_NWT_corrected_to_NT1_2016.shp",
                              destinationPath = destinationPath)
  if (!is(AD_2015_Lines, "SpatVector"))
    AD_2015_Lines <- terra::vect(AD_2015_Lines)
  
  AD_2015_Polys <- prepInputs(url = url2015,
                              alsoExtract = "similar",
                              archive = "ECCC_2015_anthro_dist_corrected_to_NT1_2016_final.zip",
                              studyArea = studyArea,
                              rasterToMatch = RTM,
                              fun = "terra::vect",
                              targetFile = "BEADpolys2015_NWT_corrected_to_NT1_2016.shp",
                              destinationPath = destinationPath)
  if (!is(AD_2015_Polys, "SpatVector"))
    AD_2015_Polys <- terra::vect(AD_2015_Polys)
  
  url_2010 <- paste0("https://www.ec.gc.ca/data_donnees/STB-DGST/003/Boreal-ecosystem",
                     "-anthropogenic-disturbance-vector-data-2008-2010.zip")
  AD_2010_Lines <- prepInputs(url = url_2010,
                              archive = paste0("Boreal-ecosystem-anthropogenic-disturbance-",
                                               "vector-data-2008-2010.zip"),
                              alsoExtract = "similar",
                              studyArea = studyArea,
                              rasterToMatch = RTM,
                              fun = "terra::vect",
                              targetFile = "EC_borealdisturbance_linear_2008_2010_FINAL_ALBERS.shp",
                              destinationPath = destinationPath)
  if (!is(AD_2010_Lines, "SpatVector"))
    AD_2010_Lines <- terra::vect(AD_2010_Lines)
  
  AD_2010_Polys <- prepInputs(url = url_2010,
                              archive = "Boreal-ecosystem-anthropogenic-disturbance-vector-data-2008-2010.zip",
                              alsoExtract = "similar",
                              studyArea = studyArea,
                              rasterToMatch = RTM,
                              fun = "terra::vect",
                              targetFile = "EC_borealdisturbance_polygonal_2008_2010_FINAL_ALBERS.shp",
                              destinationPath = destinationPath)
  if (!is(AD_2010_Polys, "SpatVector"))
    AD_2010_Polys <- terra::vect(AD_2010_Polys)
  
  bufferSize <- if (bufferedDisturbances) 500 else 30

  # First we need to buffer all lines to 30 (bufferedDisturbances = FALSE) as they are not measured 
  # otherwise, or all disturbances (bufferedDisturbances = TRUE) to 500m. The 30m come from the 
  # resolution of the original data.
  
  AD_2015_Lines <- terra::buffer(x = AD_2015_Lines, width = bufferSize)
  if (bufferedDisturbances)
    AD_2015_Polys <- terra::buffer(x = AD_2015_Polys, width = bufferSize)
  
  AD_2010_Lines <- terra::buffer(x = AD_2010_Lines, width = bufferSize)
  if (bufferedDisturbances)
    AD_2010_Polys <- terra::buffer(x = AD_2010_Polys, width = bufferSize)
  
  if (maskOutLinesFromPolys){
    # Exclude buffered polygons from lines to not double count them
    AD_2015_Lines <- terra::mask(AD_2015_Lines, mask = AD_2015_Polys, inverse = TRUE)
    AD_2010_Lines <- terra::mask(AD_2010_Lines, mask = AD_2010_Polys, inverse = TRUE)
  }
  
  # Bind both data so I can extract the area 
  AD_2015_all <- rbind(AD_2015_Lines, AD_2015_Polys)
  AD_2010_all <- rbind(AD_2010_Lines, AD_2010_Polys)
  
  # Class conversion needs to happen before aggregation
  # Cleanup the data. Some classes are not in both datasets.
  # 2015
  # TODO Implement here a cleanup method: whatever is not in one list, needs to be deleted from the other!
  # The whole NT1 needs only the "NotDisturbance" cleaned up
  # toChange <- which(AD_2015_all$Class == "Well site")
  # AD_2015_all$Class[toChange] <- "Oil/Gas"
  # toChange <- which(AD_2015_all$Class == "Airstrip")
  # AD_2015_all$Class[toChange] <- "Road"
  AD_2015_all <- subset(AD_2015_all, AD_2015_all$Class != "NotDisturbance", select = "Class")
  # 2010
  # toChange <- which(AD_2010_all$Class == "Well site")
  # AD_2010_all$Class[toChange] <- "Oil/Gas"
  # toChange <- which(AD_2010_all$Class == "Reservoir")
  # AD_2010_all$Class[toChange] <- "Oil/Gas"  
  # toChange <- which(AD_2010_all$Class == "Powerline")
  # AD_2010_all$Class[toChange] <- "Pipeline"

  if (aggregateSameDisturbances){
    AD_2010_all <- terra::aggregate(AD_2010_all, by = "Class", dissolve = TRUE)
    AD_2015_all <- terra::aggregate(AD_2015_all, by = "Class", dissolve = TRUE)
  } 
  
    AD_2015_all$Area_sqKm <- terra::expanse(AD_2015_all, transform = FALSE, unit = "km")
    AD_2010_all$Area_sqKm <- terra::expanse(AD_2010_all, transform = FALSE, unit = "km")
  
  # 3.2. Bring all to a data.table and summarize
  AD_2015_DT <- as.data.table(as.data.frame(AD_2015_all[, c("Class", "Area_sqKm")]))
  AD_2010_DT <- as.data.table(as.data.frame(AD_2010_all[, c("Class", "Area_sqKm")]))
  
  # Test if after data cleanup we have the same classes in both datasets 
  if (all(sort(unique(AD_2015_DT$Class)) != sort(unique(AD_2010_DT$Class))))
    stop("Classes of 2010 and 2015 data do not match. Please debug.")

  AD_2015_DT_summ <- AD_2015_DT[, totalArea_sqKm := sum(Area_sqKm), by = "Class"]
  AD_2015_DT_summ <- unique(AD_2015_DT_summ[, Area_sqKm := NULL]) # For when aggregateSameDisturbances = FALSE
  AD_2010_DT_summ <- AD_2010_DT[, totalArea_sqKm := sum(Area_sqKm), by = "Class"]
  AD_2010_DT_summ <- unique(AD_2010_DT_summ[, Area_sqKm := NULL]) # For when aggregateSameDisturbances = FALSE
  
  # 3.3. Now I can see how much the disturbances changed from 2010 to 2015 in sq Km
  AD_2015_DT_summ[, Year := "year2015"]
  AD_2010_DT_summ[, Year := "year2010"]
  
  AD_change <- rbind(AD_2015_DT_summ, AD_2010_DT_summ)
  
  AD_changed <- dcast(data = AD_change, formula = Class ~ Year, value.var = "totalArea_sqKm")

  AD_changed[, disturbProportionInArea2010 := year2010/totalstudyAreaVAreaSqKm]
  AD_changed[, disturbProportionInArea2015 := year2015/totalstudyAreaVAreaSqKm]

  AD_changed[, totalProportionAreaDisturbed2010 := sum(disturbProportionInArea2010)]
  AD_changed[, totalProportionAreaDisturbed2015 := sum(disturbProportionInArea2015)]

  AD_changed[, totalPercentAreaDisturbed2010 := totalProportionAreaDisturbed2010*100]
  AD_changed[, totalPercentAreaDisturbed2015 := totalProportionAreaDisturbed2015*100]

  AD_changed[, relativeChangePerClassPerYear := ((year2015-year2010)/year2010)/5]
  AD_changed[, relativeChangePerClassPerYearPerc := relativeChangePerClassPerYear*100]

  # [ NOTE: Using percent is not a good practice because the really high percentages have a really small total area.
  # And the largest change (seismic lines in smaller study area, reduction has a HUGE area but comparatively really small %)

  message(paste0("Percentage of total human disturbance ", 
                 ifelse(bufferedDisturbances, " buffered by 500m ", " unbuffered "),
                 ifelse(aggregateSameDisturbances, " aggregated by class ", " unaggregated "),
                 ifelse(maskOutLinesFromPolys, " with masked out lines from polygons ", 
                        " without masking lines out of polygons "),
                 " in the study area based on 2010 and 2015 data are: ",
                 paste0(round(unique(AD_changed[["totalPercentAreaDisturbed2010"]]), 3),
                        "% and ",
                        round(unique(AD_changed[["totalPercentAreaDisturbed2015"]]), 3),
                        "%, respectively.")))

    # Here we used all disturbances to calculate the proportion that belongs to roads and other lines and polys
    # This helps building the connection between totalDisturbanceRate and each disturbance class
    toFill <- AD_changed[Class %in% classesAvailable[["classToSearch"]]]
    toUse <- merge(toFill, classesAvailable[, c("classToSearch", "dataClass")], all.x = TRUE, 
                    by.x = "Class", by.y = "classToSearch")
    toUse <- toUse[, c("dataClass", "disturbProportionInArea2010", "disturbProportionInArea2015")]
    toUse[, proportionAreaSqKmChangedPerYear := (disturbProportionInArea2015-disturbProportionInArea2010)/5]
    toUse <- toUse[, c("dataClass", "proportionAreaSqKmChangedPerYear")]
    
    proportionTable <- dcast(toUse, dataClass ~ ., fun.agg = sum, 
                     value.var = "proportionAreaSqKmChangedPerYear")
    names(proportionTable) <- c("dataClass", "proportionAreaSqKmChangedPerYear")
    
    # If we have anything reducing from one year to the next, at this time we will remove from here
    if (any(proportionTable[["proportionAreaSqKmChangedPerYear"]] < 0)){
      whichDistReducing <- proportionTable[["dataClass"]][proportionTable[["proportionAreaSqKmChangedPerYear"]] < 0]
      # Percentage of total disturbance belonging to each sector
      warning(paste0("The disturbance(s) '", paste0(whichDistReducing, collapse = ", "), "' are reducing ",
                     "between the years 2010 and 2015 in the study area according to ECCC data. These ",
                     "will be excluded from the calculation of totalDisturbance"), immediate. = TRUE)
      proportionTable <- proportionTable[proportionAreaSqKmChangedPerYear >= 0, ]
    }
    # ttDistubanceRate: total disturbance proportion across the studyArea per year
    ttDistubanceRate <- sum(proportionTable[["proportionAreaSqKmChangedPerYear"]])
    proportionTable[, proportionOfTotalDisturbance := proportionAreaSqKmChangedPerYear/ttDistubanceRate]
    write.csv(x = AD_changed, file.path(destinationPath, 
                                        "anthropogenicDisturbance_ECCC_2010_2015.csv"))
    write.csv(x = proportionTable, file.path(destinationPath, 
                                             "proportionTable_ECCC_2010_2015.csv"))
    return(list(AD_changed = AD_changed,
                proportionTable = proportionTable))
}

