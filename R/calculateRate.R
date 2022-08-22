calculateRate <- function(disturbanceParameters,
                          disturbanceDT,
                          disturbanceList,
                          whichToUpdate,
                          RTM,
                          overwriteDisturbanceLayers2015,
                          overwriteDisturbanceLayers2010,
                          studyArea,
                          useECCCData,
                          checkChangeInDisturbance,
                          checkDisturbance2015){
  
  dP <- disturbanceParameters[whichToUpdate, ]
  
  # Total study area
  studyAreaV <- vect(studyArea)
  studyAreaV <- project(x = studyAreaV, y = terra::crs(RTM))
  uniStudyArea <- terra::aggregate(studyAreaV)
  totalstudyAreaVAreaSqKm <- terra::expanse(uniStudyArea, unit = "km")

  updatedDisturbanceParameters <- rbindlist(lapply(1:NROW(dP), function(INDEX){
    # Get the row to be calculated
    sub <- dP[INDEX, ]
    # Get the corresponding layer
    lay <- disturbanceList[[sub[["dataName"]]]][[sub[["disturbanceOrigin"]]]]
    # If the layer was empty, it will return null. That means we don't have data to calculate the 
    # size. In this case, we will use defaults for what we know: 
    if (is.null(lay)){
      # windTurbines: 6070m2 per 2-megawatt of energy. In the NWT the only big turbine is in a
      #               Diamond exploration, with 9.2MW being generated (27,900 m2). Each pixel is
      #               62,500m2, so we will have in average, one big turbine in 2.3 pixels.
      #               There is one turbine being constructed and should be done in 2023 (2021 to be 
      #               easier on the simulations) -- Size: 3.5-megawatt = less than a pixel
      #               There are current no plans for more turbines. Still, we might be able to add
      #               one medium turbine every 10 years. That would mean 1 pixel every 10 years. 
      # There are no other values that need to be inputted. If the lay is not windTurbine, set the
      # size to NULL and basically skip it with a warning
      if (sub[["dataClass"]] == "potentialWindTurbines"){
        # 1 Turbine = 0.0625 km2
        # 0.1 Turbine = 0.00625 km2 
        # areaSqKmChangedPerYear <- 0.00625
        message(crayon::yellow(paste0("There is no information on rate for ", sub[["dataClass"]],
                                      ". However, this is a potentialWindTurbines. The module will ",
                                      "return a rate of one turbine per 10 years, or 0.1 turbines ",
                                      "per year. The 0.1 pixels (or 0.00625 km2) correspond to ", 
                                      format(100*0.00625/totalstudyAreaVAreaSqKm, 
                                             scientific = TRUE), 
                                      "% of the total area, being this the value replacing NA. ",
                                      "If this is wrong, please provide both disturbanceRate and ",
                                      "disturbanceSize in the disturbanceParameters table.")))
        sub[, disturbanceRate := 0.00625/totalstudyAreaVAreaSqKm] #areaSqKmChangedPerYear/totalAreaSqKm
      } else {
        message(crayon::red(paste0("There is no information on rate for ", sub[["dataClass"]],
                                   ". The module will return NULL and this class will not be simulated. ",
                                   "If this is wrong, please provide both disturbanceRate and disturbanceSize ",
                                   "in the disturbanceParameters table.")))
        sub <- NULL
      }
    } else {
      if (useECCCData) {
        # If the table doesn't have disturbance rate, we can consider the 0.2% of the total area a 
        # year as a default (according to ECCC report -- 2019 Species at Risk Act 
        # Conservation Agreement for the Conservation of the Boreal Caribou) or calculate it based 
        # on data!
        # 1. Get the layers (poly + lines) 2015 30m and 2010 30m
        
        dist2015Filename <- file.path(Paths$inputPath, 
                                      "totalDisturbance_2015.shp")
        dist2010Filename <- file.path(Paths$inputPath, 
                                      "totalDisturbance_2010.shp")
        if (any(!file.exists(dist2015Filename),
                !file.exists(dist2010Filename),
                overwriteDisturbanceLayers2015,
                overwriteDisturbanceLayers2010)) {
          if (any(!file.exists(dist2015Filename),
                  overwriteDisturbanceLayers2015)) {
            # 2015
            # NOTE: The 2015 layer is a GDB, and it doesn't work with prepInputs. Therefore I had to download,
            # open in ArcGIS (NorthwestTerritories_30m.gdb), create a shapefile for lines 
            # (NorthwestTerritories_30m_Disturb_Perturb_Line.shp) and one for polygons 
            # (NorthwestTerritories_30m_Disturb_Perturb_Poly.shp), zip and upload to Google Drive
            
            # Original Layer URL: "https://data-donnees.ec.gc.ca/data/species/developplans/2015-anthropogenic-disturbance-footprint-within-boreal-caribou-ranges-across-canada-as-interpreted-from-2015-landsat-satellite-imagery/Anthro_Disturb_Perturb_30m_2015.zip"
            
            lineURL <- "https://drive.google.com/file/d/14j9lvuE-Y6VnYdeERF8kLs_9gTkZKT3n/view?usp=sharing"
            AD_2015_Lines <- prepInputs(url = lineURL,
                                        archive = "NorthwestTerritories_30m_Disturb_Perturb_Line.zip",
                                        alsoExtract = "similar",
                                        studyArea = studyArea,
                                        targetFile = "NorthwestTerritories_30m_Disturb_Perturb_Line.shp",
                                        destinationPath = Paths$inputPath)
            AD_2015_Lines <- vect(AD_2015_Lines)
            AD_2015_Lines <- project(AD_2015_Lines, terra::crs(studyArea)) 
            # Although SA likely has the same projection as the rest of the data, it might not. 
            # So I need to reproject everything tothe study area to match all other layers
            
            polyURL <- "https://drive.google.com/file/d/1bbgRZjLXZzo-jJgdMccqwIwFd5dp1LDt/view?usp=sharing"
            AD_2015_Polys <- prepInputs(url = polyURL,
                                        alsoExtract = "similar",
                                        studyArea = studyArea,
                                        targetFile = "NorthwestTerritories_30m_Disturb_Perturb_Poly.shp",
                                        destinationPath = Paths$inputPath)
            AD_2015_Polys <- vect(AD_2015_Polys)
            AD_2015_Polys <- project(AD_2015_Polys, terra::crs(studyArea))
            
            # 2. Clip all 4 layers to the SA
            AD_2015_Polys_clip <- terra::intersect(x = AD_2015_Polys, y = uniStudyArea)
            AD_2015_Lines_clip <- terra::intersect(x = AD_2015_Lines, y = uniStudyArea)
            # 3. Buffer both lines to 30m
            AD_2015_Lines_clipBuff <- terra::buffer(x = AD_2015_Lines_clip, width = 30)
            # 4. Union of polys and buffered lines for both 2015 and 2010, merging overlapping features by class
            # 4.1. First, though, we need to simplify the AD_2015_Lines_clipBuff to only the columns in AD_2015_Polys
            # 2015 --> NOPE.
            # As it turns out, this process is highly time consuming and the classification ends up being impossible
            # We decided instead to just:
            # 1. Clip out (mask with inverse = TRUE to keep what does NOT intersect) the polygons from the buffered lines  
            AD_2015 <- terra::mask(AD_2015_Lines_clipBuff, mask = AD_2015_Polys_clip, inverse = TRUE)
            # 2. Simply rbind the two shapefiles (polys + buffered lines with polys cropped out). Do the same for 2010.
            AD_2015_all <- rbind(AD_2015, AD_2015_Polys_clip)
            # 3. Calculate the rate of growth for each of the 12 different classes
            # 3.1. Calculate the area for each individual form:
            AD_2015_all$Area_sqKm <- expanse(AD_2015_all) / 1000000
            # Now that the calculations are done, save the objects
            writeVector(x = AD_2015_all, filename = dist2015Filename, overwrite = TRUE)
            
          } else {
            # If objects exist, just load them
            AD_2015_all <- vect(dist2015Filename)
          }
          if (any(!file.exists(dist2010Filename),
                  overwriteDisturbanceLayers2010)) {
            # 2010
            url_2010 <- "https://www.ec.gc.ca/data_donnees/STB-DGST/003/Boreal-ecosystem-anthropogenic-disturbance-vector-data-2008-2010.zip"
            AD_2010_Lines <- prepInputs(url = url_2010,
                                        archive = "Boreal-ecosystem-anthropogenic-disturbance-vector-data-2008-2010.zip",
                                        alsoExtract = "similar",
                                        studyArea = studyArea,
                                        targetFile = "EC_borealdisturbance_linear_2008_2010_FINAL_ALBERS.shp",
                                        destinationPath = Paths$inputPath)
            AD_2010_Lines <- vect(AD_2010_Lines)
            AD_2010_Lines <- project(AD_2010_Lines, terra::crs(studyArea))
            
            AD_2010_Polys <- prepInputs(url = url_2010,
                                        archive = "Boreal-ecosystem-anthropogenic-disturbance-vector-data-2008-2010.zip",
                                        alsoExtract = "similar",
                                        studyArea = studyArea,
                                        targetFile = "EC_borealdisturbance_polygonal_2008_2010_FINAL_ALBERS.shp",
                                        destinationPath = Paths$inputPath)
            AD_2010_Polys <- vect(AD_2010_Polys)
            AD_2010_Polys <- project(AD_2010_Polys, terra::crs(studyArea))
            
            # 2. Clip all 4 layers to the SA
            AD_2010_Polys_clip <- terra::intersect(x = AD_2010_Polys, y = uniStudyArea)
            AD_2010_Lines_clip <- terra::intersect(x = AD_2010_Lines, y = uniStudyArea)
            
            # 3. Buffer both lines to 30m
            AD_2010_Lines_clipBuff <- terra::buffer(x = AD_2010_Lines_clip, width = 30)
            
            # 4. Union of polys and buffered lines for both 2015 and 2010, merging overlapping features by class
            # 4.1. First, though, we need to simplify the AD_2015_Lines_clipBuff to only the columns in AD_2015_Polys
            # 2015 --> NOPE.
            # As it turns out, this process is highly time consuming and the classification ends up being impossible
            # We decided instead to just:
            # 1. Clip out (mask with inverse = TRUE to keep what does NOT intersect) the polygons from the buffered lines  
            AD_2010 <- terra::mask(AD_2010_Lines_clipBuff, mask = AD_2010_Polys_clip, inverse = TRUE)
            
            # 2. Simply rbind the two shapefiles (polys + buffered lines with polys cropped out). Do the same for 2010.
            AD_2010_all <- rbind(AD_2010, AD_2010_Polys_clip)
            
            # 3. Calculate the rate of growth for each of the 12 different classes
            # 3.1. Calculate the area for each individual form:
            AD_2010_all$Area_sqKm <- expanse(AD_2010_all) / 1000000
            
            # Now that the calculations are done, save the objects
            writeVector(x = AD_2010_all, filename = dist2010Filename, overwrite = TRUE)
          } else {
            # If objects exist, just load them
            AD_2010_all <- vect(dist2010Filename)
          }
        } else {
          # If objects exist, just load them
          AD_2015_all <- vect(dist2015Filename)
          AD_2010_all <- vect(dist2010Filename)
        }
          # 3.2. Bring all to a data.table and summarize
          AD_2015_DT <- as.data.table(as.data.frame(AD_2015_all[, c("Class", "Area_sqKm")]))
          AD_2010_DT <- as.data.table(as.data.frame(AD_2010_all[, c("Class", "Area_sqKm")]))
          # Make sure all classes match from 2010 to 2015
          testthat::expect_true(all(sort(unique(AD_2015_DT$Class)) == sort(unique(AD_2010_DT$Class))))
          AD_2015_DT_summ <- AD_2015_DT[, totalArea_sqKm := sum(Area_sqKm), by = "Class"]
          AD_2015_DT_summ <- unique(AD_2015_DT_summ[, Area_sqKm := NULL])
          AD_2010_DT_summ <- AD_2010_DT[, totalArea_sqKm := sum(Area_sqKm), by = "Class"]
          AD_2010_DT_summ <- unique(AD_2010_DT_summ[, Area_sqKm := NULL])
          
          # 3.3. Now I can see how much the disturbances changed from 2010 to 2015 in sq Km
          AD_2015_DT_summ[, Year := "year2015"]
          AD_2010_DT_summ[, Year := "year2010"]
          AD_change <- rbind(AD_2015_DT_summ, AD_2010_DT_summ)
          AD_changed <- dcast(data = AD_change, formula = Class ~ Year, value.var = "totalArea_sqKm")

          if (checkChangeInDisturbance){
            # EXTRA: Quick calculation of the % disturbance change per year over the entire uniStudyArea area:
            totAreaChangedPerYearSqKm <- sum(AD_changed$areaSqKmChangedPerYear)
            print(paste0("Percentage of area changed per year averaged across 2010 to 2015 across ",
                         "the total area: ", 
                         round((totAreaChangedPerYearSqKm/totalstudyAreaVAreaSqKm)*100, 3), "%"))
          }
          if (checkDisturbance2015){
            # EXTRA: Quick calculation of the % 500m buffered disturbance over the entire studyAreaV area:
            AD_2015_L_500m <- buffer(x = AD_2015_Lines_clip, width = 500)
            AD_2015_P_500m <- buffer(x = AD_2015_Polys_clip, width = 500)
            AD_2015_L_500m <- mask(AD_2015_L_500m, mask = AD_2015_P_500m, inverse = TRUE)
            AD_2015_500m_all <- rbind(AD_2015_L_500m, AD_2015_P_500m)
            AD_2015_500m_all$Area_sqKm <- expanse(AD_2015_500m_all) / 1000000
            AD_2015_500m_DT <- as.data.table(as.data.frame(AD_2015_500m_all[, c("Class", "Area_sqKm")]))
            totDistKm2 <- sum(AD_2015_500m_DT$Area_sqKm)
            print(paste0("Percentage of 500m buffered human disturbance in the area based on 2015 data ",
                         " is: ", 
                         round((totDistKm2/totalstudyAreaVAreaSqKm)*100, 3), "%"))
          }
          
          # 4. Aggregate disturbances that need to be aggregated (based on the table I have)
          # 4.1. The disturbances that need rate are the Enlarging and the Generating ones. Subset 
          # what I have to those
          classesAvailable <- unique(disturbanceDT[fieldToSearch == "Class", c("classToSearch", "dataName", "dataClass")])
          toLookFor <- classesAvailable[dataClass %in% dP[["disturbanceOrigin"]]]
          toFill <- AD_changed[Class %in% toLookFor[["classToSearch"]]]
          toUse <- merge(toFill, toLookFor[, c("classToSearch", "dataClass")], all.x = TRUE, 
                         by.x = "Class", by.y = "classToSearch")
          # Sum the disturbances that are the same
          toUseUp <- toUse[, c("year2010", "year2015") := list(sum(year2010), 
                                                               sum(year2015)), 
                           by = "dataClass"] 
          toUseUp[, Class := NULL]
          toUseUp <- unique(toUseUp)
          # Add the total study area area so I can calculate the rate over the total study area that we
          # need to generate
          toUseUp[, totalAreaSqKm := totalstudyAreaVAreaSqKm]
          
          toUseUp[, areaSqKmChangedPerYear := (year2015-year2010)/5]
          toUseUp[, percOfTotalAreaToIncreasePerYear := areaSqKmChangedPerYear/totalAreaSqKm]
          # If any have negative growth, we ignore
          toUseUp[percOfTotalAreaToIncreasePerYear < 0, percOfTotalAreaToIncreasePerYear := 0]
          
          # 5. Modify the table, which should be an Input.
          # This process ensures we are not double counting disturbances, and that we can calculate the change 
          # in disturbance per class, which we need to simulate the new disturbances coming.
          toUseUpSub <- toUseUp[dataClass %in% sub[["disturbanceOrigin"]], ]
          if (NROW(toUseUpSub) > 1) {
            stop("Repeated data found in the process. Please debug")
          }
          sub[, disturbanceRate := toUseUpSub[["percOfTotalAreaToIncreasePerYear"]]]
        } else {
        sub[, disturbanceRate := 0.0001790889] # Here we default to 0.0001790889 over the total area, which 
                                            # corresponds to approximately 0.2% of the current 
                                            # disturbance, as mentioned by ECCC 2019 caribou report
                                            # Using our data, this value seems to be overstated
        }
      }
    return(sub)
  }))
  
  disturbanceParametersUp <- rbind(disturbanceParameters[!whichToUpdate, ], 
                                   updatedDisturbanceParameters, use.names = TRUE)
  return(disturbanceParametersUp)
  
}