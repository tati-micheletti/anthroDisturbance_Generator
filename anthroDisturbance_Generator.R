## Everything in this file and any files in the R directory are sourced during `simInit()`;
## all functions and objects are put into the `simList`.
## To use objects, use `sim$xxx` (they are globally available to all modules).
## Functions can be used inside any function that was sourced in this module;
## they are namespaced to the module, just like functions in R packages.
## If exact location is required, functions will be: `sim$.mods$<moduleName>$FunctionName`.
defineModule(sim, list(
  name = "anthroDisturbance_Generator",
  description = paste0("This is a module to generate anthropogenic disturbances.",
                       "It's primarily intended for the Northwest Territories region",
                       " (default), but the structure is universal.",
                       " All needed to do is provide the metadata information ",
                       "required by the `disturbanceParameters` ",
                       "object and the disturbanceList containing current disturbances ",
                       "and potential disturbances (i.e., for disturbances of Generating type)."),
  keywords = "",
  authors = structure(list(list(given = "Tati", 
                                family = "Micheletti", role = c("aut", "cre"), 
                                email = "tati.micheletti@gmail.com", 
                                comment = NULL)), 
                      class = "person"),
  childModules = character(0),
  version = list(anthroDisturbance_Generator = "0.0.2"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "anthroDisturbance_Generator.Rmd"), ## same file
  reqdPkgs = list("SpaDES.core (>=1.0.10)", "ggplot2", 
                  "data.table", "PredictiveEcology/reproducible",
                  "raster", "terra", "crayon", "msm", "sf", "pik-piam/rmndt",
                  "fasterize", "stars", "nngeo", "tictoc"), #TODO review needed packages. 
  parameters = rbind(
    defineParameter(".plots", "character", "screen", NA, NA,
                    "Used by Plots function, which can be optionally used here"),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA,
                    "Describes the simulation time at which the first plot event should occur."),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "Describes the simulation time interval between plot events."),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "Describes the simulation time at which the first save event should occur."),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events."),
    defineParameter(".seed", "list", list(), NA, NA,
                    "Named list of seeds to use for each event (names)."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Should caching of events or module be used?"),
    defineParameter("runInterval", "numeric", 1, NA, NA,
                    paste0("Should the module be run every decade? This speeds up module testing as ",
                           "testing if the events need to be run at every time is time-consuming. If ",
                           "the user knows the disturbances happen every X years, X can be passed here.")),
    defineParameter("saveInitialDisturbances", "logical", TRUE, NA, NA,
                    paste0("Should the disturbance rasters be saved at each step? These are saved ",
                           "to Paths[['outputPath']] as a RasterLayer, with disturbanceLayer as prefix",
                           "the name of the industry and the year as suffix.",
                           "If TRUE, it saves the initial conditions (IC)")),
    defineParameter("generatedDisturbanceAsRaster", "logical", FALSE, NA, NA,
                    paste0("Should the new disturbances generated be in raster format? This has ",
                           "potential downsides regarding size of disturbances generated (i.e., minimum",
                           " size possible is the resolution of raster)")),
    defineParameter("checkDisturbancesForBuffer", "logical", FALSE, NA, NA,
                    paste0("Should the module check the recently generated disturbances? This means ",
                           "that the module will buffer the recently disturbed layer and compare it ",
                           "to the previous layer, outputting a message indicating how much of the total",
                           " area this represents. Note that no objects are created, just the message ",
                           "is outputted.")),
    defineParameter("disturbFirstYear", "logical", FALSE, NA, NA,
                    paste0("Should disturbances be generated already in the initial year? Normally, ",
                           "we would save the initial disturbances (i.e., 2011) as they are coming from ",
                           "data, and only start generating disturbances once we don't have the data",
                           "(i.e., post-2015). So this defaults to FALSE.",
                           "If TRUE, it will already generate disturbances in start(sim)")),
    defineParameter("saveCurrentDisturbances", "logical", TRUE, NA, NA,
                    paste0("Should the disturbance rasters be saved at each step? These are saved ",
                           "to Paths[['outputPath']] as a RasterLayer, with disturbanceLayer as prefix",
                           "the name of the industry and the year as suffix.",
                           "If TRUE, it saves at the end of each step.")),
    defineParameter("disturbanceRateRelatesToBufferedArea", "logical", TRUE, NA, NA,
                    paste0("Is the DisturbanceRate a % of already buffered (to 500m) disturbance?",
                           " This is normally what is used for caribou.")),
    defineParameter("growthStepEnlargingPolys", "numeric", 1, NA, NA,
                    paste0("Growth step used for iteratively achieving the total area growth of ",
                           "new disturbances type Enlarging for polygons. If the iterations take too",
                           " long, one should increase this number. If the summarized value is too",
                           " far from 0, one should decrease this number.")),
    defineParameter("growthStepGenerating", "numeric", 1, NA, NA,
                    paste0("Increasing factor to speed up total area of ",
                           "new disturbances type Generating. If the iterations take too",
                           " long, one should increase this number. If the summarized value is too",
                           " far from 0, one should decrease this number.",
                           " Not used if disturbanceRateRelatesToBufferedArea == TRUE")),
    defineParameter("growthStepEnlargingLines", "numeric", 1, NA, NA,
                    paste0("Growth step used for iteratively achieving the total area growth of ",
                           "new disturbances type Enlarging for lines. If the iterations take too",
                           " long, one should increase this number. If the summarized value is too",
                           " far from 0, one should decrease this number.")),
    defineParameter("connectingBlockSize", "numeric", NULL, NA, NA,
                    paste0("connectingBlockSize is default to NULL. It is used to connecting layers ",
                           "after generation. Applying blocking technique speeds up disturbance.",
                           " If too high, many lines might connect from the same place. Decreasing ",
                           " the parameter connectingBlockSize or setting it to NULL, will improve")),
    defineParameter(".runName", "character", "run1", NA, NA,
                    paste0("If you would like your simulations' results to have an appended name ",
                           "(i.e., replicate number, study area, etc) you can use this parameter")),
    defineParameter(".inputFolderFireLayer", "character", Paths[["inputPath"]], NA, NA,
                    paste0("If you have the fire (i.e., rstCurrBurn) in a folder that is NOT the ",
                           "inputs folder, you can pass it here")),
    defineParameter("totalDisturbanceRate", "numeric", NULL, NA, NA,
                    paste0("If passed, the module will use ECCC data to calculate the % each",
                           "disturbance should represent to be to achieve (as close as possible)",
                           " the total expected disturbance rate. Mainly used for early in ",
                           "simulations, when value ranges are unknown.")),
    defineParameter("aggregateSameDisturbances", "logical", TRUE, NA, NA,
                    paste0("If TRUE, when using ECCC data to calculate disturbance rates, it aggregates",
                           " the features that have the same class and dissolves their boundaries. ",
                           " This may influence disturbance rate calculations especially if ",
                           "disturbanceRateRelatesToBufferedArea = TRUE (i.e., overlapping disturbances ",
                           "from the same class will be double counted if aggregateSameDisturbances = FALSE) ",
                           "If DisturbanceRate is provided, this parameter is ignored.")),
    defineParameter("maskOutLinesFromPolys", "logical", TRUE, NA, NA,
                    paste0("If TRUE, when using ECCC data to calculate disturbance rates, it masks out",
                           " the lines from polygons when these overlap (i.e., more likely when ",
                           " disturbanceRateRelatesToBufferedArea = TRUE). This may influence ",
                           "disturbance rate calculations as these will be double counted if FALSE).",
                           " If DisturbanceRate is provided, this parameter is ignored."))
  ),
  inputObjects = bindrows(
    #expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
    expectsInput(objectName = "disturbanceList", objectClass = "list",
                 desc = paste0("List (general category) of lists (specific ",
                               "class) needed for generating ",
                               "disturbances. This is generally the output from a potential",
                               "Resources module (i.e., potentialResourcesNT_DataPrep), ",
                               "where multiple potential layers (i.e., mining",
                               " and oilGas) we replaced by only one layer with the highest",
                               "values being the ones that need to be filled ",
                               "with new developments first, or prepared potential layers",
                               " (i.e., potentialCutblocks)."),
                 sourceURL = "https://drive.google.com/file/d/1v7MpENdhspkWxHPZMlmx9UPCGFYGbbYm/view?usp=sharing"),
    
    expectsInput(objectName = "disturbanceParameters", objectClass = "data.table",
                 desc = paste0("Table with the following columns: ",
                               
                               "dataName --> this column groups the type of data ",
                               "by sector (i.e., Energy, Settlements, OilGas, ",
                               "Mining, Forestry, Roads)",
                               
                               "dataClass --> this column details the type of data ",
                               "ALWAYS with 'potential' starting (i.e., potentialSettlements ",
                               "potentialWindTurbines, potentialCutblocks, etc.)",
                               "can harmonize different ones. Potential data classes ",
                               "can be of three general disturbanceType (see below)",
                               
                               "disturbanceType --> Potential data classes ",
                               "can be of three general types:",
                               "1. Enlarging (i.e., potentialSettlements",
                               " and potentialSeismicLines): where the potential one is ",
                               "exactly the same as the current layer, and we",
                               " only buffer it with time",
                               "2. Generating (i.e., potentialWindTurbines, potentialOilGas",
                               "potentialMineral, potentialForestry): where the ",
                               "potential layers are only the potential where ",
                               "structures can appear based on a specific rate",
                               "3. Connecting (i.e., potentialPipelines, ",
                               "potentialTransmission, potentialRoads incl. ",
                               "forestry ones): where the potential layer needs ",
                               "to have the current/latest transmission, pipeline,",
                               " and road network. This process will depend on ",
                               "what is generated in point 2.",
                               
                               "disturbanceRate --> what is the rate of generation for ",
                               "disturbances PER YEAR in % of type Enlarging and Generating. For ",
                               "disturbances type Connecting, disturbanceRate is NA. ",
                               "If not specified when needed, the module will try to derive ",
                               "it from data.",
                               
                               "disturbanceSize --> if there is a specific size the disturbance in m2 ",
                               "type Generating should have, it is specified here. If not specified, ",
                               " the module will try to derive it from data. For disturbances ",
                               "type Enlarging anb Connecting, disturbanceSize is NA.",
                               "If not specified when needed, the module will try to derive ",
                               "it from data.",
                               
                               "disturbanceOrigin --> dataClass that should be used as the 'original' ",
                               "to be either modified (i.e., Enlarging, Generating) or as origin  ",
                               "point for Connecting types.",
                               
                               "disturbanceEnd --> end points for Connecting layers (i.e., ",
                               "newly created windTurbines: connect into powerLines, ",
                               "newly created windTurbines: connect into roads, ",
                               "newly created oilGas: connect into pipeline, ",
                               "newly created oilGas: connect into roads, ",
                               "newly created settlements: connect into roads, ",
                               "newly created mines: connect into roads",
                               "newly created cutblocks: connect into roads)",
                               
                               "disturbanceInterval --> interval for which ",
                               "this disturbance should happen ",
                               
                               "resolutionVector --> original resolution ",
                               "of the data that generated this vector ",
                               
                               "It defaults to an example in the Northwest ",
                               "Territories and needs to be provided if the ",
                               "study area is not in this region (i.e., union ",
                               "of BCR6 and NT1)"),
                 sourceURL = "https://drive.google.com/file/d/1xJypz-VOA_bHN0y4GYikY25UKID_oTh_/view?usp=sharing"),
    expectsInput(objectName = "disturbanceDT", objectClass = "data.table", 
                 desc = paste0("This data.table needs to contain the following",
                               " columns: ",
                               
                               "dataName --> this column groups the type of data ",
                               "by sector (i.e., Energy, Settlements, OilGas, ",
                               "Mining, Forestry, Roads)",
                               
                               "URL --> URL link for the specific dataset",
                               
                               "classToSearch --> exact polygon type/class to ",
                               "search for when picking from a dataset with ",
                               "multiple types. If this is not used (i.e., your",
                               " shapefile is alreday all the data needed), you ",
                               "should still specify this so each entry has a ",
                               "different name",
                               
                               "fieldToSearch --> where should classToSearch be ",
                               "found? If this is specified, then the function ",
                               "will subset the spatial object (most likely a ",
                               "shapefile) to classToSearch. Only provide this ",
                               "if this is necessary!",
                               
                               "dataClass --> this column details the type of data ",
                               "further (i.e., Settlements, potentialSettlements ",
                               "otherPolygons, otherLines, windTurbines, potentialWindTurbines, ",
                               "hydroStations, oilFacilities, pipelines, etc). ",
                               "Common class to rename the dataset to, so we ",
                               "can harmonize different ones. Potential data classes ",
                               "can be of three general types (that will be ",
                               "specified in the disturbanceGenerator module as ",
                               "a parameter -- ALWAYS with 'potential' starting): ",
                               "1. Enlarging (i.e., potentialSettlements",
                               " and potentialSeismicLines): where the potential one is ",
                               "exactly the same as the current layer, and we",
                               " only buffer it with time",
                               "2. Generating (i.e., potentialWind, potentialOilGas",
                               "potentialMineral, potentialForestry): where the ",
                               "potential layers are only the potential where ",
                               "structures can appear based on a specific rate",
                               "3. Connecting (i.e., potentialPipelines, ",
                               "potentialRoads incl. ",
                               "forestry ones): where the potential layer needs ",
                               "to have the current/latest transmission, pipeline,",
                               " and road network. This process will depend on ",
                               "what is generated in point 2.",
                               
                               "fileName --> If the original file is a .zip and ",
                               "the features are stored in one of more shapefiles",
                               " inside the .zip, please provide which shapefile ",
                               " to be used",
                               
                               "dataType --> please provide the data type of the ",
                               "layer to be used. These are the current accepted ",
                               " formats: 'shapefile' (.shp or .gdb), 'raster' ",
                               "(.tif, which will be converted into shapefile), ",
                               "and 'mif' (which will be read as a shapefile).",
                               
                               "It defaults to an example in the Northwest ",
                               "Territories and needs to be provided if the ",
                               "study area is not in this region (i.e., union ",
                               "of BCR6 and NT1)"),
                 sourceURL = "https://drive.google.com/file/d/1wHIz_G088T66ygLK9i89NJGuwO3f6oIu/view?usp=sharing"),
    expectsInput(objectName = "studyArea", 
                 objectClass = "SpatialPolygonDataFrame|vect", 
                 desc = paste0("Study area to which the module should be ",
                               "constrained to. Defaults to NT1+BCR6. Object ",
                               "can be of class 'vect' from terra package"), 
                 sourceURL = "https://drive.google.com/file/d/1RPfDeHujm-rUHGjmVs6oYjLKOKDF0x09/view?usp=sharing"),
    expectsInput(objectName = "rasterToMatch", 
                 objectClass = "RasterLayer|rast", 
                 desc = paste0("All spatial outputs will be reprojected and ",
                               "resampled to it. Defaults to NT1+BCR6. Object ",
                               "can be of class 'rast' from terra package"), 
                 sourceURL = "https://drive.google.com/file/d/11yCDc2_Wia2iw_kz0f0jOXrLpL8of2oM/view?usp=sharing"),
    expectsInput(objectName = "rstCurrentBurn", 
                 objectClass = "RasterLayer", 
                 desc = paste0("A binary raster with 1 values representing burned pixels. ",
                               "This raster is normally produced by either the module historicFires or ",
                               "a fire simulation module (i.e., fireSense, SCFM, LandMine)"),
                 sourceURL = NA),
    expectsInput(objectName = "DisturbanceRate", 
                 objectClass = "data.table", 
                 desc = paste0("Rate of change (disturbance) increase over the study area per year. ",
                               "Defaults to calculating the disturbance (ECCC data) over the entire area per year,",
                               " if totalDisturbanceRate is not provided. Needs to have:",
                               "dataName: settlements, oilGas, oilGas, mining, forestry, Energy",
                               "dataClass: potentialSettlements, potentialSeismicLines, potentialOilGas,",
                               "potentialMining, potentialCutblocks, potentialWindTurbines",
                               "disturbanceType: Enlarging or Generating (see disturbanceDT object for details)",
                               "disturbanceOrigin: settlements, seismicLines, oilGas, mining, cutblocks, windTurbines",
                               "disturbanceRate: representing a % of the study area to be newly disturbed per year"),
                 sourceURL = NA)
    ),
  outputObjects = bindrows(
    createsOutput(objectName = "disturbanceList", objectClass = NA, 
                  desc = paste0("Updated list (general category) of lists (specific ",
                                "class) of disturbances and the potential needed for ",
                                "generating disturbances. ")),
    createsOutput(objectName = "currentDisturbanceLayer", objectClass = "list", 
                  desc = paste0("List (per year) of rasters with all current disturbances.",
                                "Can be used  for other purposes but was created to filter potential",
                                " pixels that already have disturbances to avoid choosing new ",
                                "pixels in existing disturbed ones"))
  )
))

## event types
#   - type `init` is required for initialization

doEvent.anthroDisturbance_Generator = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      
      # Make sure RTM and studyArea projections match
      sim$studyArea <- projectInputs(sim$studyArea, 
                                     targetCRS = crs(sim$rasterToMatch))

      ### Make sure that disturbanceDT is data.table      
      if (any(class(sim$disturbanceDT) != "data.table")){
        tryCatch({
          sim$disturbanceDT <- as.data.table(sim$disturbanceDT)
        }, error = function(e){
          warning(paste0("disturbanceDT was provided with class ",
                         class(sim$disturbanceDT), " and conversion to ",
                         "data.table failed. Please provide this object as a ",
                         "data.table"), immediate. = TRUE)
          stop(e)
        })
      }
      
      sim$currentDisturbanceLayer <- list()
      
      # schedule future event(s)
        sim <- scheduleEvent(sim, time(sim), "anthroDisturbance_Generator", "calculatingSize", eventPriority = 4)
        sim <- scheduleEvent(sim, time(sim), "anthroDisturbance_Generator", "calculatingRate", eventPriority = 4)
        sim <- scheduleEvent(sim, time(sim), "anthroDisturbance_Generator", "generatingDisturbances", eventPriority = 4)
        sim <- scheduleEvent(sim, time(sim), "anthroDisturbance_Generator", "updatingDisturbanceList", eventPriority = 4)
    },
    calculatingSize = {
      # Check for needed sizes. If all provided, skip event
      mod$.whichNeedSize <- which(sim$disturbanceParameters[, disturbanceType == "Generating" &
                                                              is.na(disturbanceSize)])
      if (length(mod$.whichNeedSize) != 0)
        sim$disturbanceParameters <- calculateSize(disturbanceParameters = sim$disturbanceParameters,
                                                 whichToUpdate = mod$.whichNeedSize,
                                                 disturbanceList = sim$disturbanceList)
    },
    calculatingRate = {
      # Check for needed rates. If all provided, skip event
      mod$.whichNeedRates <- which(sim$disturbanceParameters[, disturbanceType %in% c("Generating", "Enlarging") &
                                                               is.na(disturbanceRate)])
      if (length(mod$.whichNeedRates) != 0)
        sim$disturbanceParameters <- calculateRate(disturbanceParameters = sim$disturbanceParameters,
                                                 whichToUpdate = mod$.whichNeedRates,
                                                 studyArea = sim$studyArea,
                                                 disturbanceList = sim$disturbanceList,
                                                 RTM = sim$rasterToMatch,
                                                 disturbanceDT = sim$disturbanceDT,
                                                 destinationPath = Paths[["inputPath"]],
                                                 totalDisturbanceRate = P(sim)$totalDisturbanceRate,
                                                 DisturbanceRate = sim$DisturbanceRate,
                                                 disturbanceRateRelatesToBufferedArea = P(sim)$disturbanceRateRelatesToBufferedArea,
                                                 maskOutLinesFromPolys = P(sim)$maskOutLinesFromPolys,
                                                 aggregateSameDisturbances = P(sim)$aggregateSameDisturbances)
      },
    generatingDisturbances = {

      if (all(P(sim)$saveInitialDisturbances,
              start(sim) == time(sim))){
        message(crayon::yellow(paste0("The parameter saveInitialDisturbances is TRUE.",
                                      " Saving initial disturbance layers")))
        saveDisturbances(disturbanceList = sim$disturbanceList,
                         currentTime = "IC", 
                         overwrite = TRUE,
                         runName = P(sim)$.runName)
        
        initialBufferedAnthropogenicDisturbance500m <- createBufferedDisturbances(disturbanceList = sim$disturbanceList,
                                                                                  bufferSize = 500,
                                                                                  rasterToMatch = sim$rasterToMatch,
                                                                                  studyArea = sim$studyArea,
                                                                                  currentTime = "IC",
                                                                                  convertToRaster = TRUE)
        anthroDistFilePath <- file.path(Paths[["outputPath"]],
                                        paste0("bufferedAnthDist_500m_", 
                                               time(sim), ".tif"))
        
        message(paste0("Writing buffered disturbance layer for ", time(sim)))
        terra::writeRaster(initialBufferedAnthropogenicDisturbance500m, 
                           filename = anthroDistFilePath,
                           overwrite = TRUE)
      }
      
      # Check if the time(sim) is within the interval to run the disturbances
      mod$.whichToRun <- if (all(!P(sim)$disturbFirstYear,
                                 time(sim) == start(sim))){
        NULL
      } else {
        whichDisturbancesToGenerate(startTime = start(sim),
                                    currentTime = time(sim),
                                    endTime = end(sim),
                                    disturbanceParameters = sim$disturbanceParameters)
      }
      
      if (length(mod$.whichToRun) != 0){ # If anything is scheduled
        forestryScheduled <- "forestry" %in% sim$disturbanceParameters[mod$.whichToRun, dataName]
        if (forestryScheduled) { # forestry scheduled?
          if (is.null(sim$rstCurrentBurn)){
            # Try to get rstCurrentBurn if no module is producing it!
            mod$rstCurrentBurn <- createModObject(data = "rstCurrentBurn", 
                                                  sim = sim, 
                                                  pathInput = P(sim)$.inputFolderFireLayer, 
                                                  currentTime = time(sim),
                                                  returnNULL = TRUE,
                                                  fun = raster::raster)
            if (!is.null(mod$rstCurrentBurn))
              message(crayon::green("Fire layer found in inputs folder! Using it for the simulation."))
          } else {
            mod$rstCurrentBurn <- sim$rstCurrentBurn
          }
        }
        currDis <- tryCatch({
          tail(sim$currentDisturbanceLayer,
               n = 1)[[1]]
        }, error = function(e){
          return(NULL) 
        })
        if (P(sim)$generatedDisturbanceAsRaster){
          mod$updatedLayers <- generateDisturbances(disturbanceParameters = sim$disturbanceParameters,
                                                       disturbanceList = sim$disturbanceList,
                                                       currentTime = time(sim),
                                                       studyArea = sim$studyArea,
                                                       rasterToMatch = sim$rasterToMatch,
                                                       fires = mod$rstCurrentBurn,
                                                       currentDisturbanceLayer = currDis,
                                                       disturbanceRateRelatesToBufferedArea = P(sim)$disturbanceRateRelatesToBufferedArea,
                                                       growthStepEnlargingPolys = P(sim)$growthStepEnlargingPolys,
                                                       growthStepEnlargingLines = P(sim)$growthStepEnlargingLines,
                                                       growthStepGenerating = P(sim)$growthStepGenerating,
                                                       connectingBlockSize = P(sim)$connectingBlockSize,
                                                       outputsFolder = Paths[["outputPath"]],
                                                       runName = P(sim)$.runName,
                                                       checkDisturbancesForBuffer = P(sim)$checkDisturbancesForBuffer)
        } else {
          mod$updatedLayers <- generateDisturbancesShp(disturbanceParameters = sim$disturbanceParameters,
                                                       disturbanceList = sim$disturbanceList,
                                                       currentTime = time(sim),
                                                       studyArea = sim$studyArea,
                                                       rasterToMatch = sim$rasterToMatch,
                                                       fires = mod$rstCurrentBurn,
                                                       currentDisturbanceLayer = currDis,
                                                       disturbanceRateRelatesToBufferedArea = P(sim)$disturbanceRateRelatesToBufferedArea,
                                                       growthStepEnlargingPolys = P(sim)$growthStepEnlargingPolys,
                                                       growthStepEnlargingLines = P(sim)$growthStepEnlargingLines,
                                                       growthStepGenerating = P(sim)$growthStepGenerating,
                                                       connectingBlockSize = P(sim)$connectingBlockSize,
                                                       outputsFolder = Paths[["outputPath"]],
                                                       runName = P(sim)$.runName,
                                                       checkDisturbancesForBuffer = P(sim)$checkDisturbancesForBuffer)
        }

        sim$currentDisturbanceLayer[[paste0("Year", time(sim))]] <- mod$updatedLayers$currentDisturbanceLayer
      } else {
        message("No disturbances scheduled for this year")
      }
      sim <- scheduleEvent(sim, time(sim) + P(sim)$runInterval, "anthroDisturbance_Generator", "generatingDisturbances")
    },
    updatingDisturbanceList = {
      if (length(mod$.whichToRun) != 0){
      sim$disturbanceList <- replaceList(disturbanceList = sim$disturbanceList, 
                                         updatedLayers = mod$updatedLayers$individuaLayers)
      
      if (P(sim)$saveCurrentDisturbances){
        message(paste0("Saving current disturbances for year ", time(sim)))
        saveDisturbances(disturbanceList = sim$disturbanceList,
                         currentTime = time(sim), 
                         overwrite = TRUE,
                         runName = P(sim)$.runName)
        }
      }
      sim <- scheduleEvent(sim, time(sim) + P(sim)$runInterval, "anthroDisturbance_Generator", "updatingDisturbanceList")

    },
    warning(paste("Undefined event type: \'", current(sim)[1, "eventType", with = FALSE],
                  "\' in module \'", current(sim)[1, "moduleName", with = FALSE], "\'", sep = ""))
  )
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  dPath <- dataPath(sim)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")
  if (!suppliedElsewhere(object = "studyArea", sim = sim)) {
    sim$studyArea <- prepInputs(url = extractURL("studyArea"),
                                targetFile = "NT1_BCR6.shp",
                                alsoExtract = "similar",
                                destinationPath = dPath)
    
    message(crayon::red(paste0("studyArea was not supplied. Defaulting to BCR6+NT1 in the",
                               " Northwest Territories")))
  }
  
  if (!suppliedElsewhere(object = "rasterToMatch", sim = sim)) {
    sim$rasterToMatch <- prepInputs(url = extractURL("rasterToMatch"),
                                    targetFile = "RTM.tif",
                                    destinationPath = dPath)
    
    message(crayon::red(paste0("rasterToMatch was not supplied. Defaulting to BCR6+NT1 in the",
                               " Northwest Territories")))
  }
  
  if (!suppliedElsewhere(object = "disturbanceList", sim = sim)) {
    sim$disturbanceList <- unwrapTerraList(terraList = extractURL("disturbanceList"), 
                                           generalPath = dataPath(sim))
    
    message(crayon::red(paste0("disturbanceList was not supplied. The current should only ",
                               " be used for module testing purposes ! Please run the module(s) ",
                               "`anthroDisturbance_DataPrep` and `potentialResourcesNT_DataPrep`")))
  }
  
  if (!suppliedElsewhere(object = "disturbanceParameters", sim = sim)) {
    sim$disturbanceParameters <- prepInputs(url = extractURL("disturbanceParameters"),
                                            targetFile = "disturbanceParameters.csv",
                                            destinationPath = dPath,
                                            fun = "data.table::fread",
                                            header = TRUE, 
                                            userTags = "disturbanceParameters")

    message(crayon::red(paste0("disturbanceParameters was not supplied. Defaulting to an example from ",
                   " Northwest Territories")))
  }
  
  if (!suppliedElsewhere(object = "rstCurrentBurn", sim = sim)) {
    
    message(crayon::red(paste0("rstCurrentBurn was not supplied and is not being generated. ",
                               "The module will try to recover it from the inputs folder. ",
                               " If this fails, the forestry activities will be considered ",
                               "without a fire regime!")))
  }
  
  if (!suppliedElsewhere(object = "disturbanceDT", sim = sim)) {
    sim$disturbanceDT <- prepInputs(url = extractURL("disturbanceDT"),
                                    targetFile = "disturbanceDT.csv",
                                    destinationPath = dPath,
                                    fun = "data.table::fread",
                                    header = TRUE, 
                                    userTags = "disturbanceDT")
    
    warning(paste0("disturbanceDT was not supplied. Defaulting to an example from ",
                   " Northwest Territories"), immediate. = TRUE)
  }
  
  if (!suppliedElsewhere(object = "DisturbanceRate", sim = sim)){
    sim$DisturbanceRate <- NULL
    warning(paste0("The table DisturbanceRate was not supplied. The module will use ECCC data to",
                   "calculate disturbance rates.",
                   "If other rates are desired, please provide DisturbanceRate"), 
            immediate. = TRUE)
  }
  
  return(invisible(sim))
}

