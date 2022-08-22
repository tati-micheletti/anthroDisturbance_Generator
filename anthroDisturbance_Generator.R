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
  version = list(anthroDisturbance_Generator = "0.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "anthroDisturbance_Generator.Rmd"), ## same file
  reqdPkgs = list("SpaDES.core (>=1.0.10)", "ggplot2", 
                  "data.table", "PredictiveEcology/reproducible@development",
                  "raster", "terra", "crayon", "msm", "rgdal", "sf", 
                  "fasterize", "stars", "nngeo", "tictoc"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
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
    defineParameter("saveInitialDisturbances", "logical", TRUE, NA, NA,
                    paste0("Should the disturbance rasters be saved at each step? These are saved ",
                           "to Paths[['outputPath']] as a RasterLayer, with disturbanceLayer as prefix",
                           "the name of the industry and the year as suffix.",
                           "If TRUE, it saves the initial conditions (IC)")),
    defineParameter("saveCurrentDisturbances", "logical", TRUE, NA, NA,
                    paste0("Should the disturbance rasters be saved at each step? These are saved ",
                           "to Paths[['outputPath']] as a RasterLayer, with disturbanceLayer as prefix",
                           "the name of the industry and the year as suffix.",
                           "If TRUE, it saves at the end of each step.")),
    defineParameter("useECCCData", "logical", TRUE, NA, NA,
                    paste0("If the rate of change is not provided, the module can either try to ",
                    "extract it from data if useECCCData = TRUE (i.e., ECCC human footprint for 2010 and 2015), or ",
                    "simply apply a rate of 0.2% (of the total area) increase per year")),
    defineParameter("checkChangeInDisturbance", "logical", FALSE, NA, NA,
                    paste0("Prints the change in ECCC human footprint (2010 to 2015) at 30m resolution",
                           " over the shapefile provided")),
    defineParameter("checkDisturbance2015", "logical", FALSE, NA, NA,
                    paste0("Prints the total % of ECCC human footprint (2010 to 2015) at 500m resolution",
                           " over the shapefile provided")),
    defineParameter("overwriteDisturbanceLayers2015", "logical", FALSE, NA, NA,
                    paste0("Should the disturbance layer from 2015 be overwritten?")),
    defineParameter("overwriteDisturbanceLayers2010", "logical", FALSE, NA, NA,
                    paste0("Should the disturbance layer from 2010 be overwritten?")),
    defineParameter("growthStepEnlargingPolys", "numeric", 0.0075, NA, NA,
                    paste0("Growth step used for iteratively achieving the total area growth of ",
                           "new disturbances type Enlarging for polygons. If the iterations take too",
                           " long, one should increase this number. If the summarized value is too",
                           " far from 0, one should decrease this number.")),
    defineParameter("growthStepEnlargingLines", "numeric", 0.1, NA, NA,
                    paste0("Growth step used for iteratively achieving the total area growth of ",
                           "new disturbances type Enlarging for lines. If the iterations take too",
                           " long, one should increase this number. If the summarized value is too",
                           " far from 0, one should decrease this number."))
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
                               "disturbances per year of type Enlarging and Generating. For ",
                               "disturbances type Connecting, disturbanceRate is NA. ",
                               "If not specified when needed, the module will try to derive ",
                               "it from data. If this fails, it will fall on a yearly average ",
                               " of 0.2% of the current disturbance (except for windTurbine, ",
                               "which has a default value of 1 per 10 years considering the reduced ",
                               "of potential in the region.",
                               
                               "disturbanceSize --> if there is a specific size the disturbance in m2 ",
                               "type Generating should have, it is specified here. If not specified, ",
                               " the module will try to derive it from data. For disturbances ",
                               "type Enlarging anb Connecting, disturbanceSize is NA",
                               
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
                               "potentialTransmission, potentialRoads incl. ",
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
      
      if (P(sim)$saveInitialDisturbances){
        message(crayon::yellow(paste0("The parameter saveInitialDisturbances is TRUE.",
                                      " Saving initial disturbance layers")))
        saveDisturbances(disturbanceList = sim$disturbanceList,
                         currentTime = "IC", overwrite = TRUE)
      }

      mod$.whichToRun <- whichDisturbancesToGenerate(startTime = start(sim),
                                                     currentTime = time(sim),
                                                     endTime = end(sim),
                                                     disturbanceParameters = sim$disturbanceParameters)
      sim$currentDisturbanceLayer <- list()
      # schedule future event(s)
        sim <- scheduleEvent(sim, time(sim), "anthroDisturbance_Generator", "calculatingSize")
        sim <- scheduleEvent(sim, time(sim), "anthroDisturbance_Generator", "calculatingRate")
        sim <- scheduleEvent(sim, time(sim), "anthroDisturbance_Generator", "generatingDisturbances")
        sim <- scheduleEvent(sim, time(sim), "anthroDisturbance_Generator", "updatingDisturbanceList")
    },
    calculatingSize = {
      # Check for needed sizes. If all provided, skip event
      mod$.whichNeedSize <- which(sim$disturbanceParameters[, disturbanceType == "Generating" &
                                                              is.na(disturbanceSize)])
      if (all(length(mod$.whichToRun) != 0,
              length(mod$.whichNeedSize) != 0))
        sim$disturbanceParameters <- calculateSize(disturbanceParameters = sim$disturbanceParameters,
                                                 whichToUpdate = mod$.whichNeedSize,
                                                 disturbanceList = sim$disturbanceList)
    },
    calculatingRate = {
      # Check for needed rates. If all provided, skip event
      mod$.whichNeedRates <- which(sim$disturbanceParameters[, disturbanceType %in% c("Generating", "Enlarging") &
                                                               is.na(disturbanceRate)])
      if (all(length(mod$.whichToRun) != 0,
              length(mod$.whichNeedRates) != 0))
        sim$disturbanceParameters <- calculateRate(disturbanceParameters = sim$disturbanceParameters,
                                                 whichToUpdate = mod$.whichNeedRates,
                                                 studyArea = sim$studyArea,
                                                 disturbanceList = sim$disturbanceList,
                                                 RTM = sim$rasterToMatch,
                                                 useECCCData = P(sim)$useECCCData,
                                                 disturbanceDT = sim$disturbanceDT,
                                                 checkChangeInDisturbance = P(sim)$checkChangeInDisturbance,
                                                 checkDisturbance2015 = P(sim)$checkDisturbance2015,
                                                 overwriteDisturbanceLayers2010 = P(sim)$overwriteDisturbanceLayers2010,
                                                 overwriteDisturbanceLayers2015 = P(sim)$overwriteDisturbanceLayers2015)
      },
    generatingDisturbances = {
      # Check if the time(sim) is within the interval to run the disturbances
      mod$.whichToRun <- whichDisturbancesToGenerate(startTime = start(sim),
                                                     currentTime = time(sim),
                                                     endTime = end(sim),
                                                     disturbanceParameters = sim$disturbanceParameters)
      if (length(mod$.whichToRun) != 0){ # If anything is scheduled
        forestryScheduled <- "forestry" %in% sim$disturbanceParameters[mod$.whichToRun, dataName]
        if (forestryScheduled) { # forestry scheduled?
          if (is.null(sim$rstCurrentBurn)){
            # Try to get rstCurrentBurn if no module is producing it!
            mod$rstCurrentBurn <- createModObject(data = "rstCurrentBurn", 
                                                  sim = sim, 
                                                  pathInput = inputPath(sim), 
                                                  currentTime = time(sim),
                                                  returnNULL = TRUE)
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
        mod$updatedLayers <- generateDisturbances(disturbanceParameters = sim$disturbanceParameters,
                                                  # disturbanceParameters = sim$disturbanceParameters[disturbanceOrigin != "cutblocks", ],
                                                  disturbanceList = sim$disturbanceList,
                                                  # disturbanceList = sim$disturbanceList[names(sim$disturbanceList) != "forestry"],
                                                  currentTime = time(sim),
                                                  studyArea = sim$studyArea,
                                                  rasterToMatch = sim$rasterToMatch,
                                                  fires = mod$rstCurrentBurn,
                                                  currentDisturbanceLayer = currDis,
                                                  growthStepEnlargingPolys = P(sim)$growthStepEnlargingPolys,
                                                  growthStepEnlargingLines = P(sim)$growthStepEnlargingLines)
        
        # mod$updatedLayers <- generateDisturbances_forestry(disturbanceParameters = sim$disturbanceParameters[disturbanceOrigin == "cutblocks", ],
        #                                           disturbanceList = sim$disturbanceList[names(sim$disturbanceList) == "forestry"],
        #                                           fires = mod$rstCurrentBurn,
        #                                           currentTime = time(sim),
        #                                           rasterToMatch = sim$rasterToMatch,
        #                                           currentDisturbanceLayer = mod$updatedLayers)
        
        sim$currentDisturbanceLayer[paste0("Year", time(sim))] <- mod$updatedLayers$currentDisturbanceLayer
      }
      sim <- scheduleEvent(sim, time(sim) + 1, "anthroDisturbance_Generator", "generatingDisturbances")
    },
    updatingDisturbanceList = {
      if (length(mod$.whichToRun) != 0){
      sim$disturbanceList <- replaceList(disturbanceList = sim$disturbanceList, 
                                         updatedLayers = mod$updatedLayers$individuaLayers)
      if (P(sim)$saveCurrentDisturbances){
        saveDisturbances(disturbanceList = sim$disturbanceList,
                         currentTime = time(sim))
        }
      }

      sim <- scheduleEvent(sim, time(sim) + 1, "anthroDisturbance_Generator", "updatingDisturbanceList")

    },
    warning(paste("Undefined event type: \'", current(sim)[1, "eventType", with = FALSE],
                  "\' in module \'", current(sim)[1, "moduleName", with = FALSE], "\'", sep = ""))
  )
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  # Any code written here will be run during the simInit for the purpose of creating
  # any objects required by this module and identified in the inputObjects element of defineModule.
  # This is useful if there is something required before simulation to produce the module
  # object dependencies, including such things as downloading default datasets, e.g.,
  # downloadData("LCC2005", modulePath(sim)).
  # Nothing should be created here that does not create a named object in inputObjects.
  # Any other initiation procedures should be put in "init" eventType of the doEvent function.
  # Note: the module developer can check if an object is 'suppliedElsewhere' to
  # selectively skip unnecessary steps because the user has provided those inputObjects in the
  # simInit call, or another module will supply or has supplied it. e.g.,
  # if (!suppliedElsewhere('defaultColor', sim)) {
  #   sim$map <- Cache(prepInputs, extractURL('map')) # download, extract, load file from url in sourceURL
  # }

  #cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
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
                   " be used for module testing purposes! Please run the module(s) ",
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

  return(invisible(sim))
}

