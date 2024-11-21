message("PRE-PROCESSING - HBV CHRONIC MOUSE PART")
message("***************************************")


require("CytoPipeline")         ## from Bioconductor
require("CytoPipelineUtils")    ## from https://github.com/UCLouvain-CBIO/CytoPipelineUtils

## preliminaries: path to raw data files, json input, compensation input, 
## and outputs
rawDataDir <- "./data/HBV_chronic_mouse/rawData"
jsonDir <- "./json"
compensationDir <- "./data/HBV_chronic_mouse/compensation"
resultsDir <- "./preprocessing"

## Biocparallel set up
## create appropriate BiocParallel back-end

rmCache <- TRUE
useBiocParallel <- TRUE

logDir <- "./BPLog/prep"
dir.create(logDir, showWarnings = FALSE, recursive = TRUE)

bp <- BiocParallel::SnowParam(progressbar = TRUE,
                              log = TRUE,
                              logdir = logDir,
                              stop.on.error = FALSE)

## build CytoPipeline object

expName <- "HBV_chronic_mouse"
              
jsonFile <- file.path(jsonDir,
                      paste0(expName, "_CytoPipeline_Preprocessing.json"))

(sampleFiles <- list.files(path = rawDataDir,
                           pattern = ".fcs",
                           full.names = TRUE))

## construct phenoData data.frame
day <- factor(substr(basename(sampleFiles), 1, 3))
well <- substr(basename(sampleFiles), 5, 7)
name <- sapply(basename(sampleFiles), 
                  FUN = function(nm){
                      substr(nm, 1, nchar(nm)-4)
                  } )
ffName <- paste0(name, ".fcs")
group <- substr(name, 6, 7)
group <- sapply(group,
                FUN = function(grstr) {
                    if (grstr == "01") 1 
                    else if (grstr == "03") 2
                    else if (grstr == "05") 3
                    else if (grstr == "07") 4
                    else  5
                })
group <- factor(group)
sex <- factor(c("male",                                         #D91_G well
                "female", "female", "female", "female", "male", #D91_A well
                "female", "male", "male", "male", "male",       #D91_B well
                "male", "male", "male", "male", "male",         #D91_C well
                "male", "male", "male", "male",                 #D91_D well
                "male", "male", "male",                         #D91_E well
                "male", "male", "male", "male",                 #D91_F well
                "female", "female",                             #D93_G well
                "female", "female", "female", "female", "male", #D93_A well
                "female", "female", "female", "female", "male", #D93_B well
                "female", "female", "female", "female",         #D93_C well
                "female", "female", "female", "female",         #D93_D well
                "female", "female", "female", "female",         #D93_E well
                "female", "female", "female", "female"          #D93_F well
))
pData <- data.frame(
    ffName = ffName,
    name = name,
    day = day,
    group = group,
    well = well,
    sex = sex,
    row.names = 1)

pData$compensation <- ifelse(
    pData$day == "D91",
    file.path(compensationDir, "Compensations Liver D91.csv"),
    file.path(compensationDir, "Compensations Liver D93.csv"))

pipL <- CytoPipeline(jsonFile,
                     experimentName = expName,
                     sampleFiles = sampleFiles,
                     pData = pData)

## run CytoPipeline objects
execute(
    pipL, 
    path = resultsDir,
    rmCache = rmCache,
    useBiocParallel = useBiocParallel,
    BPPARAM = bp,
    BPOPTIONS = BiocParallel::bpoptions(
        package = c("flowCore",
                    "CytoPipelineUtils")),
    saveScaleTransforms = FALSE,
    saveLastStepFF = FALSE)

message("ALL PRE-PROCESSING DONE !")
                
