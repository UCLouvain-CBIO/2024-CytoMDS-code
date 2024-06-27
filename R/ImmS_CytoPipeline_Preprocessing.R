message("PRE-PROCESSING - IMMUNOSENESCENCE HUMAN PBMC PART")
message("*************************************************")


require("CytoPipeline")         ## from Bioconductor
require("CytoPipelineUtils")    ## from https://github.com/UCLouvain-CBIO/CytoPipelineUtils

# output dir
resultsDir <- "./preprocessing"

rmCache <- TRUE
useBiocParallel <- TRUE

jsonFile <- "./json/ImmS_CytoPipeline_Preprocessing.json"
expName <- "ImmunoSenescence_human_PBMC"

rawDataDir <- "./data/ImmunoSenescence_human_PBMC/rawData"

sampleFiles <- list.files(
    rawDataDir,
    pattern = "*.fcs",
    full.names = TRUE
)

ffName <- substr(basename(sampleFiles), start = 1, stop = 3)
panel <- substr(ffName, start = 1, stop = 1)
panel <- ifelse(panel == "F", "former", "later")
panel <- factor(panel)
patientId <- as.integer(substr(ffName, start = 2, stop = 3))
patientId <- factor(patientId)
group <- sapply(basename(sampleFiles), 
                  FUN = function(nm){
                      substr(nm, 5, nchar(nm)-4)
                  } )
group <- factor(group)

phenoData <- data.frame(ffName = ffName,
                        group = group,
                        patientId = patientId,
                        panel = panel,
                        file = basename(sampleFiles))
row.names(phenoData) <- phenoData$file

phenoData$compensation <- ifelse(
    phenoData$panel == "former",
    "data/ImmunoSenescence_human_PBMC/compensation/formerPanel_compensation.csv",
    "data/ImmunoSenescence_human_PBMC/compensation/laterPanel_compensation.csv")

sel <- 1:20
#sel <- 1
sampleFiles <- sampleFiles[sel]
phenoData <- phenoData[sel,]


# creation on CytoPipeline object,
# using json file as input
pipL <- CytoPipeline(
    jsonFile,
    experimentName = expName,
    sampleFiles = sampleFiles,
    pData = phenoData)

# execute PeacoQC pipeline

logDir <- file.path(".","BPLog" ,"prep")
# create log dir if needed
dir.create(
    logDir,
    recursive = TRUE,
    showWarnings = FALSE
)

bp <- BiocParallel::SnowParam(progressbar = TRUE,
                              log = TRUE,
                              logdir = logDir,
                              stop.on.error = FALSE)
execute(pipL, path = resultsDir, 
        rmCache = rmCache, useBiocParallel = useBiocParallel,
        BPPARAM = bp,
        BPOPTIONS = BiocParallel::bpoptions(package = c("flowCore", 
                                                        "CytoPipelineUtils")),
        saveScaleTransforms = FALSE,
        saveLastStepFF = TRUE)


message("ALL PRE-PROCESSING DONE !")
