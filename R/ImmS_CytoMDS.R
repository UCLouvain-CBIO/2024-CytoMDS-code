message("FIGURES - IMMUNOSENESCENCE HUMAN PBMC PART")
message("******************************************")


require(CytoPipeline)
require(CytoPipelineGUI)
require(CytoMDS)
require(ggplot2)
require(patchwork)


message("Getting Results...")

# directory from which we read the preprocessing results
prepDir <- "./preprocessing"


expName <- "ImmunoSenescence_human_PBMC"

prepstep <- "compensate"
pipelineObj <- paste0(prepstep, "_obj")

## marker selection: either pre-processing channels or all 
#channelTypes <- c("prep", "all")
channelTypes <- c("prep")
selChannelList <- setNames(vector("list", 2), channelTypes)
selChannelList[["prep"]] <- c("FSC-A", "SSC-A", "FSC-H", "LD")

# recreate CytoPipeline objects from cache
pipL <- CytoPipeline::buildCytoPipelineFromCache(
    experimentName = expName, path = prepDir
)

inputFF <- CytoPipeline::getCytoPipelineFlowFrame(
    pipL, 
    path = prepDir,
    whichQueue = "pre-processing",
    sampleFile = 1,
    objectName = pipelineObj
)
selChannelList[["all"]] <- flowCore::colnames(inputFF)[c(2,5,8:33)]

# reading target flowSet
sampleFiles <- sampleFiles(pipL)

allPhenoData <- CytoPipeline::pData(pipL)

allPhenoData["data_acquisition"] <- allPhenoData["panel"]

N <- length(sampleFiles)

tgtSeq <- 1:N
ffList <- list()
for(i in seq_along(tgtSeq)) {
    message("Reading file #", i, "/", N, "...")
    ffList[[i]] <- getCytoPipelineFlowFrame(
        pipL,
        whichQueue = "pre-processing",
        sampleFile = tgtSeq[i],
        objectName = pipelineObj,
        path = prepDir)
}

names(ffList) <- rownames(allPhenoData)

fsAll <- as(ffList,"flowSet")

flowCore::pData(fsAll) <- allPhenoData

transList <- CytoPipeline::getCytoPipelineScaleTransform(
    pipL,
    whichQueue = "scale transform",
    objectName = "scale_transform_estimate_obj",
    path = prepDir)

message("applying scale transformations to all flow frames...")
fsAll <- CytoPipeline::applyScaleTransforms(
    fsAll,
    transList = transList,
    verbose = TRUE
)

# definition of possible external variables for biplots
statFUNs = list("median" = stats::median, "std-dev" = stats::sd)


statFUNs[["Q10"]] <- function(x, na.rm) {
    stats::quantile(x, probs = 0.1)
}

statFUNs[["Q20"]] <- function(x, na.rm) {
    stats::quantile(x, probs = 0.2)
}

statFUNs[["Q80"]] <- function(x, na.rm) {
    stats::quantile(x, probs = 0.8)
}

statFUNs[["Q90"]] <- function(x, na.rm) {
    stats::quantile(x, probs = 0.9)
}

nCompute <- length(channelTypes)

verbose <- FALSE

resultArray <- vapply(
    channelTypes,
    USE.NAMES = TRUE,
    FUN.VALUE = vector("list", 3),
    FUN = function(chType) {
        message("Launching calculations for channel type : ", chType)
        message("**************************************************")
        
        # computing channel stats
        
        message("Computing channel stats for all flow frames...")
        
        chStats <- channelSummaryStats(
            fsAll,
            channels = selChannelList[[chType]],
            statFUNs = statFUNs,
            verbose = verbose)
        
        message("Computing pairwise distances...")
        startTime <- Sys.time()
        pwDist <- 
            CytoMDS::pairwiseEMDDist(
                x = fsAll,
                channels = selChannelList[[chType]],
                verbose = verbose
            )
        endTime <- Sys.time()
        if (verbose) {
            # prints recorded time
            print(endTime - startTime)    
        }
        
        message("Computing MDS...")
        mds <- computeMetricMDS(pwDist, seed = 0)
        
        list(chStats = chStats,
             pwDist = pwDist,
             mds = mds)
    }
)

message("ALL DONE ! :-)")

## plots for preprocessing channels only

selChType <- "prep"

phenoData <- allPhenoData
chStats <- resultArray["chStats", selChType][[1]]
pwDist <- resultArray["pwDist", selChType][[1]]
mds <- resultArray["mds", selChType][[1]]

distVec <- pwDist[upper.tri(pwDist)]

# display histogram of all pairwise distances
# hist(distVec, 
#      main = paste0(
#          "Pairwise distances (step = ",
#          prepstep,
#          "; channels = ",
#          selChType,
#          ")"))

#nDim(mds)
#CytoMDS::RSq(mds)

# Figure S3

message("Generating Figure S3...")
pSh <- ggplotSampleMDSShepard(mds) + theme_bw()
ggplotResults(pSh, name = "FigS3")
message("Done!")

pDataForShape <- "data_acquisition"
pDataForLabel <- "ffName"
pDataForColour <- "group"
pDataForAdditionalLabelling <- c("file", "patientId")

# Figure 3
message("Generating Figure 3...")
p1 <- ggplotSampleMDS(
    mdsObj = mds,
    pData = phenoData,
    projectionAxes = c(1,2),
    pDataForColour = pDataForColour,
    pDataForShape = pDataForShape,
    pDataForLabel = pDataForLabel,
    pDataForAdditionalLabelling = 
        pDataForAdditionalLabelling,
    title = "MDS - projection axes 1 and 2") + 
    theme(legend.position="none") +
    theme_bw() + 
    scale_shape_manual(values = c(17, 15)) +
    scale_color_manual(values = c("red", "blue"))

p2 <- ggplotSampleMDS(
    mdsObj = mds,
    pData = phenoData,
    projectionAxes = c(2,3),
    pDataForColour = pDataForColour,
    pDataForShape = pDataForShape,
    pDataForLabel = pDataForLabel,
    pDataForAdditionalLabelling = 
        pDataForAdditionalLabelling,
    title = "MDS - projection axes 2 and 3") +
    theme_bw() + 
    scale_shape_manual(values = c(17, 15)) + 
    scale_color_manual(values = c("red", "blue"))

pRes <- p1 + p2
ggplotResults(pRes, name = "Fig3")
message("Done!")

# Figure S4

message("Generating Figure S4 - three parts...")
p01 <- ggplotSampleMDS(
    mdsObj = mds,
    pData = phenoData,
    projectionAxes = c(1,2),
    pDataForColour = pDataForColour,
    pDataForShape = pDataForShape,
    pDataForLabel = pDataForLabel,
    pDataForAdditionalLabelling = 
        pDataForAdditionalLabelling,
    title = "MDS - projection axes 1 and 2") +
    theme_bw() + 
    scale_shape_manual(values = c(17, 15)) +
    scale_color_manual(values = c("red", "blue"))

ggplotResults(p01, name = "FigS4_main_part")


bp11 <- ggplotSampleMDS(
    mdsObj = mds,
    pData = phenoData,
    biplot = TRUE,
    extVariables = chStats[["median"]],
    projectionAxes = c(1,2),
    pDataForColour = pDataForColour,
    pDataForShape = pDataForShape,
    pDataForLabel = pDataForLabel,
    pDataForAdditionalLabelling = 
        pDataForAdditionalLabelling,
    displayPointLabels = FALSE,
    title = "Bi-plot with medians") + 
    theme(legend.position="none") +
    theme_bw() + 
    scale_shape_manual(values = c(17, 15)) +
    scale_color_manual(values = c("red", "blue"))


bp12 <- ggplotSampleMDS(
    mdsObj = mds,
    pData = phenoData,
    biplot = TRUE,
    extVariables = chStats[["std-dev"]],
    projectionAxes = c(1,2),
    pDataForColour = pDataForColour,
    pDataForShape = pDataForShape,
    pDataForLabel = pDataForLabel,
    pDataForAdditionalLabelling = 
        pDataForAdditionalLabelling,
    displayPointLabel = FALSE,
    title = "Bi-plot with standard deviations") + 
    theme(legend.position="none") + 
    theme_bw() + 
    scale_shape_manual(values = c(17, 15)) +
    scale_color_manual(values = c("red", "blue"))

ggplotResults(bp11 / bp12, name = "FigS4_biplots")


file1 <- "F05_Young.fcs"
file2 <- "L05_Young.fcs"

pF11 <- CytoPipelineGUI::plotSelectedFlowFrame(
    experimentName = expName,
    whichQueue = "pre-processing",
    sampleFile = file1,
    flowFrameName = "compensate_obj",
    path = prepDir,
    xChannelLabel = "FSC-A : NA",
    yChannelLabel = "Comp-BV421-A : LD",
    useAllCells = FALSE,
    nDisplayCells = 10000,
    useFixedLinearRange = TRUE,
    linearRange = c(-300, 262144),
    transfoListName = "scale_transform_estimate_obj") + 
    labs(subtitle = "")

pF12 <- CytoPipelineGUI::plotSelectedFlowFrame(
    experimentName = expName,
    whichQueue = "pre-processing",
    sampleFile = file2,
    flowFrameName = "compensate_obj",
    path = prepDir,
    xChannelLabel = "FSC-A : NA",
    yChannelLabel = "Comp-BV421-A : LD",
    useAllCells = FALSE,
    nDisplayCells = 10000,
    useFixedLinearRange = TRUE,
    linearRange = c(-300, 262144),
    transfoListName = "scale_transform_estimate_obj") +
    labs(subtitle = "")

ggplotResults(pF11 + pF12, name = "FigS4_indiv_samples")
message("Done!")

# Figure S5
message("Generating Figure S5 - three parts...")
p02 <- ggplotSampleMDS(
    mdsObj = mds,
    pData = phenoData,
    projectionAxes = c(2,3),
    pDataForColour = pDataForColour,
    pDataForShape = pDataForShape,
    pDataForLabel = pDataForLabel,
    pDataForAdditionalLabelling = 
        pDataForAdditionalLabelling,
    title = "MDS - projection axes 2 and 3") +
    theme_bw() + 
    scale_shape_manual(values = c(17, 15)) +
    scale_color_manual(values = c("red", "blue"))

ggplotResults(p02, name = "FigS5_main_part")


bp21 <- ggplotSampleMDS(
    mdsObj = mds,
    pData = phenoData,
    biplot = TRUE,
    extVariables = chStats[["median"]],
    projectionAxes = c(2,3),
    pDataForColour = pDataForColour,
    pDataForShape = pDataForShape,
    pDataForLabel = pDataForLabel,
    pDataForAdditionalLabelling = 
        pDataForAdditionalLabelling,
    repelArrowLabels = TRUE,
    displayPointLabels = FALSE,
    title = "Bi-plot with medians") + 
    theme(legend.position="none") +
    theme_bw() + 
    scale_shape_manual(values = c(17, 15)) +
    scale_color_manual(values = c("red", "blue"))

bp22 <- ggplotSampleMDS(
    mdsObj = mds,
    pData = phenoData,
    biplot = TRUE,
    extVariables = chStats[["Q10"]],
    projectionAxes = c(2,3),
    pDataForColour = pDataForColour,
    pDataForShape = pDataForShape,
    pDataForLabel = pDataForLabel,
    pDataForAdditionalLabelling = 
        pDataForAdditionalLabelling,
    repelArrowLabels = TRUE,
    displayPointLabel = FALSE,
    title = "Bi-plot with 10th quantile") + 
    theme(legend.position="none") + 
    theme_bw() + 
    scale_shape_manual(values = c(17, 15)) +
    scale_color_manual(values = c("red", "blue"))

ggplotResults(bp21 / bp22, name = "FigS5_biplots")


#CytoPipelineGUI::CytoPipelineCheckApp(dir = prepDir)

file1 <- "L18_Old.fcs"
file2 <- "F12_Old.fcs"
file3 <- "F01_Young.fcs"
file4 <- "L01_Young.fcs"


pF21 <- CytoPipelineGUI::plotSelectedFlowFrame(
    experimentName = expName,
    whichQueue = "pre-processing",
    sampleFile = file1,
    flowFrameName = "compensate_obj",
    path = prepDir,
    xChannelLabel = "FSC-A : NA",
    yChannelLabel = " ",
    useAllCells = TRUE,
    nDisplayCells = NULL,
    useFixedLinearRange = TRUE,
    linearRange = c(-300, 262144),
    transfoListName = " ") + labs(subtitle = "") + 
    scale_y_continuous(limits = c(0, 3e-5))

pF22 <- CytoPipelineGUI::plotSelectedFlowFrame(
    experimentName = expName,
    whichQueue = "pre-processing",
    sampleFile = file2,
    flowFrameName = "compensate_obj",
    path = prepDir,
    xChannelLabel = "FSC-A : NA",
    yChannelLabel = " ",
    useAllCells = TRUE,
    nDisplayCells = NULL,
    useFixedLinearRange = TRUE,
    linearRange = c(-300, 262144),
    transfoListName = " ") + labs(subtitle = "") + 
    scale_y_continuous(limits = c(0, 3e-5))

pF23 <- CytoPipelineGUI::plotSelectedFlowFrame(
    experimentName = expName,
    whichQueue = "pre-processing",
    sampleFile = file3,
    flowFrameName = "compensate_obj",
    path = prepDir,
    xChannelLabel = "FSC-A : NA",
    yChannelLabel = " ",
    useAllCells = TRUE,
    nDisplayCells = NULL,
    useFixedLinearRange = TRUE,
    linearRange = c(-300, 262144),
    transfoListName = " ") + labs(subtitle = "") + 
    scale_y_continuous(limits = c(0, 3e-5))

pF24 <- CytoPipelineGUI::plotSelectedFlowFrame(
    experimentName = expName,
    whichQueue = "pre-processing",
    sampleFile = file4,
    flowFrameName = "compensate_obj",
    path = prepDir,
    xChannelLabel = "FSC-A : NA",
    yChannelLabel = " ",
    useAllCells = TRUE,
    nDisplayCells = NULL,
    useFixedLinearRange = TRUE,
    linearRange = c(-300, 262144),
    transfoListName = " ") + labs(subtitle = "") + 
    scale_y_continuous(limits = c(0, 3e-5))

pRes <- 
    (pF21 + pF23 ) / (pF22 + pF24  )

ggplotResults(pRes, name = "FigS5_indiv_samples")
message("Done!")

## Figure S8 - projection using channel medians
message("Figure S8 - calculating MDS with medians...")
pwDistMedians <- as.matrix(dist(chStats[["median"]]))

# display histogram of all pairwise distances
# hist(pwDistMedians,
#      main = paste0(
#          "Pairwise distances using medians (channels = ",
#          selChType,
#          ")"
#      ))

mdsMedians <- computeMetricMDS(pwDistMedians, nDim = 4)
#nDim(mdsMedians)
#RSq(mdsMedians)
message("Done!")

message("Generating Figure S10...")
p1Medians <- ggplotSampleMDS(
    mdsObj = mdsMedians,
    pData = phenoData,
    projectionAxes = c(1,2),
    pDataForColour = pDataForColour,
    pDataForShape = pDataForShape,
    pDataForLabel = pDataForLabel,
    pDataForAdditionalLabelling = 
        pDataForAdditionalLabelling,
    title = "MDS - projection axes 1 and 2") + 
    scale_shape_manual(values = c(17, 15)) +
    scale_color_manual(values = c("red", "blue")) 

p2Medians <- ggplotSampleMDS(
    mdsObj = mdsMedians,
    pData = phenoData,
    projectionAxes = c(2,3),
    pDataForColour = pDataForColour,
    pDataForShape = pDataForShape,
    pDataForLabel = pDataForLabel,
    pDataForAdditionalLabelling = 
        pDataForAdditionalLabelling,
    displayPointLabels = FALSE,
    flipYAxis = TRUE,
    max.overlaps = 100,
    title = "MDS - projection axes 2 and 3") +
    scale_shape_manual(values = c(17, 15)) + 
    scale_color_manual(values = c("red", "blue"))

pp1 <- p1 + labs(
    title = "CytoMDS - coordinates 1 and 2",
    subtitle = ""
) + theme(legend.position = "none")

pp2 <- p1Medians + labs(
    title = "MDS with medians - coordinates 1 and 2",
    subtitle = ""
) + theme_bw() 

pp3 <- p2 + labs(
    title = "CytoMDS - coordinates 2 and 3",
    subtitle = "") + 
    theme(legend.position = "none")

pp4 <- p2Medians + labs(
    title = "MDS with medians - coordinates 2 and 3",
    subtitle = ""
) + theme_bw()

pRes <- 
    (pp1 + pp2) / (pp3 + pp4)

ggplotResults(pRes, name = "FigS8")
message("Done!")



