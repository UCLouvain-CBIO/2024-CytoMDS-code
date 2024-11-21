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

names(ffList) <- rownames(allPhenoData)[tgtSeq]

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

#display histogram of all pairwise distances
hist(distVec,
     main = paste0(
         "Pairwise distances (step = ",
         prepstep,
         "; channels = ",
         selChType,
         ")"))

#nDim(mds)
#CytoMDS::RSq(mds)

# saveRDS(object = mds,
#         file = "./rds/ImmS_mdsObj.rds")
# 
# saveRDS(object = phenoData,
#         file = "./rds/ImmS_phenoData.rds")
# 
# saveRDS(object = chStats,
#         file = "./rds/ImmS_stats.rds")

pointSizeShepard <- 1.0
pointSizeMDS <- 2.5

# Figure S3

message("Generating Figure S3...")

pSh <- ggplotSampleMDSShepard(mds, 
                              lineWidth = 1.0, 
                              pointSize = pointSizeShepard,
                              title = "") + 
    scale_x_continuous(limits = c(0,2)) + 
    theme_bw() + theme(plot.title = element_text(size = 20),
                       plot.subtitle = element_text(size = 15),
                       axis.title = element_text(size = 15),
                       axis.text = element_text(size = 15),
                       legend.title = element_text(size = 15),
                       legend.text = element_text(size = 12))

ggplotResults(pSh, name = "FigS3_ImmS_Shepard", width = 480, height = 480)
message("Done!")

# Figure 3
message("Generating Figure 3...")

pDataForShape <- "data_acquisition"
pDataForLabel <- "ffName"
pDataForColour <- "group"
pDataForAdditionalLabelling <- c("file", "patientId")

p1 <- ggplotSampleMDS(
    mdsObj = mds,
    pData = phenoData,
    projectionAxes = c(1,2),
    pDataForColour = pDataForColour,
    pDataForShape = pDataForShape,
    pDataForLabel = pDataForLabel,
    pDataForAdditionalLabelling = 
        pDataForAdditionalLabelling,
    pointSize = pointSizeMDS,
    title = "MDS - projection axes 1 and 2") + 
    scale_shape_manual(values = c(17, 15)) +
    scale_color_manual(values = c("red", "blue")) +
    theme_bw()

p1Th <-  p1 + 
    theme(legend.position="none",
          plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15))

p2 <- ggplotSampleMDS(
    mdsObj = mds,
    pData = phenoData,
    projectionAxes = c(2,3),
    pDataForColour = pDataForColour,
    pDataForShape = pDataForShape,
    pDataForLabel = pDataForLabel,
    pDataForAdditionalLabelling = 
        pDataForAdditionalLabelling,
    pointSize = pointSizeMDS,
    title = "MDS - projection axes 2 and 3") +
    scale_shape_manual(values = c(17, 15)) + 
    scale_color_manual(values = c("red", "blue")) +
    theme_bw() 
    

p2Th <- p2 + 
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12))
    

pRes <- p1Th + p2Th
ggplotResults(pRes, name = "Fig3", width = 960, height = 480)
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
    pointSize = pointSizeMDS,
    title = "MDS - projection axes 1 and 2") +
    theme_bw() + 
    scale_shape_manual(values = c(17, 15)) +
    scale_color_manual(values = c("red", "blue")) + 
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12))

ggplotResults(p01, name = "FigS4_main_part", width = 540, height = 480)


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
    pointSize = pointSizeMDS,
    displayPointLabels = FALSE,
    title = "Bi-plot with medians") + 
    theme(legend.position="none") +
    theme_bw() + 
    scale_shape_manual(values = c(17, 15)) +
    scale_color_manual(values = c("red", "blue")) + 
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12))


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
    pointSize = pointSizeMDS,
    displayPointLabel = FALSE,
    title = "Bi-plot with standard deviations") + 
    theme(legend.position="none") + 
    theme_bw() + 
    scale_shape_manual(values = c(17, 15)) +
    scale_color_manual(values = c("red", "blue")) + 
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12))

ggplotResults(bp11 / bp12, name = "FigS4_biplots", width = 570, height = 960)


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
    labs(subtitle = "") + 
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 15))

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
    labs(subtitle = "") + 
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 15))

ggplotResults(pF11 + pF12, name = "FigS4_indiv_samples", width = 960, height = 480)
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
    pointSize = pointSizeMDS,
    title = "MDS - projection axes 2 and 3") +
    theme_bw() + 
    scale_shape_manual(values = c(17, 15)) +
    scale_color_manual(values = c("red", "blue")) + 
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12))

ggplotResults(p02, name = "FigS5_main_part", width = 540, height = 480)


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
    pointSize = pointSizeMDS,
    repelArrowLabels = TRUE,
    displayPointLabels = FALSE,
    title = "Bi-plot with medians") + 
    theme(legend.position="none") +
    theme_bw() + 
    scale_shape_manual(values = c(17, 15)) +
    scale_color_manual(values = c("red", "blue")) + 
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12))

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
    pointSize = pointSizeMDS,
    repelArrowLabels = TRUE,
    displayPointLabel = FALSE,
    title = "Bi-plot with 10th quantile") + 
    theme(legend.position="none") + 
    theme_bw() + 
    scale_shape_manual(values = c(17, 15)) +
    scale_color_manual(values = c("red", "blue")) + 
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12))

ggplotResults(bp21 / bp22, name = "FigS5_biplots", width = 570, height = 960)


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
    scale_y_continuous(limits = c(0, 3e-5)) + 
    labs(subtitle = "") + 
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 15))

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
    scale_y_continuous(limits = c(0, 3e-5)) + 
    labs(subtitle = "") + 
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 15))

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
    scale_y_continuous(limits = c(0, 3e-5)) + 
    labs(subtitle = "") + 
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 15))

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
    scale_y_continuous(limits = c(0, 3e-5)) + 
    labs(subtitle = "") + 
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 15))

pRes <- 
    (pF21 + pF23 ) / (pF22 + pF24  )

ggplotResults(pRes, name = "FigS5_indiv_samples", width = 960, height = 480)
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

message("Generating Figure S8...")
p1Medians <- ggplotSampleMDS(
    mdsObj = mdsMedians,
    pData = phenoData,
    projectionAxes = c(1,2),
    pDataForColour = pDataForColour,
    pDataForShape = pDataForShape,
    pDataForLabel = pDataForLabel,
    pDataForAdditionalLabelling = 
        pDataForAdditionalLabelling,
    pointSize = pointSizeMDS,
    title = "MDS - projection axes 1 and 2") + 
    scale_shape_manual(values = c(17, 15)) +
    scale_color_manual(values = c("red", "blue")) + 
    theme_bw() + 
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12))

p2Medians <- ggplotSampleMDS(
    mdsObj = mdsMedians,
    pData = phenoData,
    projectionAxes = c(2,3),
    pDataForColour = pDataForColour,
    pDataForShape = pDataForShape,
    pDataForLabel = pDataForLabel,
    pDataForAdditionalLabelling = 
        pDataForAdditionalLabelling,
    pointSize = pointSizeMDS,
    displayPointLabels = FALSE,
    flipYAxis = TRUE,
    max.overlaps = 100,
    title = "MDS - projection axes 2 and 3") +
    scale_shape_manual(values = c(17, 15)) + 
    scale_color_manual(values = c("red", "blue")) + 
    theme_bw() + 
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12))


pp1 <- p1 + labs(
    title = "CytoMDS - coordinates 1 and 2",
    subtitle = ""
) 

pp2 <- p1Medians + labs(
    title = "MDS with medians - coordinates 1 and 2",
    subtitle = "")

pp3 <- p2 + labs(
    title = "CytoMDS - coordinates 2 and 3",
    subtitle = "")

pp4 <- p2Medians + labs(
    title = "MDS with medians - coordinates 2 and 3",
    subtitle = "")

pp1Th <- pp1 + 
    theme(legend.position="none",
          plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15))

pp2Th <- pp2 + 
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12))
    
pp3Th <- pp3 + 
    theme(legend.position="none",
          plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15))

pp4Th <- pp4 + 
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12))

pRes <- 
    (pp1Th + pp2Th) / (pp3Th + pp4Th)

ggplotResults(pRes, name = "FigS8", width = 960, height = 900)
message("Done!")



