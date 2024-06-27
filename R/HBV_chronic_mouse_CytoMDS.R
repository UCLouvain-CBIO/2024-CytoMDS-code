message("FIGURES - HBV CHRONIC MOUSE PART")
message("********************************")

require(CytoPipeline)
require(CytoPipelineGUI)
require(CytoMDS)
require(ggplot2)
require(patchwork)
require(Polychrome)
require(cyCombine)      # github repo user : biosurf/cyCombine
#require(uwot)          # for umap
#require(irlba)         # for PCA



message("Getting Results...")

# directory from which we read the preprocessing results
prepDir <- "./preprocessing"

expName <- "HBV_chronic_mouse"

prepstep <- "compensate"
pipelineObj <- paste0(prepstep, "_obj")

## marker selection: either pre-processing channels or all 
# channelTypes <- c("prep", "all")
channelTypes <- "prep"
selChannelList <- setNames(vector("list", 2), channelTypes)
selChannelList[["prep"]] <- c("FSC-A", "SSC-A", "FSC-H", "Live & Dead")

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

selChannelList[["all"]] <- flowCore::colnames(inputFF)[c(2:14)]
flowCore::pData(flowCore::parameters(inputFF))

selChannelList[["all"]] <- c(
    "FSC-A", "FSC-H", "SSC-A", 
    "CD4", "CD8", "CD38", "CD44",
    "Marker1", "Marker2", "Marker3", "Marker4", "Marker5"
)

# reading target flowSet
sampleFiles <- sampleFiles(pipL)

allPhenoData <- CytoPipeline::pData(pipL)

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

verbose <- FALSE

resultArray <- vapply(
    channelTypes,
    USE.NAMES = TRUE,
    FUN.VALUE = vector("list", 4),
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
        
        list(phenoData = allPhenoData,
             chStats = chStats,
             pwDist = pwDist,
             mds = mds)
    }
)

message("ALL DONE ! :-)")

## plots for preprocessing channels only


selChType <- "prep"

phenoData <- resultArray["phenoData", selChType][[1]]
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

# nDim(mds)
# CytoMDS::RSq(mds)

# Figure S1

message("Generating Figure S1...")
pSh <- ggplotSampleMDSShepard(mds, 
                               lineWidth = 1.0, 
                               pointSize = 1.0,
                               title = "") + 
    scale_x_continuous(limits = c(0,2)) + 
    theme_bw()

ggplotResults(pSh, name = "FigS1")
message("Done!")

# Figure 2

message("Generating 3 parts of Figure 2...")
pDataForShape <- "group"
pDataForLabel <- "well"
pDataForShape <- "day"
pDataForColour <- "day"

p1 <- ggplotSampleMDS(
    mdsObj = mds,
    pData = phenoData,
    projectionAxes = c(1,2),
    pDataForColour = pDataForColour,
    pDataForShape = pDataForShape,
    pDataForLabel = pDataForLabel,
    max.overlaps = 100,
    title = "MDS - coordinates 1 and 2"
) + scale_shape_manual(values=c(17, 15))+ 
    scale_color_manual(values = c("red", "blue")) + 
    theme_bw()

ggplotResults(p1, name = "Fig2_main_part")

b1 <- ggplotSampleMDS(
    mdsObj = mds,
    pData = phenoData,
    projectionAxes = c(1,2),
    biplot = TRUE,
    extVariables = chStats[["median"]],
    pDataForColour = pDataForColour,
    pDataForShape = pDataForShape,
    pDataForLabel = pDataForLabel,
    displayPointLabels = FALSE,
    arrowThreshold = 0.5,
    repelArrowLabels = TRUE,
    title = "Bi-plot with medians" 
    ) + theme(legend.position = "none") + 
    scale_shape_manual(values=c(17, 15)) + 
    scale_color_manual(values = c("red", "blue")) + 
    theme_bw()


b2 <- ggplotSampleMDS(
    mdsObj = mds,
    pData = phenoData,
    projectionAxes = c(1,2),
    biplot = TRUE,
    extVariables = chStats[["std-dev"]],
    pDataForColour = pDataForColour,
    pDataForShape = pDataForShape,
    pDataForLabel = pDataForLabel,
    displayPointLabels = FALSE,
    arrowThreshold = 0.5,
    repelArrowLabels = TRUE,
    title = "Bi-plot with standard deviations" 
) + theme(legend.position = "none") + 
    scale_shape_manual(values=c(17, 15))+ 
    scale_color_manual(values = c("red", "blue")) + 
    theme_bw()

ggplotResults(b1 / b2, name = "Fig2_biplots")


#CytoPipelineGUI::CytoPipelineCheckApp(dir = prepDir)


pF1 <- CytoPipelineGUI::plotSelectedFlowFrame(
    experimentName = expName,
    whichQueue = "pre-processing",
    sampleFile = "D93_G03.fcs",
    flowFrameName = "compensate_obj",
    path = prepDir,
    xChannelLabel = "FSC-A : NA",
    yChannelLabel = "SSC-A : NA",
    useAllCells = FALSE,
    nDisplayCells = 10000,
    useFixedLinearRange = TRUE,
    linearRange = c(-300, 262144),
    transfoListName = " ")

pF2 <- CytoPipelineGUI::plotSelectedFlowFrame(
    experimentName = expName,
    whichQueue = "pre-processing",
    sampleFile = "D93_A05.fcs",
    flowFrameName = "compensate_obj",
    path = prepDir,
    xChannelLabel = "FSC-A : NA",
    yChannelLabel = "SSC-A : NA",
    useAllCells = FALSE,
    nDisplayCells = 10000,
    useFixedLinearRange = TRUE,
    linearRange = c(-300, 262144),
    transfoListName = " ")

ggplotResults(pF1 + pF2, name = "Fig2_indiv_samples" )

message("Done!")

# Figure S2

message("Generating Figure S2...")
pFacet1 <- CytoPipeline::ggplotEvents(
    obj = fsAll[1:27],
    xChannel = "FSC-A",
    yChannel = "SSC-A",
    nDisplayCells = 10000,
    seed = 0,
    xLinearRange = c(0, 2.),
    yLinearRange = c(0, 2.)) + facet_wrap(~name, ncol = 3)

pFacet2 <- CytoPipeline::ggplotEvents(
    obj = fsAll[28:55],
    xChannel = "FSC-A",
    yChannel = "SSC-A",
    nDisplayCells = 10000,
    seed = 0,
    xLinearRange = c(0, 2.),
    yLinearRange = c(0, 2.)) + facet_wrap(~name, ncol = 3)

ggplotResults(pFacet1 + pFacet2, name = "FigS2" )
message("Done!")

# Figure S1b (not shown in article)

message("Generating Figure S1b...")
pF3 <- CytoPipelineGUI::plotSelectedFlowFrame(
    experimentName = expName,
    whichQueue = "pre-processing",
    sampleFile = "D93_G03.fcs",
    flowFrameName = "compensate_obj",
    path = prepDir,
    xChannelLabel = "FSC-A : NA",
    yChannelLabel = " ",
    useAllCells = TRUE,
    useFixedLinearRange = TRUE,
    linearRange = c(-300, 262144),
    transfoListName = " ") + labs(subtitle = "") + 
    scale_y_continuous(limits = c(0, 0.0001))

pF4 <- CytoPipelineGUI::plotSelectedFlowFrame(
    experimentName = expName,
    whichQueue = "pre-processing",
    sampleFile = "D93_A05.fcs",
    flowFrameName = "compensate_obj",
    path = prepDir,
    xChannelLabel = "FSC-A : NA",
    yChannelLabel = " ",
    useAllCells = TRUE,
    useFixedLinearRange = TRUE,
    linearRange = c(-300, 262144),
    transfoListName = " ") + labs(subtitle = "") + 
    scale_y_continuous(limits = c(0, 0.0001))

pF5 <- CytoPipelineGUI::plotSelectedFlowFrame(
    experimentName = expName,
    whichQueue = "pre-processing",
    sampleFile = "D93_G03.fcs",
    flowFrameName = "compensate_obj",
    path = prepDir,
    xChannelLabel = "SSC-A : NA",
    yChannelLabel = " ",
    useAllCells = TRUE,
    useFixedLinearRange = TRUE,
    linearRange = c(-300, 262144),
    transfoListName = " ") + labs(subtitle = "") + 
    scale_y_continuous(limits = c(0, 0.00005))

pF6 <- CytoPipelineGUI::plotSelectedFlowFrame(
    experimentName = expName,
    whichQueue = "pre-processing",
    sampleFile = "D93_A05.fcs",
    flowFrameName = "compensate_obj",
    path = prepDir,
    xChannelLabel = "SSC-A : NA",
    yChannelLabel = " ",
    useAllCells = TRUE,
    useFixedLinearRange = TRUE,
    linearRange = c(-300, 262144),
    transfoListName = " ") + labs(subtitle = "") + 
    scale_y_continuous(limits = c(0, 0.00005))

pF7 <- CytoPipelineGUI::plotSelectedFlowFrame(
    experimentName = expName,
    whichQueue = "pre-processing",
    sampleFile = "D93_G03.fcs",
    flowFrameName = "compensate_obj",
    path = prepDir,
    xChannelLabel = "Comp-APC-Cy7-A : Live & Dead",
    yChannelLabel = " ",
    useAllCells = TRUE,
    useFixedLinearRange = TRUE,
    linearRange = c(-300, 262144),
    transfoListName = "scale_transform_estimate_obj") + labs(subtitle = "") + 
    scale_y_continuous(limits = c(0, 3))

pF8 <- CytoPipelineGUI::plotSelectedFlowFrame(
    experimentName = expName,
    whichQueue = "pre-processing",
    sampleFile = "D93_A05.fcs",
    flowFrameName = "compensate_obj",
    path = prepDir,
    xChannelLabel = "Comp-APC-Cy7-A : Live & Dead",
    yChannelLabel = " ",
    useAllCells = TRUE,
    useFixedLinearRange = TRUE,
    linearRange = c(-300, 262144),
    transfoListName = "scale_transform_estimate_obj") + labs(subtitle = "") + 
    scale_y_continuous(limits = c(0, 3))




layout <- "
AABB
CCDD
EEFF
"
pRes <- 
    pF3 + pF4 +
    pF5 + pF6 + 
    pF7 + pF8 + 
    plot_layout(design = layout) 

ggplotResults(pRes, name = "FigS1b" )
message("Done!")
    

# Figure S1c (not shown in article)

message("Generating Figure S1c...")
file1 <- "D93_B09.fcs"
file2 <- "D93_C01.fcs"

pF9 <- CytoPipelineGUI::plotSelectedFlowFrame(
    experimentName = expName,
    whichQueue = "pre-processing",
    sampleFile = file1,
    flowFrameName = "compensate_obj",
    path = prepDir,
    xChannelLabel = "FSC-A : NA",
    yChannelLabel = " ",
    useAllCells = TRUE,
    useFixedLinearRange = TRUE,
    linearRange = c(-300, 262144),
    transfoListName = " ") + labs(subtitle = "") + 
    scale_y_continuous(limits = c(0, 0.0002))

pF10 <- CytoPipelineGUI::plotSelectedFlowFrame(
    experimentName = expName,
    whichQueue = "pre-processing",
    sampleFile = file2,
    flowFrameName = "compensate_obj",
    path = prepDir,
    xChannelLabel = "FSC-A : NA",
    yChannelLabel = " ",
    useAllCells = TRUE,
    useFixedLinearRange = TRUE,
    linearRange = c(-300, 262144),
    transfoListName = " ") + labs(subtitle = "") + 
    scale_y_continuous(limits = c(0, 0.0002))

pF11 <- CytoPipelineGUI::plotSelectedFlowFrame(
    experimentName = expName,
    whichQueue = "pre-processing",
    sampleFile = file1,
    flowFrameName = "compensate_obj",
    path = prepDir,
    xChannelLabel = "SSC-A : NA",
    yChannelLabel = " ",
    useAllCells = TRUE,
    useFixedLinearRange = TRUE,
    linearRange = c(-300, 262144),
    transfoListName = " ") + labs(subtitle = "") + 
    scale_y_continuous(limits = c(0, 0.00006))

pF12 <- CytoPipelineGUI::plotSelectedFlowFrame(
    experimentName = expName,
    whichQueue = "pre-processing",
    sampleFile = file2,
    flowFrameName = "compensate_obj",
    path = prepDir,
    xChannelLabel = "SSC-A : NA",
    yChannelLabel = " ",
    useAllCells = TRUE,
    useFixedLinearRange = TRUE,
    linearRange = c(-300, 262144),
    transfoListName = " ") + labs(subtitle = "") + 
    scale_y_continuous(limits = c(0, 0.00006))

pF13 <- CytoPipelineGUI::plotSelectedFlowFrame(
    experimentName = expName,
    whichQueue = "pre-processing",
    sampleFile = file1,
    flowFrameName = "compensate_obj",
    path = prepDir,
    xChannelLabel = "Comp-APC-Cy7-A : Live & Dead",
    yChannelLabel = " ",
    useAllCells = TRUE,
    useFixedLinearRange = TRUE,
    linearRange = c(-300, 262144),
    transfoListName = "scale_transform_estimate_obj") + labs(subtitle = "") + 
    scale_y_continuous(limits = c(0, 1))

pF14 <- CytoPipelineGUI::plotSelectedFlowFrame(
    experimentName = expName,
    whichQueue = "pre-processing",
    sampleFile = file2,
    flowFrameName = "compensate_obj",
    path = prepDir,
    xChannelLabel = "Comp-APC-Cy7-A : Live & Dead",
    yChannelLabel = " ",
    useAllCells = TRUE,
    useFixedLinearRange = TRUE,
    linearRange = c(-300, 262144),
    transfoListName = "scale_transform_estimate_obj") + labs(subtitle = "") + 
    scale_y_continuous(limits = c(0, 1))

pRes <- 
    pF9 + pF10 +
    pF11 + pF12 + 
    pF13 + pF14 + 
    plot_layout(design = layout) 


ggplotResults(pRes, name = "FigS1c" )
message("Done!")

# UMAP, TSNE, PCA
message("Calculating tSNE data for Fig S7...")
panel <- flowCore::pData(flowCore::parameters(fsAll[[1]]))[, c("name", "desc")]
colnames(panel) <- c("channel", "marker")


fsAllSub <- flowCore::fsApply(
    fsAll,
    FUN = CytoPipeline::subsample,
    nEvents = 10000,
    seed = 0
)


# convert flowSet to dataframe
uncorrected <- convert_flowset(
    flowset = fsAllSub,
    metadata = as.data.frame(phenoData),
    filename_col = "name",
    sample_ids = "name",
    batch_ids = "day",
    condition = "group",
    down_sample = FALSE,
    clean_colnames = FALSE,
    panel = panel,
    panel_channel = "channel",
    panel_antigen = "marker"
)

colNames <- colnames(uncorrected)
colNames[2:5] <- panel$channel[1:4]
colNames[16:19] <- c("Live/Dead", "sample", "day", "group")
colnames(uncorrected) <- colNames

#table(uncorrected$sample)
lowQuality <- 
    c("D93_C01", "D93_D05", "D93_A01", "D93_D07", "D93_E05",
      "D93_B01", "D93_F03", "D93_G05", "D93_E03", "D93_B03",
      "D93_A05", "D93_C05", "D93_D01", "D93_D03", "D93_B07",
      "D93_A07", "D93_C07", "D93_C03", "D93_A09", "D93_B09")

unClear <- 
    c("D93_E07", "D93_F07", "D91_F05", "D91_D03", "D91_A03")


uncorrected$quality <- factor(ifelse(
    uncorrected$sample %in% lowQuality, "low",
    "ok"),
    levels = c("ok", "low"))


nEvents <- nrow(uncorrected)
nSamples <- 50000

# important for reproducibility !! 
set.seed(0)

sam <- sample(1:nrow(uncorrected), nSamples)

chMarkers <- selChannelList[["all"]]

#run tsne
message("Calculating tSNE...")
start <- Sys.time()

tsne <- Rtsne::Rtsne(
    X = uncorrected[sam, chMarkers],
    perplexity = 100,
    verbose = TRUE
)

end <- Sys.time()

message("Done!")
#message(difftime(end, start)) 


# RUN umap
# start <- Sys.time()
# 
# umap <- uwot::umap(
#     uncorrected[sam, chMarkers],
#     n_neighbors = 100,
#     min_dist = 0.01,
#     metric = "euclidean",
#     init = "spca",
#     seed = 0,
#     n_epochs = 1000,
#     verbose = TRUE)
#                    
# end <- Sys.time()
# 
# print( difftime(end, start)) 

# RUN PCA
# 
# pca <- irlba::prcomp_irlba(x = uncorrected[sam, chMarkers],
#                            n = 2)
# 
# pcaXLabel <- paste0(
#     "PC1 (",
#     round(100 * pca$sdev[1] / pca$totalvar, 2),
#     "% var.)")
# 
# pcaYLabel <- paste0(
#     "PC2 (",
#     round(100 * pca$sdev[2] / pca$totalvar, 2),
#     "% var.)")


DF <- cbind(uncorrected[sam, ], 
            data.frame(
                #UC1 = umap[, 1],
                #UC2 = umap[, 2],
                TC1 = tsne$Y[, 1],
                TC2 = tsne$Y[, 2]#,
                #PC1 = pca$x[, 1],
                #PC2 = pca$x[, 2])
            ))

# build-in color palette
P55 = setNames(
    createPalette(55,  c("#ff0000", "#00ff00", "#0000ff")),
    nm = unique(DF$sample))

pTSNE1 <- ggplot(DF, 
                mapping = aes(x = TC1, y = TC2, colour = day)) + 
    geom_point(size = 0.5) + 
    labs(title = "TSNE") + 
    scale_color_manual(values = c("red", "blue")) 

# pUMAP1 <- ggplot(DF, 
#                  mapping = aes(x = UC1, y = UC2, colour = day)) + 
#     geom_point(size = 0.5) +
#     labs(title = "UMAP") + 
#     scale_color_manual(values = c("red", "blue")) 

# pPCA1 <- ggplot(DF, 
#                mapping = aes(x = PC1, y = PC2, colour = day)) + 
#     geom_point(size = 0.5) + 
#     labs(title = "PCA", x = pcaXLabel, y = pcaYLabel) + 
#     scale_color_manual(values = c("red", "blue")) 

pTSNE2 <- ggplot(DF, 
                 mapping = aes(x = TC1, y = TC2, colour = sample)) + 
    geom_point(size = 0.5) + 
    labs(title = "TSNE") + 
    scale_color_manual(values = P55) + 
    theme(legend.key.size = unit(0.3, 'cm'),
          legend.text = element_text(size = 8)) 

# pUMAP2 <- ggplot(DF, 
#                  mapping = aes(x = UC1, y = UC2, colour = sample)) + 
#     geom_point(size = 0.5) +
#     labs(title = "UMAP") + 
#     scale_color_manual(values = P55) + 
#     theme(legend.key.size = unit(0.3, 'cm'),
#           legend.text = element_text(size = 8)) 
# 
# pPCA2 <- ggplot(DF, 
#                 mapping = aes(x = PC1, y = PC2, colour = sample)) + 
#     geom_point(size = 0.5) + 
#     labs(title = "PCA", x = pcaXLabel, y = pcaYLabel) + 
#     scale_color_manual(values = P55) + 
#     theme(legend.key.size = unit(0.3, 'cm'),
#           legend.text = element_text(size = 8)) 



pTSNE3 <- ggplot(DF, 
                 mapping = aes(x = TC1, y = TC2, colour = quality)) + 
    geom_point(size = 0.5) +
    labs(title = "TSNE",
         colour = "sample quality") + 
    scale_color_manual(values = c("red", "blue")) 

# pUMAP3 <- ggplot(DF, 
#                  mapping = aes(x = UC1, y = UC2, colour = quality)) + 
#     geom_point(size = 0.5) +
#     labs(title = "UMAP",
#          colour = "sample quality") + 
#     scale_color_manual(values = c("red", "blue")) 
# 
# pPCA3 <- ggplot(DF, 
#                 mapping = aes(x = PC1, y = PC2, colour = quality)) + 
#     geom_point(size = 0.5) +
#     labs(title = "PCA", x = pcaXLabel, y = pcaYLabel, 
#          colour = "sample quality") + 
#     scale_color_manual(values = c("red", "blue")) 


# layout <- "
# AABBCC
# DDEEFF
# GGHHII
# "
# 
# (bigPlot <- 
#     (pUMAP1 + theme(legend.position = "none")) +
#     (pTSNE1 + theme(legend.position = "none")) +
#     pPCA1 + 
#     (pUMAP2 + theme(legend.position = "none")) + 
#     (pTSNE2 + theme(legend.position = "none")) + 
#     pPCA2 + 
#     (pUMAP3 + theme(legend.position = "none")) + 
#     (pTSNE3 + theme(legend.position = "none")) + 
#     pPCA3 + plot_layout(design = layout))
# 

# Figure S7
message("Generating Fig S7...")
pRes <- 
    (pTSNE1 + labs(title = "")) / (pTSNE2 + labs(title = "")) / 
    (pTSNE3 + labs(title = ""))

ggplotResults(pRes, name = "FigS7" )
message("Done!")



