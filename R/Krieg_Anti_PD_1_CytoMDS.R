require(HDCytoData)
require(CytoMDS)
require(ggplot2)
require(patchwork)
require(cyCombine)      # github repo user : biosurf/cyCombine
#require(uwot)          # for umap
#require(irlba)         # for PCA


message("FIGURES - KRIEG ANTI PD1 PART")
message("*****************************")

message("Loading Dataset...")

# load flow Set

Krieg_fs <- Krieg_Anti_PD_1_flowSet()
Krieg_fs

# Krieg_Anti_PD_1 is a mass cytometry (CyTOF) dataset from Krieg et al. (2018), 
# who used mass cytometry to characterize immune cell subsets in 
# peripheral blood from melanoma skin cancer patients treated 
# with anti-PD-1 immunotherapy. 
# This study found that the frequency of CD14+CD16-HLA-DRhi monocytes 
# in baseline samples (taken from patients prior to treatment) was a strong 
# predictor of survival in response to immunotherapy treatment. 
# In particular, the frequency of a small subpopulation of 
# CD14+CD33+HLA-DRhiICAM-1+CD64+CD141+CD86+CD11c+CD38+PD-L1+CD11b+ monocytes 
# in baseline samples was strongly associated with responder status 
# following immunotherapy treatment. 
# Note that this dataset contains a strong batch effect, 
# due to sample acquisition on two different days (Krieg et al., 2018).

nSamples <- length(Krieg_fs)

# update phenoData structure

chLabels <- 
    keyword(Krieg_fs[[1]], "MARKER_INFO")$MARKER_INFO$channel_name
chMarkers <- 
    keyword(Krieg_fs[[1]], "MARKER_INFO")$MARKER_INFO$marker_name
marker_class <- 
    keyword(Krieg_fs[[1]], "MARKER_INFO")$MARKER_INFO$marker_class

chLabels <- chLabels[marker_class != "none"]
chMarkers <- chMarkers[marker_class != "none"]
# marker_class all equal to "type"

phenoData <- flowCore::pData(Krieg_fs)
additionalPhenoData <- 
    keyword(Krieg_fs[[1]], "EXPERIMENT_INFO")$EXPERIMENT_INFO
phenoData <- cbind(phenoData, additionalPhenoData)
    
flowCore::pData(Krieg_fs) <- phenoData


# transform flow set (arcsinh(cofactor = 5))
trans <- arcsinhTransform(
    transformationId="ArcsinhTransform", 
    a = 0, 
    b = 1/5, 
    c = 0)

message("Applying arcsinh() transformation...")
Krieg_fs_trans <- transform(
    Krieg_fs,
    transformList(chLabels, trans))


message("Calculating Sample distances...")
pwDist <- pairwiseEMDDist(
    Krieg_fs_trans,
    channels = chMarkers,
    verbose = FALSE
)

#hist(pwDist)


statFUNs <- c("median" = stats::median,
              "std-dev" = stats::sd,
              "mean" = base::mean,
              "Q10" = function(x, na.rm) 
                  stats::quantile(x, probs = 0.1, na.rm = na.rm),
              "Q20" = function(x, na.rm) 
                  stats::quantile(x, probs = 0.2, na.rm = na.rm),
              "Q80" = function(x, na.rm) 
                  stats::quantile(x, probs = 0.8, na.rm = na.rm),
              "Q90" = function(x, na.rm) 
                  stats::quantile(x, probs = 0.9, na.rm = na.rm))

message("Calculating summary stats...")
chStats <- CytoMDS::channelSummaryStats(
    Krieg_fs_trans,
    channels = chMarkers,
    statFUNs = statFUNs,
    verbose = TRUE)

mds <- computeMetricMDS(
    pwDist,
    seed = 0)

#nDim(mds)
#RSq(mds)

message("ALL DONE ! :-)")

#ggplotSampleMDSShepard(mds)


p1 <- ggplotSampleMDS(
    mds,
    pData = phenoData,
    projectionAxes = c(1,2),
    pDataForColour = "batch_id",
    pDataForLabel = "sample_id",
    repelPointLabels = FALSE
) + theme_bw()

p2 <- ggplotSampleMDS(
    mds,
    pData = phenoData,
    projectionAxes = c(1,2),
    pDataForColour = "group_id",
    pDataForShape = "batch_id",
    pDataForLabel = "sample_id",
    repelPointLabels = FALSE,
    flipXAxis = TRUE,
    title = "Before batch correction") + 
    scale_shape_manual(values = c(17, 15)) +
    scale_color_manual(values = c("red", "blue")) + 
    scale_x_continuous(limits = c(-5, 7)) + 
    scale_y_continuous(limits = c(-5, 7)) + 
    theme_bw()

bp <- ggplotSampleMDS(
    mds,
    biplot = TRUE,
    extVariables = chStats[["median"]],
    pData = phenoData,
    projectionAxes = c(1,2),
    pDataForColour = "group_id",
    pDataForShape = "batch_id",
    pDataForLabel = "sample_id",
    pDataForAdditionalLabelling = c("batch_id"),
    displayPointLabels = FALSE,
    repelArrowLabels = TRUE,
    flipXAxis = TRUE,
    title = "After batch correction: bi-plot with medians") + 
    theme(legend.position = "none") + 
    scale_shape_manual(values = c(17, 15)) + 
    scale_color_manual(values = c("red", "blue")) + 
    scale_x_continuous(limits = c(-5, 7)) + 
    scale_y_continuous(limits = c(-5, 7)) + theme_bw()

#Visualization of channel distributions before batch correction

# pMargDens <- CytoMDS::ggplotMarginalDensities(
#     Krieg_fs_trans,
#     channels = chLabels,
#     pDataForColour = "group_id",
#     pDataForGroup = "sample_id"
# )

# CytoMDS::ggplotMarginalDensities(
#     Krieg_fs_trans,
#     channels = chLabels,
#     pDataForColour = "batch_id",
#     pDataForGroup = "sample_id"
# )

# apply a simple batch correction using cyCombine (Combat) with one label

message("Applying batch correction...")

# convert flowSet to dataframe
uncorrected <- convert_flowset(
    flowset = Krieg_fs_trans,
    metadata = as.data.frame(phenoData),
    filename_col = "name",
    sample_ids = "sample_id",
    batch_ids = "batch_id",
    condition = "group_id",
    down_sample = FALSE,
    clean_colnames = FALSE,
    panel = keyword(Krieg_fs[[1]], "MARKER_INFO")$MARKER_INFO,
    panel_channel = "channel_name",
    panel_antigen = "marker_name"
)

# apply batch correction
corrected <- batch_correct(
    df = uncorrected,
    covar = "condition",
    markers = chMarkers,
    norm_method = "scale",
    xdim = 1,
    ydim = 1,
    rlen = 10,
    seed = 101
)

# norm_ff <- normalize(
#     uncorrected,
#     markers = chMarkers,
#     norm_method = "scale", 
#     ties.method = "average")
# 
# labels <- create_som(
#     norm_ff,
#     markers = chMarkers,
#     rlen = 10, # If results are not convincing, consider using a higher value (e.g. 100)
#     seed = 101,
#     xdim = 1,
#     ydim = 1) # setting xdim & y dim = 1 forces the existence of one! label 
# 
# # Batch correct
# corrected <- correct_data(
#     df = uncorrected,
#     label = labels, # Add custom labels here, if desired
#     covar = "condition",
#     markers = chMarkers,
#     parametric = TRUE
# )

# convert back into flowSet
correctedDF <- as.data.frame(corrected)

ffList <- lapply(
    seq_len(nSamples),
    FUN = function(i) {
        ff <- Krieg_fs_trans[[i]]
        sampleId <- phenoData[i, "sample_id"]
        
        correctedExprs <- flowCore::exprs(ff)
        
        for (j in seq_along(chLabels)) {
            correctedExprCol <- 
                correctedDF[correctedDF$sample == sampleId, chMarkers[j]]
            correctedExprs[, chLabels[j]] <- correctedExprCol
        }
        
        flowCore::exprs(ff) <- correctedExprs
        ff
    }
)

names(ffList) <- rownames(phenoData)

Krieg_fs_corrected <- as(ffList, "flowSet")
flowCore::pData(Krieg_fs_corrected) <- phenoData

# distance calculation after batch correction

message("Calculating sample distances after batch correction...")
pwDistCorr <- pairwiseEMDDist(
    Krieg_fs_corrected,
    channels = chMarkers,
    verbose = FALSE
)

#hist(pwDistCorr)

message("Calculating summary stats after batch corrections...")

chStatsCorr <- CytoMDS::channelSummaryStats(
    Krieg_fs_corrected,
    channels = chMarkers,
    statFUNs = statFUNs,
    verbose = TRUE)

mdsCorr <- computeMetricMDS(
    pwDistCorr,
    seed = 0,
    targetPseudoRSq = 0.95)

#nDim(mdsCorr)
#RSq(mdsCorr)

message("ALL DONE ! :-)")

#ggplotSampleMDSShepard(mdsCorr)


p1Corr <- ggplotSampleMDS(
    mdsCorr,
    pData = phenoData,
    projectionAxes = c(1,2),
    pDataForColour = "batch_id",
    pDataForLabel = "sample_id",
    repelPointLabels = FALSE,
    flipYAxis = TRUE
) + theme_bw()

p2Corr <- ggplotSampleMDS(
    mdsCorr,
    pData = phenoData,
    projectionAxes = c(1,2),
    pDataForColour = "group_id",
    pDataForShape = "batch_id",
    pDataForLabel = "sample_id",
    repelPointLabels = FALSE,
    flipYAxis = TRUE,
    title = "After batch correction") +
    scale_shape_manual(values = c(17, 15)) + 
    scale_color_manual(values = c("red", "blue")) + 
    scale_x_continuous(limits = c(-5, 7)) + 
    scale_y_continuous(limits = c(-5, 7)) + theme_bw()

# Figure 4
message("Generating Figure 4...")
pRes <- 
    (p2 + theme(legend.position = "none")) + p2Corr
ggplotResults(pRes, "Fig4")
message("Done!")



markerSel <- c("CD14", "CD33", "HLA-DR", "ICAM-1", "CD64", 
               "CD141", "CD86", "CD11c", "CD38", 
               "CD274_PDL1", "CD11b")

bpCorr <- ggplotSampleMDS(
    mdsCorr,
    biplot = TRUE,
    extVariables = chStatsCorr[["median"]][, markerSel],
    pData = phenoData,
    projectionAxes = c(1,2),
    pDataForColour = "group_id",
    pDataForShape = "batch_id",
    pDataForLabel = "sample_id",
    pDataForAdditionalLabelling = c("batch_id"),
    displayPointLabels = FALSE,
    displayArrowLabels = TRUE,
    arrowLabelSize = 3,
    repelArrowLabels = TRUE,
    min.segment.length = 0.,
    flipYAxis = TRUE,
    arrowThreshold = 0.7,
    title = "After batch correction: bi-plot with medians") +
    scale_shape_manual(values = c(17, 15)) + 
    scale_color_manual(values = c("red", "blue")) + 
    scale_x_continuous(limits = c(-5, 7)) + 
    scale_y_continuous(limits = c(-5, 7)) + theme_bw()

bpCorr2 <- ggplotSampleMDS(
    mdsCorr,
    biplot = TRUE,
    extVariables = chStatsCorr[["median"]][, !colnames(chStatsCorr[["median"]]) %in% markerSel],
    pData = phenoData,
    projectionAxes = c(1,2),
    pDataForColour = "group_id",
    pDataForShape = "batch_id",
    pDataForLabel = "sample_id",
    pDataForAdditionalLabelling = c("batch_id"),
    displayPointLabels = FALSE,
    displayArrowLabels = TRUE,
    arrowLabelSize = 3,
    repelArrowLabels = TRUE,
    min.segment.length = 0.,
    flipYAxis = TRUE,
    arrowThreshold = 0.,
    title = " ") +
    labs(subtitle = " ") + 
    scale_shape_manual(values = c(17, 15)) + 
    scale_color_manual(values = c("red", "blue")) + 
    scale_x_continuous(limits = c(-5, 7)) + 
    scale_y_continuous(limits = c(-5, 7)) + theme_bw()

# Figure S6

message("Generating Figure S6...")
pRes <- 
    (bpCorr + theme(legend.position = "none")) + bpCorr2
ggplotResults(pRes, "FigS6")
message("Done!")

# Figure S6b (not shown in article)

message("Generating Figure S6b...")
bpUnCorr <- ggplotSampleMDS(
    mds,
    biplot = TRUE,
    extVariables = chStatsCorr[["median"]][, markerSel],
    pData = phenoData,
    projectionAxes = c(1,2),
    pDataForColour = "group_id",
    pDataForShape = "batch_id",
    pDataForLabel = "sample_id",
    pDataForAdditionalLabelling = c("batch_id"),
    displayPointLabels = FALSE,
    displayArrowLabels = TRUE,
    arrowLabelSize = 3,
    repelArrowLabels = TRUE,
    min.segment.length = 0.,
    flipXAxis = TRUE,
    flipYAxis = FALSE,
    arrowThreshold = 0.7,
    title = "Before batch correction: bi-plot with medians") +
    scale_shape_manual(values = c(17, 15)) +
    scale_color_manual(values = c("red", "blue")) +
    scale_x_continuous(limits = c(-5, 7)) +
    scale_y_continuous(limits = c(-5, 7)) + theme_bw()

bpUnCorr2 <- ggplotSampleMDS(
    mds,
    biplot = TRUE,
    extVariables = chStatsCorr[["median"]][, !colnames(chStatsCorr[["median"]]) %in% markerSel],
    pData = phenoData,
    projectionAxes = c(1,2),
    pDataForColour = "group_id",
    pDataForShape = "batch_id",
    pDataForLabel = "sample_id",
    pDataForAdditionalLabelling = c("batch_id"),
    displayPointLabels = FALSE,
    displayArrowLabels = TRUE,
    arrowLabelSize = 3,
    repelArrowLabels = TRUE,
    min.segment.length = 0.,
    flipXAxis = TRUE,
    flipYAxis = FALSE,
    arrowThreshold = 0.,
    title = " ") +
    labs(subtitle = " ") +
    scale_shape_manual(values = c(17, 15)) +
    scale_color_manual(values = c("red", "blue")) +
    scale_x_continuous(limits = c(-5, 7)) +
    scale_y_continuous(limits = c(-5, 7)) + theme_bw()

# Figure S6b

pRes <-
    ((bpUnCorr + labs(title = "Before batch correction") + theme(legend.position = "none")) + bpUnCorr2) /
    ((bpCorr + labs(title = "After batch correction") + theme(legend.position = "none")) + bpCorr2)

ggplotResults(pRes, "FigS6b")
message("Done!")


# Figure 1, two small plots...
message("Generating two small plots part of Figure 1...")
pCorrSmall <- ggplotSampleMDS(
    mdsCorr,
    pData = phenoData,
    projectionAxes = c(1,2),
    pDataForColour = "group_id",
    pDataForShape = "batch_id",
    pDataForLabel = "sample_id",
    repelPointLabels = FALSE,
    flipYAxis = TRUE,
    displayPointLabels = FALSE,
    title = "After batch correction") +
    scale_shape_manual(values = c(17, 15)) + 
    scale_color_manual(values = c("red", "blue")) + 
    scale_x_continuous(limits = c(-5, 7)) + 
    scale_y_continuous(limits = c(-5, 7)) + 
    theme_classic() + 
    theme(axis.text.x=element_blank()) + 
    theme(axis.text.y=element_blank()) +
    labs(x = "axis 1", y = "axis 2")

pCorrSmall

ggplotResults(pCorrSmall, "Fig1_small_MDS")

bpCorrSmall <- ggplotSampleMDS(
    mdsCorr,
    biplot = TRUE,
    extVariables = chStatsCorr[["median"]][, markerSel],
    pData = phenoData,
    projectionAxes = c(1,2),
    pDataForColour = "group_id",
    pDataForShape = "batch_id",
    pDataForLabel = "sample_id",
    pDataForAdditionalLabelling = c("batch_id"),
    displayPointLabels = FALSE,
    displayArrowLabels = FALSE,
    arrowLabelSize = 3,
    repelArrowLabels = TRUE,
    min.segment.length = 0.,
    flipYAxis = TRUE,
    arrowThreshold = 0.,
    title = " ") +
    labs(subtitle = " ") + 
    scale_shape_manual(values = c(17, 15)) + 
    scale_color_manual(values = c("red", "blue")) + 
    scale_x_continuous(limits = c(-5, 7)) + 
    scale_y_continuous(limits = c(-5, 7)) + 
    theme_classic() + 
    theme(axis.text.x=element_blank()) + 
    theme(axis.text.y=element_blank()) +
    labs(x = "axis 1", y = "axis 2")

ggplotResults(bpCorrSmall, "Fig1_small_biplot")

message("Done!")

# UMAP, PCA, t-SNE

# nEvents <- nrow(uncorrected)
# nSamples <- 50000
# 
# # important for reproducibility !! 
# set.seed(0)
# 
# sam <- sample(1:nrow(uncorrected), nSamples)
# 
# start <- Sys.time()
# 
# umapUncorr <- uwot::umap(
#     uncorrected[sam, chMarkers],
#     n_neighbors = 100,
#     min_dist = 0.01,
#     metric = "euclidean",
#     init = "spca",
#     seed = 0,
#     n_epochs = 500,
#     verbose = TRUE)
# 
# end <- Sys.time()
# 
# print( difftime(end, start)) 
# 
# start <- Sys.time()
# 
# tsneUncorr <- Rtsne::Rtsne(
#     X = uncorrected[sam, chMarkers],
#     perplexity = 100,
#     verbose = TRUE
# )
# 
# end <- Sys.time()
# 
# print( difftime(end, start)) 
# 
# pcaUncorr <- irlba::prcomp_irlba(
#     x = uncorrected[sam, chMarkers],
#     n = 2)
# 
# DFUncorr <- cbind(uncorrected[sam, ], 
#                   data.frame(UC1 = umapUncorr[,1],
#                              UC2 = umapUncorr[,2],
#                              TC1 = tsneUncorr$Y[, 1],
#                              TC2 = tsneUncorr$Y[, 2],
#                              PC1 = pcaUncorr$x[,1],
#                              PC2 = pcaUncorr$x[,2]))
# 
# DFUncorr$cond_bat <- as.factor(
#     paste(DFUncorr$condition, DFUncorr$batch, sep = "_"))
# 
# start <- Sys.time()
# 
# umapCorr <- uwot::umap(
#     corrected[sam, chMarkers],
#     n_neighbors = 100,
#     min_dist = 0.01,
#     metric = "euclidean",
#     init = "spca",
#     seed = 0,
#     n_epochs = 500,
#     verbose = TRUE)
# 
# end <- Sys.time()
# 
# print( difftime(end, start)) 
# 
# start <- Sys.time()
# 
# tsneCorr <- Rtsne::Rtsne(
#     X = corrected[sam, chMarkers],
#     perplexity = 100,
#     verbose = TRUE
# )
# 
# end <- Sys.time()
# 
# print( difftime(end, start)) 
# 
# pcaCorr <- irlba::prcomp_irlba(
#     x = corrected[sam, chMarkers],
#     n = 2)
# 
# DFCorr <- cbind(corrected[sam, ], 
#                 data.frame(UC1 = umapCorr[,1],
#                            UC2 = umapCorr[,2],
#                            TC1 = tsneCorr$Y[, 1],
#                            TC2 = tsneCorr$Y[, 2],
#                            PC1 = pcaCorr$x[,1],
#                            PC2 = pcaCorr$x[,2]))
# 
# DFCorr$cond_bat <- as.factor(
#     paste(DFCorr$condition, DFCorr$batch, sep = "_"))
# 
# 
# #saveRDS(DFUncorr, "RDS/Krieg_uncorr_tsne_results.rds")
# #saveRDS(DFCorr, "RDS/Krieg_corr_tsne_results.rds")
# DFUncorr <- readRDS("RDS/Krieg_uncorr_tsne_results.rds")
# DFCorr <- readRDS("RDS/Krieg_corr_tsne_results.rds")
# 
# pUncorrUMAP <- ggplot(DFUncorr, 
#                       mapping = aes(x = UC1, y = UC2, colour = batch)) + 
#     geom_point(size = 0.5) + 
#     labs(title = "UMAP before batch correction")
#     
# pUncorrTSNE <- ggplot(DFUncorr, 
#                       mapping = aes(x = TC1, y = TC2, colour = batch)) + 
#     geom_point(size = 0.5) + 
#     labs(title = "TSNE before batch correction") 
# 
# pUncorrPCA <- ggplot(DFUncorr, 
#                      mapping = aes(x = PC1, y = PC2, colour = batch)) + 
#     geom_point(size = 0.5) + 
#     labs(title = "PCA before batch correction")
# 
# pCorrUMAP <- ggplot(DFCorr, 
#                     mapping = aes(x = UC1, y = UC2, colour = batch)) + 
#     geom_point(size = 0.5) + 
#     labs(title = "UMAP after batch correction") 
# 
# pCorrTSNE <- ggplot(DFCorr, 
#                     mapping = aes(x = TC1, y = TC2, colour = batch)) + 
#     geom_point(size = 0.5) + 
#     labs(title = "TSNE after batch correction") 
# 
# pCorrPCA <- ggplot(DFCorr, 
#                    mapping = aes(x = PC1, y = PC2, colour = batch)) + 
#     geom_point(size = 0.5) + 
#     labs(title = "PCA after batch correction")
# 
# ( (pUncorrUMAP + theme(legend.position = "none")) + 
#   (pUncorrTSNE + theme(legend.position = "none")) + 
#   pUncorrPCA ) /
# ( (pCorrUMAP + theme(legend.position = "none")) + 
#   (pCorrTSNE + theme(legend.position = "none")) + 
#   pCorrPCA )
# 
# ggsave("plots/Krieg_umap_tsne_pca.pdf")

