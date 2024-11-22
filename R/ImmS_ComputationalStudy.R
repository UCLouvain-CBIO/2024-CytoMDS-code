message("RUN COMPUTATIONAL PERFORMANCE BENCHMARK")
message("***************************************")

require(CytoPipeline)
require(CytoMDS)
require(ggplot2)
require(patchwork)
require(purrr)

# Set the following flag to TRUE to run the computational benchmarking
# WARNING: This is likely to take between one and two hours on a standard computer
computeStudy <- FALSE
resultFileMarker <- "./rds/resMarker.rds"
resultFileSampling <- "./rds/resSampling.rds"

if (computeStudy) {
    
    message("Read two compensated samples from pre-processing store...")
    
    # directory from which we read the preprocessing results
    prepDir <- "./preprocessing"
    
    expName <- "ImmunoSenescence_human_PBMC"
    
    prepstep <- "compensate"
    pipelineObj <- paste0(prepstep, "_obj")
    
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
    allChannels <- flowCore::colnames(inputFF)[c(2,5,8:33)]
    
    # reading target flowSet
    
    allPhenoData <- CytoPipeline::pData(pipL)
   
    
    tgtSeq <- c(1, 9)
    ffList <- list()
    for(i in seq_along(tgtSeq)) {
        message("Reading file #", i, "/", length(tgtSeq), "...")
        ffList[[i]] <- getCytoPipelineFlowFrame(
            pipL,
            whichQueue = "pre-processing",
            sampleFile = tgtSeq[i],
            objectName = pipelineObj,
            path = prepDir)
    }
    
    names(ffList) <- rownames(allPhenoData)[tgtSeq]
    
    fsBoth <- as(ffList,"flowSet")
    
    flowCore::pData(fsBoth) <- allPhenoData[tgtSeq,]
    
    transList <- CytoPipeline::getCytoPipelineScaleTransform(
        pipL,
        whichQueue = "scale transform",
        objectName = "scale_transform_estimate_obj",
        path = prepDir)
    
    message("applying scale transformations to both flow frames...")
    fsBoth <- CytoPipeline::applyScaleTransforms(
        fsBoth,
        transList = transList,
        verbose = TRUE
    )
    
    NCells <- 10000000
    rNoiseSD <- 0.001
    message(paste0("Sample ", NCells, " events from the two flow frames ", 
                   "(with random noise)..."))
    
    # important for reproducibility !! 
    set.seed(0)
    
    expr1 <- flowCore::exprs(fsBoth[[1]])[, allChannels]
    expr1 <- expr1[sample(nrow(expr1), NCells, replace = TRUE), ]
    expr2 <- flowCore::exprs(fsBoth[[2]])[, allChannels]
    expr2 <- expr2[sample(nrow(expr2), NCells, replace = TRUE), ]
    
    expr1 <- expr1 + matrix(data = rnorm(NCells * ncol(expr1), sd = rNoiseSD),
                            ncol = ncol(expr1))
    expr2 <- expr2 + matrix(data = rnorm(NCells * ncol(expr1), sd = rNoiseSD),
                            ncol = ncol(expr1))
    
    # startTime <- Sys.time()
    # refDist <- CytoMDS::EMDDist(x1 = expr1, x2 = expr2)
    # endTime <- Sys.time()
    # print(endTime - startTime)
    # 
    # refDist

    calculateEMDWithNbMarkers <- function(nMarker, expr1, expr2, markers) {
        message(paste0("calculating distance with ", nMarker, " markers..."))
        selMarkers <- markers[1:nMarker]
        exprA <- expr1[, selMarkers, drop = FALSE]
        exprB <- expr2[, selMarkers, drop = FALSE]
        startTime <- Sys.time()
        dis <- CytoMDS::EMDDist(x1 = exprA,
                                x2 = exprB)
        endTime <- Sys.time()
        dd = endTime - startTime
        message(paste0("It took ", as.numeric(dd, units = "secs"), " seconds"))
        list(nMarker = nMarker,
             dis = dis,
             elapsed = as.numeric(dd, units = "secs"))
    }

    calculateEMDWithSampling <- function(nSampling, expr1, expr2) {
        message(paste0("calculating distance with ", nSampling, " samplings..."))
        samples1 <- sample(nrow(expr1), nSampling, replace = FALSE)
        samples2 <- sample(nrow(expr2), nSampling, replace = FALSE)
        exprA <- expr1[samples1, , drop = FALSE]
        exprB <- expr2[samples2, , drop = FALSE]
        startTime <- Sys.time()
        dis <- CytoMDS::EMDDist(
            x1 = exprA,
            x2 = exprB)
        endTime <- Sys.time()
        dd = endTime - startTime
        message(paste0("It took ", as.numeric(dd, units = "secs"), " seconds"))
        list(nSampling = nSampling,
             dis = dis,
             elapsed = as.numeric(dd, units = "secs"))
    }

    message(paste0("Running distance calculation with a varying nb of markers..."))

    nMarkers <- seq(27, 3, -3)
    
    nRuns <- 20
    
    expr1_100 <- expr1[sample(nrow(expr1), 100000, replace = FALSE),]
    expr2_100 <- expr1[sample(nrow(expr1), 100000, replace = FALSE),]

    allRes1 <- lapply(1:nRuns,
                      FUN = function(run) {
                          message("STARTING RUN #", run, "/", nRuns)
                          message("*******************************")
                          markers <- sample(
                              ncol(expr1_100),
                              size = ncol(expr1_100),
                              replace = FALSE)
                          runResults <- lapply(nMarkers,
                                               FUN = calculateEMDWithNbMarkers,
                                               expr1 = expr1_100,
                                               expr2 = expr2_100,
                                               markers = markers)
                          list(run = run, res = runResults)
                      })

    message("ALL DONE!")

    message(paste0("Running distance calculation with various sub-sampling..."))

    #nSamplings <- c(100, 1000)
    nSamplings <- c(10000000, 1000000, 100000, 10000, 1000, 100)
    nRuns <- 20

    allRes2 <- lapply(1:nRuns,
                      FUN = function(run) {
                          message("STARTING RUN #", run, "/", nRuns)
                          message("*******************************")
                          runResults <- lapply(nSamplings,
                                               FUN = calculateEMDWithSampling,
                                               expr1 = expr1,
                                               expr2 = expr2)
                          list(run = run, res = runResults)
                      })
    message("ALL DONE!")
    
    saveRDS(allRes1, file = resultFileMarker)
    saveRDS(allRes2, file = resultFileSampling)
    
} else {
    allRes1 <- readRDS(file = resultFileMarker)
    allRes2 <- readRDS(file = resultFileSampling)
}


DFSampling <- data.frame(
    run = purrr::map_int(allRes2, "run"),
    nSampling = unlist(lapply(allRes2, function(el) {
        purrr::map_int(el$res, "nSampling")})),
    dis = unlist(lapply(allRes2, function(el) {
        purrr::map_dbl(el$res, "dis")})),
    elapsed = unlist(lapply(allRes2, function(el) {
        purrr::map_dbl(el$res, "elapsed")})))

DFSampling$nSamplingFct <- factor(DFSampling$nSampling)

visSampling1 <- c(100, 1000, 10000, 100000, 1000000, 10000000)

pSampling1 <- ggplot(DFSampling[DFSampling$nSampling %in% visSampling1, ],
                     mapping = aes(x = nSamplingFct, y = elapsed)) +
    geom_boxplot() +
    scale_x_discrete(
        labels = format(visSampling1, big.mark = ",", scientific = FALSE)) + 
    scale_y_continuous(trans='log',
                       limits = c(0.05, 100),
                       breaks = c(0.06, 0.6, 6, 60)) + 
    xlab("number of sub-sampled events") +
    ylab("CPU time (s)") +
    theme_bw() +
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15, angle = 90),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12))

visSampling2 <- c(100, 1000, 10000, 100000, 1000000, 10000000)

pSampling2 <- ggplot(DFSampling[DFSampling$nSampling %in% visSampling2, ],
                     mapping = aes(x = nSamplingFct, y = dis)) +
    geom_boxplot() +
    scale_x_discrete(
        labels = format(visSampling2, big.mark = ",", scientific = FALSE)) +
    scale_y_continuous(limits = c(0, 10)) +
    xlab("number of sub-sampled events") +
    ylab("computed EMD") +
    theme_bw() +
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15, angle = 90),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12))

ggplotResults(pSampling1 + pSampling2,
              name = "FigS9_comp_study_events", width = 960, height = 480)

DFMarkers <- data.frame(
    run = purrr::map_int(allRes1, "run"),
    nMarker = unlist(lapply(allRes1, function(el) {
        purrr::map_int(el$res, "nMarker")})),
    dis = unlist(lapply(allRes1, function(el) {
        purrr::map_dbl(el$res, "dis")})),
    elapsed = unlist(lapply(allRes1, function(el) {
        purrr::map_dbl(el$res, "elapsed")})))

DFMarkers$nMarkerFct <- factor(DFMarkers$nMarker)

DFMarkers <- DFMarkers[DFMarkers$nMarker %in% c(3,9,15,21,27),]

pMarker <- ggplot(DFMarkers, mapping = aes(x = nMarkerFct, y = elapsed)) +
    geom_boxplot() +
    #labs(title = "CPU time vs. number of markers") +
    xlab("number of markers") +
    ylab("CPU time (s)") +
    theme_bw() +
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12))

ggplotResults(pMarker, name = "FigS10_comp_study_markers", width = 480, 
              height = 480)

message("Computational study done!")







