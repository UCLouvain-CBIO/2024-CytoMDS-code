plotDir <- file.path(getwd(), "plots")
dir.create(plotDir, showWarnings = FALSE)

plotInFile <- TRUE
plotFileRoot <- file.path(plotDir, 
                          "Plot_")

ggplotResults <- function(ggplotObj, name) {
    if (plotInFile) {
        png(filename = paste0(plotFileRoot, name, ".png"))
        print(ggplotObj)
        dev.off()
        invisible(ggplotObj)   
    } else {
        ggplotObj    
    }
}
