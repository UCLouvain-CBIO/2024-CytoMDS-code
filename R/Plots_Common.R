plotDir <- file.path(getwd(), "plots")
dir.create(plotDir, showWarnings = FALSE)

plotInFile <- TRUE

ggplotResults <- function(ggplotObj, name, ...) {
    if (plotInFile) {
        png(filename = file.path(plotDir, 
                                 paste0(name, ".png")), ...)
        print(ggplotObj)
        dev.off()
        invisible(ggplotObj)   
    } else {
        ggplotObj    
    }
}
