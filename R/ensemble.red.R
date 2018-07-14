`ensemble.red` <- function(
    x=NULL
)
{
# Function inspired from red::map.sdm function
# However, AOO and EOO calculated for different thresholds of 'count' suitability maps
# AOO and EOO thresholds from http://www.iucnredlist.org/static/categories_criteria_3_1 (April 2018)
#
#   if (! require(dismo)) {stop("Please install the dismo package")}
    if (! requireNamespace("red")) {stop("Please install the red package")}
    if (is.null(x) == T) {stop("value for parameter x is missing (RasterLayer object)")}
    if (inherits(x, "RasterLayer") == F) {stop("x is not a RasterLayer object")}
#
    cat(paste("Calculation of Area of Occupancy (AOO, km2) and Extent of Occurrence (EOO, km2) (package: red)", "\n\n", sep=""))
#
    freqs <- raster::freq(x)
    freqs <- freqs[complete.cases(freqs),]
    freqs <- freqs[freqs[, 1] > 0, ]
    for (i in (nrow(freqs)-1):1) {freqs[i, 2] <- freqs[i, 2] + freqs[i+1, 2]}
#
    results <- data.frame(array(dim=c(nrow(freqs), 6)))
    names(results) <- c("threshold", "cells", "AOO", "AOO.Type", "EOO", "EOO.Type")
    results[, c(1:2)] <- freqs
#
    reclassify.matrix <- array(dim=c(2, 3))
    reclassify.matrix[1, ] <- c(0.0, 1.0, 0)
    reclassify.matrix[2, ] <- c(0.0, 1.0, 1)
    reclassify.matrix[2, 2] <- results[nrow(results), 1]
#
    for (i in 1:nrow(results)) {
        AOO.Type <- EOO.Type <- c("")
        reclassify.matrix[1, 2] <- reclassify.matrix[2, 1] <- results[i, 1]-0.1
        pres.count <- raster::reclassify(x, reclassify.matrix)
        AOO <- red::aoo(pres.count)
        results[i, "AOO"] <- AOO
        if (AOO < 2000) {AOO.Type <- "possibly Vulnerable (VU)"}
        if (AOO < 500) {AOO.Type <- "possibly Endangered (EN)"}
        if (AOO < 10) {AOO.Type  <- "possibly Critically Endangered (CR)"}
        results[i, "AOO.Type"] <- AOO.Type
        EOO <- red::eoo(pres.count)
        results[i, "EOO"] <- EOO
        if (EOO < 20000) {EOO.Type <- "Vulnerable (VU)"}
        if (EOO < 5000) {EOO.Type <- "possibly Endangered (EN)"}
        if (EOO < 100) {EOO.Type  <- "possibly Critically Endangered (CR)"}
        results[i, "EOO.Type"] <- EOO.Type
    }
    return(results)
}

