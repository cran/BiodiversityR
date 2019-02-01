`ensemble.outlierThin` <- function(
    x, predictors.stack=NULL, k=10, quant=0.95, pca.var=0.95,
    return.outliers=FALSE
) 
{

    .BiodiversityR <- new.env()
#
# create background data
    background.data <- raster::extract(predictors.stack, x)
    background.data <- data.frame(background.data)
    TrainValid <- complete.cases(background.data)
    x <- x[TrainValid,]
    background.data <- background.data[TrainValid,]

# PCA of scaled variables
    rda.result <- vegan::rda(X=background.data, scale=T)
# select number of axes
    ax <- 2
    while ( (sum(vegan::eigenvals(rda.result)[c(1:ax)])/sum(vegan::eigenvals(rda.result))) < pca.var ) {ax <- ax+1}
    cat(paste("\n", "Percentage of variance of the selected axes (1 to ", ax, ") of principal components analysis: ", 100*sum(vegan::eigenvals(rda.result)[c(1:ax)])/sum(vegan::eigenvals(rda.result)), "\n", sep = ""))
    rda.scores <- vegan::scores(rda.result, display="sites", scaling=1, choices=c(1:ax))
#
    lof.result <- Rlof::lof(rda.scores, k=10, method="euclidean")
    outliers.limit <- quantile(lof.result, probs=quant)
#
    inliers <- x[lof.result < outliers.limit, ]
    outliers <- x[lof.result >= outliers.limit, ]
    cat(paste(quant, " quantile limit for local outliers: ", outliers.limit, "\n", sep=""))
#
    if (return.outliers == T) {
        return(list(inliers=inliers, outliers=outliers))      
    }else{
        return(inliers)
    }
}

