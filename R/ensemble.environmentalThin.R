`ensemble.environmentalThin` <- function(
    x, predictors.stack=NULL, thin.n=50, runs=100, pca.var=0.95, verbose=FALSE,
    return.notRetained=FALSE
) 
{

    .BiodiversityR <- new.env()

# distEnviro.thin operates on stacked distances (do not recalculate distance each run)
    'distEnviro.thin' <- function(x, x2, thin.dist=0.1) {
        while(min(x2, 3) < thin.dist) {
            p <- nrow(x2)
            x2 <- x2[sample(p), ]
            first.min <- which(x2[, 3] < thin.dist)
            first.min <- first.min[1]
            random.col <- as.numeric(runif(1) > 0.5)+1
            selected <- x2[first.min, random.col]
            rows1 <- x2[, 1] != selected
            x2 <- x2[rows1, , drop=F]
            rows2 <- x2[, 2] != selected
            x2 <- x2[rows2, , drop=F]
        }
        retained <- unique(c(x2[, 1], x2[, 2]))
        x3 <- x[retained, ]
        dist.min <- min(x2[, 3])
        return(list(x3=x3, dist.min=dist.min, retained=retained))
        }
#
    if(thin.n >= nrow(x)) {
        cat(paste("WARNING: thinning parameter larger or equal to number of available locations", "\n"))
        cat(paste("therefore all locations selected", "\n\n"))
        return(list(retained=x, not.retained=NULL))
    }
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
    rda.dist <- as.matrix(vegdist(rda.scores, method="euc"))
#
# stack
    n <- nrow(x)
    pairs <- utils::combn(n, 2)
    p <- ncol(pairs)
    pairs <- cbind(t(pairs), numeric(p))
    for (i in 1:p) {
        pairs[i, 3] <- rda.dist[pairs[i, 1], pairs[i, 2]]
    }
#
    runs <- max(runs, 1)
    dist.all <- 0
#
# algorithm 1 removes from pairs with smallest distance
    x2 <- pairs
    retained.n <- length(unique(c(x2[, 1], x2[, 2])))
    while(retained.n > thin.n) {
        first.min <- which.min(x2[, 3])
        first.min <- first.min[1]
        random.col <- as.numeric(runif(1) > 0.5)+1
        selected <- x2[first.min, random.col]
        rows1 <- x2[, 1] != selected
        x2 <- x2[rows1, , drop=F]
        rows2 <- x2[, 2] != selected
        x2 <- x2[rows2, , drop=F]
        retained <- unique(c(x2[, 1], x2[, 2]))
        retained.n <- length(retained)
    }
    x3 <- x[retained, ]
    dist.min1 <- min(x2[, 3])
    cat(paste("Environmentally thinned point location data set obtained with first algorithm", "\n", sep=""))
    cat(paste("number of locations: ", nrow(x3), "\n"))
    cat(paste("minimum distance: ", dist.min1, "\n"))
#
# algorithm 2 uses minimum distance of previous algorithm
# now algorithm attempts to maximize the number of retained locations, similar to ensemble.spatialThin
#
    cat(paste("\n", "Environmentally thinned point location data set obtained with second algorithm", "\n", sep=""))
    dist.all <- 0
    dist.n.all <- 0
    for (i in 1:runs) {
        dist1 <- distEnviro.thin(x, x2=pairs, thin.dist=dist.min1)
        dist.min2 <- dist1$dist.min
        dist.n <- nrow(dist1$x3)
        if (verbose == T) {
            if (dist.min2 > dist.all) {cat(paste("run ", i, " (", dist.n, " locations with minimum distance: ", dist.min2, " > ", dist.all, " [previous minimum distance])", "\n", sep=""))}
            if (dist.min2 == dist.all) {cat(paste("run ", i, " (", dist.n, " locations with minimum distance: ", dist.min2, " = ", dist.all, " [previous minimum distance])", "\n", sep=""))}
            if (dist.min2 < dist.all) {cat(paste("run ", i, " (", dist.n, " locations with minimum distance: ", dist.min2, ")", "\n", sep=""))}
        }
        if (dist.min2 > dist.all) { 
            dist.all <- dist.min2
            loc.out <- dist1$x3
            dist.n.all <- dist.n
            retained <- dist1$retained
        }
        if (dist.min2 == dist.all && dist.n > dist.n.all) { 
            dist.all <- dist.min2
            loc.out <- dist1$x3
            dist.n.all <- dist.n
            retained <- dist1$retained
        }
    }
    if (verbose == T) {cat(paste("\n"))}
    cat(paste("number of locations: ", nrow(loc.out), "\n"))
    cat(paste("minimum distance: ", dist.all, "\n"))
    if (return.notRetained == T) {
        x.not <- x[(c(1:nrow(x)) %in% retained) == F, ]
        return(list(retained=loc.out, not.retained=x.not))
    }else{
        return(loc.out)
    }
}

