`ensemble.environmentalThin` <- function(
    x, predictors.stack=NULL, thin.n=50, runs=100, pca.var=0.95, 
    silent=FALSE, verbose=FALSE,
    return.notRetained=FALSE
) 
{

    .BiodiversityR <- new.env()

# distEnviro.thin operates on stacked distances (do not recalculate distance each run)
    'distEnviro.thin' <- function(x, x2, thin.dist=0.1, thin.n=50) {
        while(min(x2[, 3]) < thin.dist && nrow(x2) > 1) {
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
# special case where the remaining 2 locations are closer than minimum distance
        if (nrow(x3)==2 && x2[1, 3] < thin.dist) {
            retained <- sample(retained, size=1)
            x3 <- x[retained, , drop=F]
        }
# rather than maximizing number of locations, keep increasing distance among them until target thin.n
# now again use the first algorithm
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
        dist.min <- min(x2[, 3])
        return(list(x3=x3, dist.min=dist.min, retained=retained))
        }
#
    if (verbose == T) {silent <- FALSE}
#
    if(thin.n >= nrow(x)) {
        if (silent == F) {
            cat(paste("WARNING: thinning parameter larger or equal to number of available locations", "\n"))
            cat(paste("therefore all locations selected", "\n\n"))
        }
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
    ax <- 1
    while ( (sum(vegan::eigenvals(rda.result)[c(1:ax)])/sum(vegan::eigenvals(rda.result))) < pca.var ) {ax <- ax+1}
    if (silent == F) {cat(paste("\n", "Percentage of variance of the selected axes (1 to ", ax, ") of principal components analysis: ", 100*sum(vegan::eigenvals(rda.result)[c(1:ax)])/sum(vegan::eigenvals(rda.result)), "\n", sep = ""))}
    rda.scores <- vegan::scores(rda.result, display="sites", scaling=1, choices=c(1:ax))
    rda.dist <- as.matrix(vegdist(rda.scores, method="euc"))
    rda.dist <- signif(rda.dist, digits=6)
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
    if (silent == F) {
        cat(paste("Environmentally thinned point location data set obtained with first algorithm", "\n", sep=""))
        cat(paste("number of locations: ", nrow(x3), "\n"))
        cat(paste("minimum distance: ", dist.min1, "\n"))
    }
#
# algorithm 2 uses minimum distance of previous algorithm
# now algorithm attempts to maximize the number of retained locations, similar to ensemble.spatialThin
#
    if (silent == F) {cat(paste("\n", "Environmentally thinned point location data set obtained with second algorithm", "\n", sep=""))}
    dist.all <- 0
    dist.n.all <- 0
    for (i in 1:runs) {
        dist1 <- distEnviro.thin(x, x2=pairs, thin.dist=dist.min1, thin.n=thin.n)
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
    if (silent == F) {
        cat(paste("number of locations: ", nrow(loc.out), "\n"))
        cat(paste("minimum distance: ", dist.all, "\n"))
    }
    if (return.notRetained == T) {
        x.not <- x[(c(1:nrow(x)) %in% retained) == F, ]
        return(list(retained=loc.out, not.retained=x.not))
    }else{
        return(loc.out)
    }
}


`ensemble.environmentalThin.clara` <- function(
    x, predictors.stack=NULL, thin.n=20, runs=100, pca.var=0.95, 
    silent=FALSE, verbose=FALSE,
    clara.k=100
) 
{
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
    ax <- 1
    while ( (sum(vegan::eigenvals(rda.result)[c(1:ax)])/sum(vegan::eigenvals(rda.result))) < pca.var ) {ax <- ax+1}
    if (silent == F) {cat(paste("\n", "Percentage of variance of the selected axes (1 to ", ax, ") of principal components analysis: ", 100*sum(vegan::eigenvals(rda.result)[c(1:ax)])/sum(vegan::eigenvals(rda.result)), "\n", sep = ""))}
    rda.scores <- vegan::scores(rda.result, display="sites", scaling=1, choices=c(1:ax))
#
# cluster and thin if more locations in each cluster
    clara.result <- cluster::clara(rda.scores, k=clara.k, metric="euclidean", medoids.x=F)$clustering
#
    loc.first <- TRUE
    for (lo in 1:clara.k) {
        x.bin <- x[clara.result == lo,  , drop=F]
        x.env <- x.bin
        if (nrow(as.data.frame(x.bin)) > thin.n) {
            x.env <- ensemble.environmentalThin(x=x.bin, predictors.stack=predictors.stack, thin.n=thin.n, runs=runs, pca.var=pca.var, verbose=verbose, silent=silent)
        }
        if (loc.first == T) {
            x.out <- x.env
            loc.first <- FALSE
        }else{
            x.out <- rbind(x.out, x.env)
        }
    }
    return(x.out)
}

