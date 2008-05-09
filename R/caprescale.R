`caprescale` <-
function(x, verbose=FALSE) {
    sumev <- x$tot.chi
    if (is.null(x$CCA)) {
        nr <- nrow(x$CA$u)
        neiv <- length(x$CA$eig)
        sumeiv <- sum(x$CA$eig)
    }
    else {
        nr <- nrow(x$CCA$u)
        nceiv <- length(x$CCA$eig)
        nunceiv <- length(x$CA$eig)
        neiv <- nceiv + nunceiv
        sumeiv <- sum(x$CA$eig) + sum(x$CCA$eig)
    }
    if (substr(x$inertia, 1, 4) == "mean") {
        adjust <- 1
    }
    else {
        adjust <- (nr - 1)
    }
    const <- sqrt(sqrt((nr - 1) * sumev))/adjust
    if (is.null(x$CCA) == F) {
        x$CCA$v <- x$CCA$v * const
        x$CCA$wa <- x$CCA$wa * const
        x$CCA$u <- x$CCA$u * const
        x$CCA$biplot <- x$CCA$biplot * const
        x$CCA$centroids <- x$CCA$centroids * const
    }
    x$CA$v <- x$CA$v * const
    x$CA$u <- x$CA$u * const
    if (substr(x$inertia, 1, 4) == "mean") {
        sstot <- sumeiv * (nr - 1)
    }
    else {
        sstot <- sumeiv/(nr - 1)
    }
    if(verbose==T) {
    cat("SSTot obtained from all eigenvalues:", 
        sstot, "\n")
    distmat <- as.matrix(vegdist(summary(x, axes = neiv, scaling = 1)$sites, 
        method = "euc"))
    sstot <- sum((distmat)^2)/(2 * nrow(distmat))
    cat("SSTot reflected by distances among site scores on all axes:", 
        sstot, "\n")
    if (is.null(x$CCA) == F) {
        if (substr(x$inertia, 1, 4) == "mean") {
            sstot <- sum(x$CCA$eig) * (nr - 1)
        }
        else {
            sstot <- sum(x$CCA$eig)/(nr - 1)
        }
        cat("SSExpl obtained from eigenvalues of constrained axes:", 
            sstot, "\n")
        distmat <- as.matrix(vegdist(summary(x, axes = nceiv, 
            scaling = 1)$constraints, method = "euc"))
        sstot <- sum((distmat)^2)/(2 * nrow(distmat))
        cat("SSExpl reflected by distances among fitted site scores on constrained axes (scaling 1):", 
            sstot, "\n")
        if (substr(x$inertia, 1, 4) == "mean") {
            sstot <- sum(x$CA$eig) * (nr - 1)
        }
        else {
            sstot <- sum(x$CA$eig)/(nr - 1)
        }
        cat("SSRes obtained from eigenvalues of unconstrained axes:", 
            sstot, "\n")
        distmat <- as.matrix(vegdist(summary(x, axes = neiv, 
            scaling = 1)$sites[, ((nceiv + 1):neiv)], method = "euc"))
        sstot <- sum((distmat)^2)/(2 * nrow(distmat))
        cat("SSRes reflected by distances among site scores on unconstrained axes (scaling 1):", 
            sstot, "\n")
    }
    }
    return(x)
}

