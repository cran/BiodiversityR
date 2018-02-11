`ensemble.spatialThin` <- function(
    x, thin.km=0.1, runs=100, verbose=FALSE, return.matrix=FALSE
) 
{
    'distGeo.stack' <- function(x) {
        n <- nrow(x)
        pairs <- utils::combn(n, 2)
        p <- ncol(pairs)
        pairs <- cbind(t(pairs), numeric(p))
        for (i in 1:p) {
            pairs[i, 3] <- geosphere::distGeo(x[pairs[i, 1], ], x[pairs[i, 2], ]) / 1000
        }
        return(pairs)
    }
    'distGeo.matrix.km' <- function(x) {
        dim1 <- nrow(x)
        DistMat <- array(dim=c(dim1, dim1))
        for (i in 1:dim1) {
            for (j in 1:dim1) {
    # distGeo is in m and not km
                DistMat[i, j] <- geosphere::distGeo(x[i, ], x[j,]) / 1000
            }
        }
        return(DistMat)
    }
# distGeo.thin operates on stacked distances (do not recalculate distance each run)
    'distGeo.thin' <- function(x, x2, thin.km=0.1) {
        while(min(x2[, 3]) < thin.km) {
            p <- nrow(x2)
            x2 <- x2[sample(p), ]
            first.min <- which(x2[, 3] < thin.km)
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
        return(x3)
        }
#
    if (raster::couldBeLonLat(sp::SpatialPoints(x)) == F) {
        cat(paste("WARNING: locations not in longitude - latitude format", "\n"))
        cat(paste("therefore spatial thinning not executed", "\n\n"))
        return(x)
    }
    x2 <- distGeo.stack(x)
    if(max(x2[, 3]) <= thin.km) {
        cat(paste("WARNING: thinning parameter larger or equal to maximum distance among locations", "\n"))
        cat(paste("therefore only one location randomly selected", "\n\n"))
        p <- nrow(x)
        x3 <- x[sample(p, 1), ]
        return(x3)
    }
    if(min(x2[, 3]) >= thin.km) {
        cat(paste("WARNING: thinning parameter smaller or equal to minimum distance among locations", "\n"))
        cat(paste("therefore all locations selected", "\n\n"))
        return(x)
    }
    runs <- max(runs, 1)
    locs <- 0
    for (i in 1:runs) {
        loc1 <- distGeo.thin(x, x2, thin.km=thin.km)
        if (verbose == T) {
            if (nrow(loc1) > locs) {cat(paste("run ", i, " (locations: ", nrow(loc1), " > ", locs, " [previous maximum number of locations])", "\n", sep=""))}
            if (nrow(loc1) == locs) {cat(paste("run ", i, " (locations: ", nrow(loc1), " = ", locs, " [previous maximum number of locations])", "\n", sep=""))}
            if (nrow(loc1) < locs) {cat(paste("run ", i, " (locations: ", nrow(loc1), ")", "\n", sep=""))}
        }
        if (nrow(loc1) > locs) { 
            locs <- nrow(loc1)
            loc.out <- loc1
        }
    }
    if (verbose == T) {cat(paste("\n"))}
    loc.matrix <- distGeo.matrix.km(loc.out)
    diag(loc.matrix) <- NA
    cat(paste("Spatially thinned point location data set obtained for target minimum distance of ", thin.km,"\n", sep=""))
    cat(paste("number of locations: ", nrow(loc.out), "\n"))
    cat(paste("minimum distance: ", min(loc.matrix, na.rm=T), "\n"))
    if (return.matrix == T) {return(list(loc.out=loc.out, loc.matrix=distGeo.matrix.km(loc.out)))}
    return(loc.out)
}

