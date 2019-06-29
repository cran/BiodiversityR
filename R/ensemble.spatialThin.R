`ensemble.spatialThin` <- function(
    x, thin.km=0.1, runs=100, silent=FALSE, verbose=FALSE, 
    return.notRetained=FALSE
)
{
    'distGeo.stack' <- function(x) {
        n <- nrow(x)
        pairs <- utils::combn(n, 2)
        p <- ncol(pairs)
        pairs <- cbind(t(pairs), numeric(p))
        for (i in 1:p) {
            pairs[i, 3] <- geosphere::distGeo(x[pairs[i, 1], ], x[pairs[i, 2], ]) / 1000
            pairs[i, 3] <- round(pairs[i, 3], 2)
        }
        return(pairs)
    }

# distGeo.thin operates on stacked distances (do not recalculate distance each run)
    'distGeo.thin' <- function(x, x2, thin.km=0.1) {
        while(min(x2[, 3]) < thin.km  && nrow(x2) > 1) {
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
        x3 <- x[retained, , drop=F]
# special case where the remaining 2 locations are closer than minimum distance
        if (nrow(x3)==2 && distGeo.stack(x3)[3] < thin.km) {
            retained <- sample(retained, size=1)
            x3 <- x[retained, , drop=F]
        }
        return(list(x3=x3, retained=retained))
        }
#
    if (verbose == T) {silent <- FALSE}
    if (raster::couldBeLonLat(sp::SpatialPoints(x)) == F) {
        cat(paste("WARNING: locations not in longitude - latitude format", "\n"))
        cat(paste("therefore spatial thinning not executed", "\n\n"))
        return(x)
    }
    if(nrow(x) == 1) {
        if (silent == F) {cat(paste("NOTE: only one location was provided", "\n"))}
        return(x)
    }
    x2 <- distGeo.stack(x)
    if(max(x2[, 3]) <= thin.km) {
        if (silent == F) {
            cat(paste("WARNING: thinning parameter larger or equal to maximum distance among locations", "\n"))
            cat(paste("therefore only one location randomly selected", "\n\n"))
        }
        p <- nrow(x)
        x3 <- x[sample(p, 1), ]
        return(x3)
    }
    if(min(x2[, 3]) >= thin.km) {
        if (silent == F) {
            cat(paste("WARNING: thinning parameter smaller or equal to minimum distance among locations", "\n"))
            cat(paste("therefore all locations selected", "\n\n"))
        }
        return(x)
    }
    runs <- max(runs, 1)
    locs <- 0
    for (i in 1:runs) {
        loc.l1 <- distGeo.thin(x, x2, thin.km=thin.km)
        loc1 <- loc.l1$x3
        if (verbose == T) {
            if (nrow(loc1) > locs) {cat(paste("run ", i, " (locations: ", nrow(loc1), " > ", locs, " [previous maximum number of locations])", "\n", sep=""))}
            if (nrow(loc1) == locs) {cat(paste("run ", i, " (locations: ", nrow(loc1), " = ", locs, " [previous maximum number of locations])", "\n", sep=""))}
            if (nrow(loc1) < locs) {cat(paste("run ", i, " (locations: ", nrow(loc1), ")", "\n", sep=""))}
        }
        if (nrow(loc1) > locs) { 
            locs <- nrow(loc1)
            loc.out <- loc1
            retained.final <- loc.l1$retained
        }
    }
    if (verbose == T) {cat(paste("\n"))}
    loc.matrix <- geosphere::distm(loc.out)
    diag(loc.matrix) <- NA
    if (silent == F) {
        cat(paste("Spatially thinned point location data set obtained for target minimum distance of ", thin.km, "\n", sep=""))
        cat(paste("number of locations: ", nrow(loc.out), "\n"))
        cat(paste("minimum distance: ", min(loc.matrix, na.rm=T), "\n"))
    }
    if (return.notRetained == T) {
        x.not <- x[(c(1:nrow(x)) %in% retained.final) == F, ]
        return(list(loc.out=loc.out, not.retained=x.not))
    }
    if (return.notRetained == F) {return(loc.out)}
}


`ensemble.spatialThin.quant` <- function(
    x, thin.km=0.1, runs=100, silent=FALSE, verbose=FALSE, 
    LON.length=21, LAT.length=21
) 
{
    LON.bins <- quantile(x[, 1], probs=seq(from=0, to=1, length=LON.length))
    LAT.bins <- quantile(x[, 2], probs=seq(from=0, to=1, length=LAT.length))
#
    LON.bins[length(LON.bins)] <- +Inf
    LAT.bins[length(LAT.bins)] <- +Inf
#
    loc.first <- TRUE
    for (lo in 1:(length(LON.bins)-1)) {
        x.bin <- x[(x[, 1] >= LON.bins[lo] & x[, 1] < LON.bins[lo+1]),  ]
        for (la in 1:(length(LAT.bins)-1)) {
            x.bin2 <- x.bin[(x.bin[, 2] >= LAT.bins[la] & x.bin[, 2] < LAT.bins[la+1]),  ]     
            if (nrow(as.data.frame(x.bin2)) > 0) {
                x.spat <- ensemble.spatialThin(x=x.bin2, thin.km=thin.km, runs=runs, verbose=verbose, silent=silent)
                if (loc.first == T) {
                    x.out <- x.spat
                    loc.first <- FALSE
                }else{
                    x.out <- rbind(x.out, x.spat)
                }
            }
        }
    }
    return(x.out)
}
