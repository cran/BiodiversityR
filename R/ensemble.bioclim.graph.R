`ensemble.bioclim.graph.data` <- function(
    x=NULL, p=NULL, fraction = 0.9, 
    species.climate.name="Species001_base", factors = NULL
)
{
    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if (is.null(x) == T) {stop("value for parameter xn is missing (RasterStack object)")}
    if(inherits(x, "RasterStack") == F) {stop("x is not a RasterStack object")}
    if(fraction < 0 || fraction > 1) {stop("fraction should be in range 0-1")}

    bioclim.object <- ensemble.bioclim.object(x=x, p=p, fraction=fraction, quantiles=T, species.name=species.climate.name, factors=factors)
    vars <- names(bioclim.object$means)
    nv <- length(vars)
    range.data <- data.frame(array(dim=c(nv, 8)) )
    names(range.data) <- c("species_climate", "biovar", "mean", "median", "min", "max", "lower.limits", "upper.limits")

    range.data[, "species_climate"] <- rep(species.climate.name, nv)
    range.data[, "biovar"] <- names(bioclim.object$means)
    range.data[ , "mean"] <- bioclim.object$means
    range.data[ , "median"] <- bioclim.object$medians
    range.data[ , "min"] <- bioclim.object$minima
    range.data[ , "max"] <- bioclim.object$maxima
    range.data[ , "lower.limits"] <- bioclim.object$lower.limits
    range.data[ , "upper.limits"] <- bioclim.object$upper.limits

    return(range.data)
}

`ensemble.bioclim.graph` <- function(
    graph.data=NULL, focal.var=NULL, species.climates.subset=NULL, cols=NULL,
    var.multiply=1.0, ref.lines=TRUE
)
{
    if (is.null(species.climates.subset) == T) {
        species.climate.subset <- as.character(graph.data[, "species_climate"])
        species.climate.subset <- species.climate.subset[duplicated(species.climate.subset) == F]
    }

    if (is.null(cols) == F) {
        if (length(cols) != length(species.climate.subset)) {stop("different number of colours than number of species and climates to be plotted")}
    }else{
        cols <- grDevices::rainbow(n=length(species.climate.subset), start = 0, end = 5/6)
    }

    graph.data <- graph.data[which(graph.data[, "biovar"] == focal.var), , drop=F]
    graph.data[, "min"] <- graph.data[, "min"] * var.multiply
    graph.data[, "max"] <- graph.data[, "max"] * var.multiply
    graph.data[, "median"] <- graph.data[, "median"] * var.multiply
    graph.data[, "mean"] <- graph.data[, "mean"] * var.multiply
    graph.data[, "lower.limits"] <- graph.data[, "lower.limits"] * var.multiply
    graph.data[, "upper.limits"] <- graph.data[, "upper.limits"] * var.multiply

    nsc <- length(species.climate.subset)
    x1pos <- length(species.climate.subset)+0.5

    graphics::plot(graph.data[, "min"] ~ c(1:nsc), main=focal.var, xlim=c(0.5, x1pos), ylim=c(min(graph.data[, "min"]), max(graph.data[, "max"])), axes=F, xlab="", ylab="", type="n", pch=4, cex=1.2)
    graphics::axis(1, pos=min(graph.data[, "min"]), labels=species.climate.subset, at=c(1:nsc), las=2, cex.axis=1)
    graphics::segments(x0=0.5, y0=min(graph.data[, "min"]), x1=x1pos, y1=min(graph.data[, "min"]))
    graphics::axis(2, pos=0.5, labels=T, las=1, cex.axis=1)
    graphics::segments(x0=0.5, y0=min(graph.data[, "min"]), x1=0.5, y1=max(graph.data[, "max"]))

    if (ref.lines == T) {
        graphics::segments(x0=0.5, y0=graph.data[1, "lower.limits"], x1=x1pos, y1=graph.data[1, "lower.limits"], lty=2, lwd=1, col="grey")
        graphics::segments(x0=0.5, y0=graph.data[1, "upper.limits"], x1=x1pos, y1=graph.data[1, "upper.limits"], lty=2, lwd=1, col="grey")
    }

    # different colours for the species.climate subsets
    for (i in 1:nsc) {
        graphics::points(graph.data[i, "mean"] ~ i, pch=8, cex=2.5, col=cols[i])
        graphics::points(graph.data[i, "median"] ~ i, pch=1, cex=2.5, col=cols[i])
        graphics::points(graph.data[i, "min"] ~ i, pch=1, cex=1.5, col=cols[i])
        graphics::points(graph.data[i, "max"] ~ i, pch=1, cex=1.5, col=cols[i])
        graphics::segments(x0=i, y0=graph.data[i, "lower.limits"], x1=i, y1=graph.data[i, "upper.limits"], lty=1, col=cols[i])
        graphics::segments(x0=i-0.4, y0=graph.data[i, "lower.limits"], x1=i+0.4, y1=graph.data[i, "lower.limits"], lty=1, col=cols[i])
        graphics::segments(x0=i-0.4, y0=graph.data[i, "upper.limits"], x1=i+0.4, y1=graph.data[i, "upper.limits"], lty=1, col=cols[i])
    }
}

