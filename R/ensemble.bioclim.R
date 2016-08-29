`ensemble.bioclim.object` <- function(
    x=NULL, p=NULL, fraction=0.9,
    quantiles=FALSE, 
    species.name="Species001", 
    factors=NULL, dummy.vars=NULL
)
{
    if(is.null(x) == T) {stop("value for parameter x is missing (data.frame or RasterStack object)")}
    if(inherits(x, "RasterStack")==F && inherits(x, "data.frame")==F) {stop("x should be a data.frame or RasterStack object")}
    if(fraction < 0 || fraction > 1) {stop("fraction should be in range 0-1")}

    factors <- as.character(factors)
    dummy.vars <- as.character(dummy.vars)

# cutoff parameter based on the normal distribution
    cutoff <- qnorm(0.5+fraction/2)   
    probs <- c(0.5-fraction/2, 0.5+fraction/2)

    if(inherits(x, "RasterStack")==T  && is.null(p)==F) {
        clim.values <- data.frame(raster::extract(x, y=p))
        names(clim.values) <- names(x)
        x <- clim.values        
    }

    if(inherits(x, "data.frame") == F) {
        vars <- names(x)
        if (length(factors) > 0) {for (i in 1:length(factors)) {vars <- vars[which(names(x) != factors[i])]}}
        if (length(dummy.vars) > 0) {for (i in 1:length(dummy.vars)) {vars <- vars[which(names(x) != dummy.vars[i])]}}
        nv <- length(vars)
        lower.limitsq <- upper.limitsq <- lower.limits <- upper.limits <- minima <- maxima <- clim.sd <- clim.mean <- numeric(length=nv)
        names(lower.limitsq) <- names(upper.limitsq) <- names(lower.limits) <- names(upper.limits) <- names(minima) <- names(maxima) <- names(clim.sd) <- names(clim.mean) <- vars
        for (i in 1:nv) {
            vari <- vars[which(vars == names(x)[i])]
            raster.focus <- x[[i]]
            raster::setMinMax(raster.focus)
            meanV <- raster::cellStats(raster.focus, 'mean')
            sdV <- raster::cellStats(raster.focus, 'sd')
            minV <- raster::minValue(raster.focus)
            maxV <- raster::maxValue(raster.focus)
            lowerV <- as.numeric(quantile(raster.focus, probs[1]))
            upperV <- as.numeric(quantile(raster.focus, probs[2]))

            lower.limitsq[which(names(lower.limitsq) == vari)] <- lowerV
            upper.limitsq[which(names(upper.limitsq) == vari)] <- upperV
            clim.mean[which(names(clim.mean) == vari)] <- meanV
            clim.sd[which(names(clim.sd) == vari)] <- sdV
            minima[which(names(minima) == vari)] <- minV
            maxima[which(names(maxima) == vari)] <- maxV
        }

    }else{
        clim.values <- x
        for (i in 1:length(names(clim.values))) {if (is.factor(clim.values[, i]) == T) {factors <- c(factors, names(clim.values)[i])} }        
        factors <- unique(factors)
        if (length(factors) > 0) {for (i in 1:length(factors)) {clim.values <- clim.values[, which(names(clim.values) != factors[i])]}}
        if (length(dummy.vars) > 0) {for (i in 1:length(dummy.vars)) {clim.values <- clim.values[, which(names(clim.values) != dummy.vars[i])]}}

        clim.mean <- apply(clim.values, 2, "mean", na.rm=T)
        clim.sd <- apply(clim.values, 2, "sd", na.rm=T)

        lower.limitsq <- upper.limitsq <- lower.limits <- upper.limits <- minima <- maxima <- numeric(length=length(clim.mean))
        names(lower.limitsq) <- names(upper.limitsq) <- names(lower.limits) <- names(upper.limits) <- names(minima) <- names(maxima) <- names(clim.values)
        minima <- apply(clim.values, 2, "min", na.rm=T)
        maxima <- apply(clim.values, 2, "max", na.rm=T)          
        lower.limitsq <- apply(clim.values, 2, "quantile", probs[1], na.rm=T)
        upper.limitsq <- apply(clim.values, 2, "quantile", probs[2], na.rm=T)
    }

    if (quantiles == F){  
        lower.limits <- clim.mean - cutoff*clim.sd
        upper.limits <- clim.mean + cutoff*clim.sd
    }else{
        lower.limits <- lower.limitsq
        upper.limits <- upper.limitsq
    }

# deal with asymmetrical distributions
    for (i in 1:length(lower.limits)) {
        if (lower.limits[i] < minima[i]) {
            cat(paste("\n", "WARNING: lower limit of ", lower.limits[i], " for ", names(lower.limits)[i], " was smaller than minimum of ", minima[i], sep = ""))
            cat(paste("\n", "lower limit therefore replaced by quantile value of ", lower.limitsq[i], "\n", sep = ""))
            lower.limits[i] <- lower.limitsq[i]
        }
        if (upper.limits[i] > maxima[i]) {
            cat(paste("\n", "WARNING: upper limit of ", upper.limits[i], " for ", names(upper.limits)[i], " was larger than maximum of ", maxima[i], sep = ""))
            cat(paste("\n", "upper limit therefore replaced by quantile value of ", upper.limitsq[i], "\n", sep = ""))
            upper.limits[i] <- upper.limitsq[i]
        }
    }

    return(list(lower.limits=lower.limits, upper.limits=upper.limits, minima=minima, maxima=maxima, 
        means=clim.mean, sds=clim.sd, cutoff=cutoff, fraction=fraction, species.name=species.name))
}


`ensemble.bioclim` <- function(
    x=NULL, bioclim.object=NULL, 
    RASTER.object.name=bioclim.object$species.name, RASTER.stack.name = x@title,
    RASTER.format="raster",
    KML.out=TRUE, KML.blur=10, KML.maxpixels=100000 
)
{
    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if(is.null(x) == T) {stop("value for parameter x is missing (RasterStack object)")}
    if(inherits(x, "RasterStack") == F) {stop("x is not a RasterStack object")}
    if (is.null(bioclim.object) == T) {stop("value for parameter bioclim.object is missing (hint: use the ensemble.bioclim.object function)")}
# 
    predict.bioclim <- function(object=bioclim.object, newdata=newdata) {
        lower.limits <- object$lower.limits
        upper.limits <- object$upper.limits
        minima <- object$minima
        maxima <- object$maxima

        newdata <- newdata[, which(names(newdata) %in% names(lower.limits)), drop=F]
        result <- as.numeric(rep(NA, nrow(newdata)))
        varnames <- names(newdata)
        nvars <- ncol(newdata)

        for (i in 1:nrow(newdata)) {
            datai <- newdata[i,,drop=F]
            resulti <- 1
            j <- 0
            while (resulti > 0 && j <= (nvars-1)) {
                j <- j+1
                focal.var <- varnames[j]
                if (resulti == 1) {
                    lowerj <- lower.limits[which(names(lower.limits) == focal.var)]
                    if (datai[, j]  < lowerj) {resulti <- 0.5}
                    upperj <- upper.limits[which(names(upper.limits) == focal.var)]
                    if (datai[, j]  > upperj) {resulti <- 0.5}
                }
                minj <- minima[which(names(minima) == focal.var)]
                if (datai[, j]  < minj) {resulti <- 0}
                maxj <- maxima[which(names(maxima) == focal.var)]
                if (datai[, j]  > maxj) {resulti <- 0}
            }
            result[i] <- resulti
        }

        p <- as.numeric(result)
        return(p)
    }
  
# avoid problems with non-existing directories and prepare for output
    dir.create("ensembles", showWarnings = F)
    if (KML.out == T) {dir.create("kml", showWarnings = F)}
    if(length(x@title) == 0) {x@title <- "stack1"}
    stack.title <- RASTER.stack.name
    rasterfull <- paste("ensembles//", RASTER.object.name, "_", stack.title , "_BIOCLIM_orig", sep="")
    kmlfull <- paste("kml//", RASTER.object.name, "_", stack.title , "_BIOCLIM_orig", sep="")
  
#
# predict
    tryCatch(bioclim.raster <- raster::predict(object=x, model=bioclim.object, fun=predict.bioclim, na.rm=TRUE, 
                                           filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format),
           error= function(err) {print(paste("prediction of bioclim failed"))},
           silent=F)
#    bioclim.raster <- trunc(1000*bioclim.raster)
#    cat(paste("\n", "raster layer created (probabilities multiplied by 1000)", "\n", sep = ""))
    raster::setMinMax(bioclim.raster)
    print(bioclim.raster)

#
# avoid possible problems with saving of names of the raster layers
    raster::writeRaster(bioclim.raster, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- paste(RASTER.object.name, "_", stack.title , "_BIOCLIM_orig", sep="")
    raster::writeRaster(working.raster, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format)
#  
    if (KML.out == T) {
#        working.raster <- trunc(1000*working.raster)
        raster::KML(working.raster, filename=kmlfull, col = c("grey", "orange", "green"), colNA = 0, 
            blur=KML.blur, maxpixels=KML.maxpixels, overwrite=T, breaks = c(-0.1, 0, 0.5, 1.0))
    }
  
    cat(paste("\n", "bioclim raster provided in folder: ", getwd(), "//ensembles", "\n", sep=""))
    return(bioclim.raster)
}


