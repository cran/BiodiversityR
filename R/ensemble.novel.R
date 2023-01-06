`ensemble.novel.object` <- function(
    x=NULL, name="reference1", mask.raster=NULL, 
    quantiles=FALSE, probs=c(0.05, 0.95),
    factors=NULL
)
{
    vars <- names(x)
    if (length(factors) > 0) {for (i in 1:length(factors)) {vars <- vars[which(names(x) != factors[i])]}}
    nv <- length(vars)
    minima <- maxima <- numeric(nv)
    names(minima) <- names(maxima) <- vars
    if(inherits(x, "RasterStack") == T) {
        for (i in c(1:nv)) {
            vari <- vars[which(names(x) == vars[i])]
            raster.focus <- x[[which(names(x) == vars[i])]]
            if (is.null(mask.raster) == F) {raster.focus <- raster::mask(x[[which(names(x) == vars[i])]], mask=mask.raster)}
            raster::setMinMax(raster.focus)
            print(raster.focus)
            if (quantiles == F) {
                minV <- raster::minValue(raster.focus)
                maxV <- raster::maxValue(raster.focus)
            }else{
                minV <- as.numeric(raster::quantile(raster.focus, probs=probs[1], na.rm=T))
                maxV <- as.numeric(raster::quantile(raster.focus, probs=probs[2], na.rm=T))
            }
            minima[which(names(minima) == vari)] <- minV
            maxima[which(names(maxima) == vari)] <- maxV
        }
    }
    if(inherits(x, "data.frame") == T) {
        for (i in 1:nv) {
            vari <- vars[i]
            xdata <- x[, which(names(x)==vari), drop=F]
            if (quantiles == F) {
                minV <- min(xdata)
                maxV <- max(xdata)
            }else{
                minV <- as.numeric(stats::quantile(xdata, probs[1], na.rm=T))
                maxV <- as.numeric(stats::quantile(xdata, probs[2], na.rm=T))
            }
            minima[which(names(minima) == vari)] <- minV
            maxima[which(names(maxima) == vari)] <- maxV
        }
    }
    novel.object <- list(minima=minima, maxima=maxima, name=name)
    return(novel.object)
}

`ensemble.novel` <- function(
  x=NULL, novel.object=NULL, 
  RASTER.object.name=novel.object$name, RASTER.stack.name = x@title,
  RASTER.format="GTiff", RASTER.datatype="INT2S", RASTER.NAflag=-32767,
#  KML.out=FALSE, KML.maxpixels=100000, KML.blur=10,
  CATCH.OFF=FALSE
)
{
  .BiodiversityR <- new.env()
  #    if (! require(dismo)) {stop("Please install the dismo package")}
  if(is.null(x) == T) {stop("value for parameter x is missing (RasterStack object)")}
  if(inherits(x, "RasterStack") == F) {stop("x is not a RasterStack object")}
  if (is.null(novel.object) == T) {stop("value for parameter novel.object is missing (hint: use the ensemble.novel.object function)")}
  if (all.equal(names(novel.object$minima), names(novel.object$maxima)) == F) {{stop("different variable names for maxima and minima")}}
  # 
# 
#    if (KML.out==T && raster::isLonLat(x)==F) {
#        cat(paste("\n", "NOTE: not possible to generate KML files as Coordinate Reference System (CRS) of stack ", x@title , " is not longitude and latitude", "\n", sep = ""))
#        KML.out <- FALSE
#    }
#
  predict.novel <- function(object=novel.object, newdata=newdata) {
    minima <- object$minima
    maxima <- object$maxima
    
    newdata <- newdata[, which(names(newdata) %in% names(minima)), drop=F]
    result <- as.numeric(rep(NA, nrow(newdata)))
    varnames <- names(newdata)
    nvars <- ncol(newdata)
    
    #        for (i in 1:nrow(newdata)) {
    #            datai <- newdata[i,,drop=F]
    #            if (any(is.na(datai)) == F) {
    #                datai1 <- rbind(minima, datai)
    #                mins <- sum(as.numeric(apply(datai1, 2, which.min)))
    #                datai2 <- rbind(maxima, datai)
    #                maxs <- sum(as.numeric(apply(datai2, 2, which.max)))
    #                if ((mins+maxs) == 2*nv) {
    #                    result[i] <- 0
    #                }else{
    #                    result[i] <- 1
    #                }
    #            }
    #        }
    
    for (i in 1:nrow(newdata)) {
      datai <- newdata[i,,drop=F]
      resulti <- 0
      j <- 0
      while (resulti < 1 && j <= (nvars-1)) {
        j <- j+1
        focal.var <- varnames[j]
        minj <- minima[which(names(minima) == focal.var)]
        if (datai[, j]  < minj) {resulti <- 1}
        maxj <- maxima[which(names(maxima) == focal.var)]
        if (datai[, j]  > maxj) {resulti <- 1}
      }
      result[i] <- resulti
    }
    
    p <- as.numeric(result)
    p <- trunc(p)
    return(p)
  } 
  
  # check if all variables are present
  vars <- names(novel.object$minima)
  vars.x <- names(x)
  nv <- length(vars) 
  for (i in 1:nv) {
    if (any(vars.x==vars[i]) == F) {stop("explanatory variable '", vars[i], "' not among grid layers of RasterStack x \n", sep = "")}
  }
  nv <- length(vars.x) 
  for (i in 1:nv) {
    if (any(vars==vars.x[i]) == F) {
      cat(paste("\n", "NOTE: RasterStack layer '", vars.x[i], "' was not documented in the novel object data set and will be ignored", "\n", sep = ""))
      x <- raster::dropLayer(x, which(names(x) %in% c(vars.x[i]) ))
      x <- raster::stack(x)
    }
  }
  
  # same order of variables in stack as in novel object
  minima <- novel.object$minima
  minima <- minima[as.numeric(na.omit(match(names(x), names(minima))))]
  novel.object$minima <- minima
  maxima <- novel.object$maxima
  maxima <- maxima[as.numeric(na.omit(match(names(x), names(maxima))))]
  novel.object$maxima <- maxima
  
  # avoid problems with non-existing directories and prepare for output
  dir.create("ensembles", showWarnings = F)
  dir.create("ensembles/novel", showWarnings = F)
#  if (KML.out == T) {
#    dir.create("kml", showWarnings = F)
#    dir.create("kml/novel", showWarnings = F)
#  }
  if(length(x@title) == 0) {x@title <- "stack1"}
  stack.title <- RASTER.stack.name
  if (gsub(".", "_", stack.title, fixed=T) != stack.title) {cat(paste("\n", "WARNING: title of stack (", stack.title, ") contains '.'", "\n\n", sep = ""))}
  rasterfull <- paste("ensembles/novel/", RASTER.object.name, "_", stack.title , "_novel", sep="")
#  kmlfull <- paste("kml/novel/", RASTER.object.name, "_", stack.title , "_novel", sep="")
  
  #
  # predict
  if (CATCH.OFF == F) {  
      tryCatch(novel.raster <- raster::predict(object=x, model=novel.object, fun=predict.novel, na.rm=TRUE, 
               filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format),
           error= function(err) {print(paste("prediction of novel zones failed"))},
           silent=F)
  }else{
      novel.raster <- raster::predict(object=x, model=novel.object, fun=predict.novel, na.rm=TRUE, 
           filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format)
  }
  # save as integer
#  novel.raster <- trunc(novel.raster)
  raster::setMinMax(novel.raster)
  cat(paste("\n", "raster layer with novel areas (1) created (0 indicates not novel)", "\n", sep = ""))
  print(novel.raster)
  print(raster::freq(novel.raster))
  
  #
  # avoid possible problems with saving of names of the raster layers
  # no longer used with default GTiff format since DEC-2022
#  raster::writeRaster(novel.raster, filename="working.grd", overwrite=T)
#  working.raster <- raster::raster("working.grd")
#  names(working.raster) <- paste(RASTER.object.name, "_", stack.title , "_novel", sep="")
#  raster::writeRaster(working.raster, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
  
  #
#  if (KML.out == T) {
#    raster::KML(working.raster, filename=kmlfull, col = c("grey", "green"), colNA = 0, 
#                blur=KML.blur, maxpixels=KML.maxpixels, overwrite=TRUE, breaks = c(-1, 0, 1))
#  }
  
  cat(paste("\n", "novel climate raster provided in folder: ", getwd(), "//ensembles//novel", "\n", sep=""))
#  novel.raster <- raster::raster(rasterfull)
  return(novel.raster)
}

