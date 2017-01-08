`ensemble.raster` <- function(
    xn=NULL, 
    models.list=NULL, 
    input.weights=models.list$output.weights,
    thresholds=models.list$thresholds,
    RASTER.species.name=models.list$species.name, 
    RASTER.stack.name=xn@title, 
    RASTER.format="raster", RASTER.datatype="INT2S", RASTER.NAflag=-32767,
    RASTER.models.overwrite=TRUE,
    KML.out=FALSE, KML.maxpixels=100000, KML.blur=10,
    evaluate=FALSE, SINK=FALSE,
    p=models.list$p, a=models.list$a,
    pt=models.list$pt, at=models.list$at
)
{
    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if (is.null(xn) == T) {stop("value for parameter xn is missing (RasterStack object)")}
    if(inherits(xn, "RasterStack") == F) {stop("xn is not a RasterStack object")}
    if (is.null(models.list) == T) {stop("provide 'models.list' as models will not be recalibrated and retested")}
    if (is.null(input.weights) == T) {input.weights <- models.list$output.weights}
    if (is.null(thresholds) == T) {stop("provide 'thresholds' as models will not be recalibrated and retested")}
    thresholds.raster <- thresholds
    if(is.null(p) == F) {names(p) <- c("x", "y")}
    if(is.null(a) == F) {names(a) <- c("x", "y")}
    if(is.null(pt) == F) {names(pt) <- c("x", "y")}
    if(is.null(at) == F) {names(at) <- c("x", "y")}
#
    if (KML.out==T && raster::isLonLat(xn)==F) {
        cat(paste("\n", "NOTE: not possible to generate KML files as Coordinate Reference System (CRS) of stack ", xn@title , " is not longitude and latitude", "\n", sep = ""))
        KML.out <- FALSE
    }
#
    retest <- FALSE
    if (evaluate == T) { 
        if (is.null(p)==T || is.null(a)==T) {
            cat(paste("\n", "NOTE: not possible to evaluate the models since locations p and a are not provided", "\n", sep = ""))
            evaluate <- FALSE
        }else{
            threshold.method <- models.list$threshold.method
            threshold.sensitivity <- models.list$threshold.sensitivity
            threshold.PresenceAbsence <- models.list$threshold.PresenceAbsence <- FALSE
        }
        if (is.null(pt)==F && is.null(at)==F) {
            if(identical(pt, p) == F || identical(at, a) == F)  {retest <- TRUE}
        }
    }
#
# create output file
    dir.create("outputs", showWarnings = F)
    paste.file <- paste(getwd(), "/outputs/", RASTER.species.name, "_output.txt", sep="")
    OLD.SINK <- TRUE
    if (sink.number(type="output") == 0) {OLD.SINK <- F}
    if (SINK==T && OLD.SINK==F) {
        if (file.exists(paste.file) == F) {
            cat(paste("\n", "NOTE: results captured in file: ", paste.file, "\n", sep = ""))
        }else{
            cat(paste("\n", "NOTE: results appended in file: ", paste.file, "\n", sep = ""))
        }
        cat(paste("\n\n", "RESULTS (ensemble.raster function)", "\n\n", sep=""), file=paste.file, append=T)
        sink(file=paste.file, append=T)
        cat(paste(date(), "\n", sep=""))
        print(match.call())
    }

#
# check if all variables are present
    vars <- models.list$vars
    vars.xn <- names(xn)
    nv <- length(vars) 
    for (i in 1:nv) {
        if (any(vars.xn==vars[i]) == F) {stop("explanatory variable '", vars[i], "' not among grid layers of RasterStack xn",  "\n", sep = "")}
    }
    for (i in 1:length(vars.xn) ) {
        if (any(vars==vars.xn[i]) == F) {
            cat(paste("\n", "NOTE: RasterStack layer '", vars.xn[i], "' was not calibrated as explanatory variable", "\n", sep = ""))
            xn <- raster::dropLayer(xn, which(names(xn) == vars.xn[i]))
            xn <- raster::stack(xn)
        }
    }
#
# set minimum and maximum values for xn
    for (i in 1:raster::nlayers(xn)) {
        xn[[i]] <- raster::setMinMax(xn[[i]])
    }

# declare categorical layers for xn
    factors <- models.list$factors
    if(is.null(factors) == F) {
        for (i in 1:length(factors)) {
            j <- which(names(xn) == factors[i])
            xn[[j]] <- raster::as.factor(xn[[j]])
        }
    }
    if(length(factors) == 0) {factors <- NULL}
    dummy.vars <- models.list$dummy.vars
    dummy.vars.noDOMAIN <- models.list$dummy.vars.noDOMAIN 
#
    KML.blur <- trunc(KML.blur)
    if (KML.blur < 1) {KML.blur <- 1}
#
    if (is.null(input.weights) == F) {
        MAXENT <- max(c(input.weights["MAXENT"], -1), na.rm=T)
        MAXLIKE <- max(c(input.weights["MAXLIKE"], -1), na.rm=T)
        GBM <- max(c(input.weights["GBM"], -1), na.rm=T)
        GBMSTEP <- max(c(input.weights["GBMSTEP"], -1), na.rm=T)
        RF <- max(c(input.weights["RF"], -1), na.rm=T)
        GLM <- max(c(input.weights["GLM"], -1), na.rm=T)
        GLMSTEP <- max(c(input.weights["GLMSTEP"], -1), na.rm=T)
        GAM <- max(c(input.weights["GAM"], -1), na.rm=T)
        GAMSTEP <- max(c(input.weights["GAMSTEP"], -1), na.rm=T)
        MGCV <- max(c(input.weights["MGCV"], -1), na.rm=T)
        MGCVFIX <- max(c(input.weights["MGCVFIX"], -1), na.rm=T)
        EARTH <- max(c(input.weights["EARTH"], -1), na.rm=T)
        RPART <- max(c(input.weights["RPART"], -1), na.rm=T)
        NNET <- max(c(input.weights["NNET"], -1), na.rm=T)
        FDA <- max(c(input.weights["FDA"], -1), na.rm=T)
        SVM <- max(c(input.weights["SVM"], -1), na.rm=T)
        SVME <- max(c(input.weights["SVME"], -1), na.rm=T)
        GLMNET <- max(c(input.weights["GLMNET"], -1), na.rm=T)
        BIOCLIM.O <- max(c(input.weights["BIOCLIM.O"], -1), na.rm=T)
        BIOCLIM <- max(c(input.weights["BIOCLIM"], -1), na.rm=T)
        DOMAIN <- max(c(input.weights["DOMAIN"], -1), na.rm=T)
        MAHAL <- max(c(input.weights["MAHAL"], -1), na.rm=T)
        MAHAL01 <- max(c(input.weights["MAHAL01"], -1), na.rm=T)
    }
#
    MAXENT.OLD <- MAXLIKE.OLD <- GBM.OLD <- GBMSTEP.OLD <- RF.OLD <- GLM.OLD <- GLMSTEP.OLD <- GAM.OLD <- GAMSTEP.OLD <- MGCV.OLD <- NULL
    MGCVFIX.OLD <- EARTH.OLD <- RPART.OLD <- NNET.OLD <- FDA.OLD <- SVM.OLD <- SVME.OLD <- GLMNET.OLD <- BIOCLIM.O.OLD <- BIOCLIM.OLD <- DOMAIN.OLD <- MAHAL.OLD <- MAHAL01.OLD <- NULL
# probit models, NULL if no probit model fitted
    MAXENT.PROBIT.OLD <- MAXLIKE.PROBIT.OLD <- GBM.PROBIT.OLD <- GBMSTEP.PROBIT.OLD <- RF.PROBIT.OLD <- GLM.PROBIT.OLD <- GLMSTEP.PROBIT.OLD <- GAM.PROBIT.OLD <- GAMSTEP.PROBIT.OLD <- MGCV.PROBIT.OLD <- NULL
    MGCVFIX.PROBIT.OLD <- EARTH.PROBIT.OLD <- RPART.PROBIT.OLD <- NNET.PROBIT.OLD <- FDA.PROBIT.OLD <- SVM.PROBIT.OLD <- SVME.PROBIT.OLD <- GLMNET.PROBIT.OLD <- BIOCLIM.O.PROBIT.OLD <- BIOCLIM.PROBIT.OLD <- DOMAIN.PROBIT.OLD <- MAHAL.PROBIT.OLD <- MAHAL01.PROBIT.OLD <- NULL
    if (is.null(models.list) == F) {
        if (is.null(models.list$MAXENT) == F) {MAXENT.OLD <- models.list$MAXENT}
        if (is.null(models.list$MAXLIKE) == F) {GBM.OLD <- models.list$MAXLIKE}
        MAXLIKE.formula <- models.list$formulae$MAXLIKE.formula
        if (is.null(models.list$GBM) == F) {GBM.OLD <- models.list$GBM}
        if (is.null(models.list$GBMSTEP) == F) {GBMSTEP.OLD <- models.list$GBMSTEP}
        if (is.null(models.list$RF) == F) {RF.OLD <- models.list$RF}
        if (is.null(models.list$GLM) == F) {GLM.OLD <- models.list$GLM}
        if (is.null(models.list$GLMSTEP) == F) {GLMSTEP.OLD <- models.list$GLMSTEP}
        if (is.null(models.list$GAM) == F) {GAM.OLD <- models.list$GAM}
        if (is.null(models.list$GAMSTEP) == F) {GAMSTEP.OLD <- models.list$GAMSTEP}
        if (is.null(models.list$MGCV) == F) {MGCV.OLD <- models.list$MGCV}
        if (is.null(models.list$MGCVFIX) == F) {MGCVFIX.OLD <- models.list$MGCVFIX}
        if (is.null(models.list$EARTH) == F) {EARTH.OLD <- models.list$EARTH}
        if (is.null(models.list$RPART) == F) {RPART.OLD <- models.list$RPART}
        if (is.null(models.list$NNET) == F) {NNET.OLD <- models.list$NNET}
        if (is.null(models.list$FDA) == F) {FDA.OLD <- models.list$FDA}
        if (is.null(models.list$SVM) == F) {SVM.OLD <- models.list$SVM}
        if (is.null(models.list$SVME) == F) {SVME.OLD <- models.list$SVME}
        if (is.null(models.list$GLMNET) == F) {GLMNET.OLD <- models.list$GLMNET}
        GLMNET.class <- models.list$formulae$GLMNET.class       
        if (is.null(models.list$BIOCLIM.O) == F) {BIOCLIM.O.OLD <- models.list$BIOCLIM.O}
        if (is.null(models.list$BIOCLIM) == F) {BIOCLIM.OLD <- models.list$BIOCLIM}
        if (is.null(models.list$DOMAIN) == F) {DOMAIN.OLD <- models.list$DOMAIN}
        if (is.null(models.list$MAHAL) == F) {MAHAL.OLD <- models.list$MAHAL}
        if (is.null(models.list$MAHAL01) == F) {MAHAL01.OLD <- models.list$MAHAL01}
        MAHAL.shape <- models.list$formulae$MAHAL.shape
# probit models
        if (is.null(models.list$MAXENT.PROBIT) == F) {MAXENT.PROBIT.OLD <- models.list$MAXENT.PROBIT}
        if (is.null(models.list$MAXLIKE.PROBIT) == F) {MAXLIKE.PROBIT.OLD <- models.list$MAXLIKE.PROBIT}
        if (is.null(models.list$GBM.PROBIT) == F) {GBM.PROBIT.OLD <- models.list$GBM.PROBIT}
        if (is.null(models.list$GBMSTEP.PROBIT) == F) {GBMSTEP.PROBIT.OLD <- models.list$GBMSTEP.PROBIT}
        if (is.null(models.list$RF.PROBIT) == F) {RF.PROBIT.OLD <- models.list$RF.PROBIT}
        if (is.null(models.list$GLM.PROBIT) == F) {GLM.PROBIT.OLD <- models.list$GLM.PROBIT}
        if (is.null(models.list$GLMSTEP.PROBIT) == F) {GLMSTEP.PROBIT.OLD <- models.list$GLMSTEP.PROBIT}
        if (is.null(models.list$GAM.PROBIT) == F) {GAM.PROBIT.OLD <- models.list$GAM.PROBIT}
        if (is.null(models.list$GAMSTEP.PROBIT) == F) {GAMSTEP.PROBIT.OLD <- models.list$GAMSTEP.PROBIT}
        if (is.null(models.list$MGCV.PROBIT) == F) {MGCV.PROBIT.OLD <- models.list$MGCV.PROBIT}
        if (is.null(models.list$MGCVFIX.PROBIT) == F) {MGCVFIX.PROBIT.OLD <- models.list$MGCVFIX.PROBIT}
        if (is.null(models.list$EARTH.PROBIT) == F) {EARTH.PROBIT.OLD <- models.list$EARTH.PROBIT}
        if (is.null(models.list$RPART.PROBIT) == F) {RPART.PROBIT.OLD <- models.list$RPART.PROBIT}
        if (is.null(models.list$NNET.PROBIT) == F) {NNET.PROBIT.OLD <- models.list$NNET.PROBIT}
        if (is.null(models.list$FDA.PROBIT) == F) {FDA.PROBIT.OLD <- models.list$FDA.PROBIT}
        if (is.null(models.list$SVM.PROBIT) == F) {SVM.PROBIT.OLD <- models.list$SVM.PROBIT}
        if (is.null(models.list$SVME.PROBIT) == F) {SVME.PROBIT.OLD <- models.list$SVME.PROBIT}
        if (is.null(models.list$GLMNET.PROBIT) == F) {GLMNET.PROBIT.OLD <- models.list$GLMNET.PROBIT}
        if (is.null(models.list$BIOCLIM.O.PROBIT) == F) {BIOCLIM.O.PROBIT.OLD <- models.list$BIOCLIM.O.PROBIT}
        if (is.null(models.list$BIOCLIM.PROBIT) == F) {BIOCLIM.PROBIT.OLD <- models.list$BIOCLIM.PROBIT}
        if (is.null(models.list$DOMAIN.PROBIT) == F) {DOMAIN.PROBIT.OLD <- models.list$DOMAIN.PROBIT}
        if (is.null(models.list$MAHAL.PROBIT) == F) {MAHAL.PROBIT.OLD <- models.list$MAHAL.PROBIT}
        if (is.null(models.list$MAHAL01.PROBIT) == F) {MAHAL01.PROBIT.OLD <- models.list$MAHAL01.PROBIT}
    }
    if (MAXENT > 0) {
	jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
	if (!file.exists(jar)) {stop('maxent program is missing: ', jar, '\nPlease download it here: http://www.cs.princeton.edu/~schapire/maxent/')}
    }
    if (MAXLIKE > 0) {
        if (! requireNamespace("maxlike")) {stop("Please install the maxlike package")}
        MAXLIKE.formula <- ensemble.formulae(xn, factors=factors)$MAXLIKE.formula
        environment(MAXLIKE.formula) <- .BiodiversityR
    }
    if (GBM > 0) {
        if (! requireNamespace("gbm")) {stop("Please install the gbm package")}
	requireNamespace("splines")
    }
    if (RF > 0) {
#  get the probabilities from RF
        predict.RF <- function(object, newdata) {
            p <- predict(object=object, newdata=newdata, type="response")
            return(as.numeric(p))
        }
    }
    if (MGCV > 0 || MGCVFIX > 0) {
# get the probabilities from MGCV
        predict.MGCV <- function(object, newdata, type="response") {
            p <- mgcv::predict.gam(object=object, newdata=newdata, type=type)
            return(as.numeric(p))
        }
    }
    if (EARTH > 0) {
# get the probabilities from earth
        predict.EARTH <- function(object, newdata, type="response") {
            p <- predict(object=object, newdata=newdata, type=type)
            return(as.numeric(p))
        }
    }
    if (NNET > 0) {
# get the probabilities from nnet
        predict.NNET <- function(object, newdata, type="raw") {
            p <- predict(object=object, newdata=newdata, type=type)
            return(as.numeric(p))
        }
    }
    if (SVME > 0) {
#  get the probabilities from svm
        predict.SVME <- function(model, newdata) {
            p <- predict(model, newdata, probability=T)
            return(attr(p, "probabilities")[,1])
         }
    }
    if (GLMNET > 0) {
        if (! requireNamespace("glmnet")) {stop("Please install the glmnet package")}
#  get the mean probabilities from glmnet
        predict.GLMNET <- function(model, newdata, GLMNET.class=FALSE) {
            newdata <- as.matrix(newdata)
            if (GLMNET.class == TRUE) {
                p <- predict(model, newx=newdata, type="class", exact=T)
                n.obs <- nrow(p)
                nv <- ncol(p)
                result <- numeric(n.obs)
                for (i in 1:n.obs) {
                    for (j in 1:nv) {
                        if(p[i, j] == 1) {result[i] <- result[i] + 1}
                    }
                }
                result <- result/nv
                return(result)
             }else{
                p <- predict(model, newx=newdata, type="response", exact=T)
                n.obs <- nrow(p)
                nv <- ncol(p)
                result <- numeric(n.obs)
                for (i in 1:n.obs) {
                    for (j in 1:nv) {
                        result[i] <- result[i] + p[i, j]
                    }
                }
                result <- result/nv
                return(result)
            }
        }
    }
    if (BIOCLIM.O > 0) {
#  get the probabilities for original BIOCLIM
        predict.BIOCLIM.O <- function(object, newdata) {
            lower.limits <- object$lower.limits
            upper.limits <- object$upper.limits
            minima <- object$minima
            maxima <- object$maxima
#
            newdata <- newdata[, which(names(newdata) %in% names(lower.limits)), drop=F]
            result <- as.numeric(rep(NA, nrow(newdata)))
            varnames <- names(newdata)
            nvars <- ncol(newdata)
#
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
    }
    if (MAHAL > 0) {
# get the probabilities from mahal
        predict.MAHAL <- function(model, newdata, PROBIT) {
            p <- dismo::predict(object=model, x=newdata)
            if (PROBIT == F) {
                p[p<0] <- 0
                p[p>1] <- 1
            }
            return(as.numeric(p))
         }
    }
    if (MAHAL01 > 0) {
#  get the probabilities from transformed mahal
        predict.MAHAL01 <- function(model, newdata, MAHAL.shape) {
            p <- dismo::predict(object=model, x=newdata)
            p <- p - 1 - MAHAL.shape
            p <- abs(p)
            p <- MAHAL.shape / p
            return(p)
         }
    }
# 
    output.weights <- input.weights
    prediction.failures <- FALSE
#
# avoid problems with non-existing directories
    dir.create("models", showWarnings = F)
    dir.create("ensembles", showWarnings = F)
    dir.create("ensembles/suitability", showWarnings = F)
    dir.create("ensembles/count", showWarnings = F)
    dir.create("ensembles/presence", showWarnings = F)
    if(KML.out == T) {
        dir.create("kml", showWarnings = F)    
        dir.create("kml/suitability", showWarnings = F)
        dir.create("kml/count", showWarnings = F)
        dir.create("kml/presence", showWarnings = F)
    }
#
    stack.title <- RASTER.stack.name
    if (gsub(".", "_", stack.title, fixed=T) != stack.title) {cat(paste("\n", "WARNING: title of stack (", stack.title, ") contains '.'", "\n\n", sep = ""))}
#
    raster.title <- paste(RASTER.species.name, "_", stack.title , sep="")
    rasterfull <- paste("ensembles//suitability//", raster.title , sep="")
    kmlfull <- paste("kml//suitability//", raster.title , sep="")
    rastercount <- paste("ensembles//count//", raster.title , sep="")
    kmlcount <- paste("kml//count//", raster.title , sep="")
    rasterpresence <- paste("ensembles//presence//", raster.title, sep="")
    kmlpresence <- paste("kml//presence//", raster.title, sep="")
#
    RASTER.species.orig <- RASTER.species.name
    if (RASTER.models.overwrite==T) {
        RASTER.species.name <- "working"
    }else{
        RASTER.species.name <- paste(RASTER.species.name, "_", stack.title, sep="")
    }
#
    cat(paste("\n", "Start of predictions for organism: ", RASTER.species.orig, "\n", sep = ""))
    cat(paste("Predictions for RasterStack: ", stack.title, "\n", sep = ""))
    ensemble.statistics <- NULL
    cat(paste("ensemble raster layers will be saved in folder ", getwd(), "//ensembles", "\n\n", sep = ""))
    statistics.names <- c("n.models", "ensemble.threshold", "ensemble.min", "ensemble.max", "count.min", "count.max") 
    ensemble.statistics <- numeric(6)
    names(ensemble.statistics) <- statistics.names
#
# sometimes still error warnings for minimum and maximum values of the layers
# set minimum and maximum values for xn
    for (i in 1:raster::nlayers(xn)) {
        xn[[i]] <- raster::setMinMax(xn[[i]])
    }

# count models
    mc <- 0
#
# start raster layer creations
    if (output.weights["MAXENT"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Maximum entropy algorithm (package: dismo)\n", sep=""))
        # Put the file 'maxent.jar' in the 'java' folder of dismo
        # the file 'maxent.jar' can be obtained from from http://www.cs.princeton.edu/~schapire/maxent/.
        jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
        results <- MAXENT.OLD
        pmaxent <- NULL
        fullname <- paste("models/", RASTER.species.name, "_MAXENT", sep="")
        tryCatch(pmaxent <- raster::predict(object=results, x=xn, na.rm=TRUE, 
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("MAXENT prediction failed"))},
            silent=F)
        if (is.null(pmaxent) == F) {
            results2 <- MAXENT.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pmaxent, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "MAXENT"
                pmaxent <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pmaxent <- trunc(1000*pmaxent)
            raster::writeRaster(x=pmaxent, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pmaxent, p)/1000
                abs1 <- raster::extract(pmaxent, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["MAXENT"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["MAXENT"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pmaxent, pt)/1000
                abs1 <- raster::extract(pmaxent, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: MAXENT prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["MAXENT"] <- -1
        }
    }
    if (output.weights["MAXLIKE"] > 0) {
        results <- MAXLIKE.OLD
        results$rasters <- models.list$x
        if (all.equal(results$rasters, xn) == F) {
             cat(paste("\n", "WARNING: Maxlike prediction function is problematic with new data", sep = ""))    
             cat(paste("\n", "Maxlike predictions will not be done", sep = "")) 
             cat(paste("\n", "You should recalibrate the ensemble model without MAXLIKE", "\n\n", sep = ""))        
             output.weights["MAXLIKE"]  <- 0
        }
    }
    if (output.weights["MAXLIKE"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Maxlike algorithm (package: maxlike)\n", sep=""))
        results <- MAXLIKE.OLD
#        results$call$formula <- MAXLIKE.formula
#        results$call$rasters <- as.name("x")
#        results$rasters <- models.list$x
#        class(results) <- "maxlikeFit"
#        pmaxlike <- NULL
#         tryCatch(pmaxlike <- predict(results, rasters=xn),
#             error= function(err) {print(paste("MAXLIKE prediction failed"))},
#             silent=F)

# use raster layer from calibration
        pmaxlike <- raster::raster("MAXLIKE_raster")
        fullname <- paste("models//", RASTER.species.name, "_MAXLIKE", sep="")
        results2 <- MAXLIKE.PROBIT.OLD
        if (is.null(results2) == F) {
            cat(paste("Probit transformation", "\n", sep=""))
            fullname2 <- paste("models//", "MAXLIKE_step1", sep="")
            raster::writeRaster(x=pmaxlike, filename=fullname2, progress='text', overwrite=TRUE)
            explan.stack <- raster::stack(fullname2)
            names(explan.stack) <- "MAXLIKE"
            pmaxlike <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                filename=fullname, progress='text', overwrite=TRUE)                
        }
        pmaxlike <- trunc(1000*pmaxlike)
        raster::writeRaster(x=pmaxlike, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        if(evaluate == T) {
            eval1 <- pres1 <- abs1 <- NULL
            cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
            pres1 <- raster::extract(pmaxlike, p)/1000
            abs1 <- raster::extract(pmaxlike, a)/1000
            eval1 <- evaluate(p=pres1, a=abs1)
            print(eval1)
            thresholds.raster["MAXLIKE"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds.raster["MAXLIKE"]))
        }
        if(retest == T) {
            eval1 <- pres1 <- abs1 <- NULL
            cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
            pres1 <- raster::extract(pmaxlike, pt)/1000
            abs1 <- raster::extract(pmaxlike, at)/1000
            eval1 <- evaluate(p=pres1, a=abs1)
            print(eval1)
        }
    }
    if (output.weights["GBM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Generalized boosted regression modeling (package: gbm) \n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "WARRNING: not certain whether the correct factor levels will be used", "\n", sep=""))
        } 
        results <- GBM.OLD
        pgbm <- NULL
        fullname <- paste("models/", RASTER.species.name, "_GBM", sep="")
        tryCatch(pgbm <- raster::predict(object=xn, model=results, fun=gbm::predict.gbm, na.rm=TRUE, factors=factors,
                n.trees=results$n.trees, type="response", 
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("GBM prediction failed"))},
            silent=F)
        if (is.null(pgbm) == F) {
            results2 <- GBM.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pgbm, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "GBM"
                pgbm <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pgbm <- trunc(1000*pgbm)
            raster::writeRaster(x=pgbm, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pgbm, p)/1000
                abs1 <- raster::extract(pgbm, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["GBM"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["GBM"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pgbm, pt)/1000
                abs1 <- raster::extract(pgbm, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: GBM prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["GBM"] <- -1 
        }
    }
    if (output.weights["GBMSTEP"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". gbm step algorithm (package: dismo)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "WARRNING: not certain whether the correct factor levels will be used", "\n", sep=""))
        } 
        results <- GBMSTEP.OLD
        pgbms <- NULL
        fullname <- paste("models/", RASTER.species.name, "_GBMSTEP", sep="")
        tryCatch(pgbms <- raster::predict(object=xn, model=results, fun=gbm::predict.gbm, na.rm=TRUE, factors=factors,
                n.trees=results$n.trees, type="response", 
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("stepwise GBM prediction failed"))},
            silent=F)
        if (is.null(pgbms) == F) {
            results2 <- GBMSTEP.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pgbms, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "GBMSTEP"
                pgbms <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pgbms <- trunc(1000*pgbms)
            raster::writeRaster(x=pgbms, filename=fullname, progress='text', overwrite=TRUE)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pgbms, p)/1000
                abs1 <- raster::extract(pgbms, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["GBMSTEP"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["GBMSTEP"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pgbms, pt)/1000
                abs1 <- raster::extract(pgbms, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: stepwise GBM prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["GBMSTEP"] <- -1 
        }
    }
    if (output.weights["RF"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Random forest algorithm (package: randomForest)\n", sep=""))
        results <- RF.OLD
        prf <- NULL
        fullname <- paste("models/", RASTER.species.name, "_RF", sep="")
        tryCatch(prf <- raster::predict(object=xn, model=results, fun=predict.RF, na.rm=TRUE, factors=factors,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("random forest prediction failed"))},
            silent=F)
        if (is.null(prf) == F) {
            results2 <- RF.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=prf, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "RF"
                prf <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            prf <- trunc(1000*prf)
            raster::writeRaster(x=prf, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(prf, p)/1000
                abs1 <- raster::extract(prf, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["RF"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["RF"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(prf, pt)/1000
                abs1 <- raster::extract(prf, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: random forest prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["RF"] <- -1 
        }
    } 
    if (output.weights["GLM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Generalized Linear Model \n", sep=""))
        results <- GLM.OLD
        pglm <- NULL
        fullname <- paste("models/", RASTER.species.name, "_GLM", sep="")
        tryCatch(pglm <- raster::predict(object=xn, model=results, na.rm=TRUE, type="response", factors=factors,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("GLM prediction failed"))},
            silent=F)
        if (is.null(pglm) == F) {
            results2 <- GLM.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pglm, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "GLM"
                pglm <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pglm <- trunc(1000*pglm)
            raster::writeRaster(x=pglm, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pglm, p)/1000
                abs1 <- raster::extract(pglm, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["GLM"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["GLM"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pglm, pt)/1000
                abs1 <- raster::extract(pglm, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: GLM prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["GLM"] <- -1 
        }
    }
    if (output.weights["GLMSTEP"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Stepwise Generalized Linear Model \n", sep=""))
        results <- GLMSTEP.OLD
        pglms <- NULL
        fullname <- paste("models/", RASTER.species.name, "_GLMSTEP", sep="")
        tryCatch(pglms <- raster::predict(object=xn, model=results, na.rm=TRUE, type="response", factors=factors,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("stepwise GLM prediction failed"))},
            silent=F)
        if (is.null(pglms) == F) {
            results2 <- GLMSTEP.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pglms, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "GLMSTEP"
                pglms <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pglms <- trunc(1000*pglms)
            raster::writeRaster(x=pglms, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pglms, p)/1000
                abs1 <- raster::extract(pglms, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["GLMSTEP"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["GLMSTEP"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pglms, pt)/1000
                abs1 <- raster::extract(pglms, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: stepwise GLM prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["GLMSTEP"] <- -1 
        } 
    }
    if (output.weights["GAM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Generalized Additive Model (package: gam)\n", sep=""))
        results <- GAM.OLD
        pgam <- NULL
        fullname <- paste("models/", RASTER.species.name, "_GAM", sep="")
        tryCatch(pgam <- raster::predict(object=xn, model=results, fun=gam::predict.gam, na.rm=TRUE, type="response", factors=factors,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("GAM (package: gam) prediction failed"))},
            silent=F)
        if (is.null(pgam) == F) {
            results2 <- GAM.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pgam, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "GAM"
                pgam <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pgam <- trunc(1000*pgam)
            raster::writeRaster(x=pgam, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pgam, p)/1000
                abs1 <- raster::extract(pgam, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["GAM"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["GAM"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pgam, pt)/1000
                abs1 <- raster::extract(pgam, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: GAM prediction (gam package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["GAM"] <- -1 
        } 
    }
    if (output.weights["GAMSTEP"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Stepwise Generalized Additive Model (package: gam)\n", sep=""))
        results <- GAMSTEP.OLD
        pgams <- NULL
        fullname <- paste("models/", RASTER.species.name, "_GAMSTEP", sep="")
        tryCatch(pgams <- raster::predict(object=xn, model=results, fun=gam::predict.gam, type="response", na.rm=TRUE, factors=factors,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("stepwise GAM (package: gam) prediction failed"))},
            silent=F)
        if (is.null(pgams) == F) {
            results2 <- GAMSTEP.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pgams, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "GAMSTEP"
                pgams <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pgams <- trunc(1000*pgams)
            raster::writeRaster(x=pgams, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pgams, p)/1000
                abs1 <- raster::extract(pgams, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["GAMSTEP"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["GAMSTEP"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pgams, pt)/1000
                abs1 <- raster::extract(pgams, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: stepwise GAM prediction (gam package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["GAMSTEP"] <- -1 
        } 
    }
    if (output.weights["MGCV"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Generalized Additive Model (package: mgcv)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "WARRNING: not certain whether the correct factor levels will be used", "\n", sep=""))
        } 
        results <- MGCV.OLD
        pmgcv <- NULL
        fullname <- paste("models/", RASTER.species.name, "_MGCV", sep="")
        tryCatch(pmgcv <- raster::predict(object=xn, model=results, fun=predict.MGCV, na.rm=TRUE, type="response", factors=factors,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("GAM (package: mgcv) prediction failed"))},
            silent=F)
        if (is.null(pmgcv) == F) {
            results2 <- MGCV.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pmgcv, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "MGCV"
                pmgcv <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pmgcv <- trunc(1000*pmgcv)
            raster::writeRaster(x=pmgcv, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pmgcv, p)/1000
                abs1 <- raster::extract(pmgcv, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["MGCV"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["MGCV"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pmgcv, pt)/1000
                abs1 <- raster::extract(pmgcv, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: GAM prediction (mgcv package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["MGCV"] <- -1 
        } 
    }
    if (output.weights["MGCVFIX"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". GAM with fixed d.f. regression splines (package: mgcv)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "WARRNING: not certain whether the correct factor levels will be used", "\n", sep=""))
        } 
        results <- MGCVFIX.OLD
        pmgcvf <- NULL
        fullname <- paste("models/", RASTER.species.name, "_MGCVFIX", sep="")
        tryCatch(pmgcvf <- raster::predict(object=xn, model=results, fun=predict.MGCV, na.rm=TRUE, type="response", factors=factors,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("GAM with fixed d.f. regression splines (package: mgcv) prediction failed"))},
            silent=F)
        if (is.null(pmgcvf) == F) {
            results2 <- MGCVFIX.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pmgcvf, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "MGCVFIX"
                pmgcvf <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pmgcvf <- trunc(1000*pmgcvf)
            raster::writeRaster(x=pmgcvf, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pmgcvf, p)/1000
                abs1 <- raster::extract(pmgcvf, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["MGCVFIX"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["MGCVFIX"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pmgcvf, pt)/1000
                abs1 <- raster::extract(pmgcvf, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: GAM prediction (mgcv package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["MGCVFIX"] <- -1 
        } 
    }
    if (output.weights["EARTH"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Multivariate Adaptive Regression Splines (package: earth)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "NOTE: MARS (earth package) with factors may require explicit dummy variables", "\n", sep=""))
        } 
        results <- EARTH.OLD
        pearth <- NULL
        fullname <- paste("models/", RASTER.species.name, "_EARTH", sep="")
        tryCatch(pearth <- raster::predict(object=xn, model=results, fun=predict.EARTH, na.rm=TRUE, type="response", factors=factors,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("MARS (package: earth) prediction failed"))},
            silent=F)
        if (is.null(pearth) == F) {
            results2 <- EARTH.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pearth, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "EARTH"
                pearth <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pearth <- trunc(1000*pearth)
            raster::writeRaster(x=pearth, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pearth, p)/1000
                abs1 <- raster::extract(pearth, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["EARTH"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["EARTH"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pearth, pt)/1000
                abs1 <- raster::extract(pearth, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: MARS prediction (earth package) failed", "\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["EARTH"] <- -1 
        } 
    }
    if (output.weights["RPART"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Recursive Partitioning And Regression Trees (package: rpart)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "WARRNING: not certain whether the correct factor levels will be used", "\n", sep=""))
        } 
        results <- RPART.OLD
        prpart <- NULL
        fullname <- paste("models/", RASTER.species.name, "_RPART", sep="")
        tryCatch(prpart <- raster::predict(object=xn, model=results, na.rm=TRUE, type="prob", index=2, factors=factors,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("RPART prediction failed"))},
            silent=F)
        if (is.null(prpart) == F) {
            results2 <- RPART.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=prpart, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "RPART"
                prpart <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            prpart <- trunc(1000*prpart)
            raster::writeRaster(x=prpart, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(prpart, p)/1000
                abs1 <- raster::extract(prpart, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["RPART"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["RPART"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(prpart, pt)/1000
                abs1 <- raster::extract(prpart, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: RPART prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["RPART"] <- -1 
        } 
    }
    if (output.weights["NNET"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Artificial Neural Network (package: nnet)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "WARRNING: not certain whether the correct factor levels will be used", "\n", sep=""))
        } 
        results <- NNET.OLD
        pnnet <- NULL
        fullname <- paste("models/", RASTER.species.name, "_NNET", sep="")
        tryCatch(pnnet <- raster::predict(object=xn, model=results, fun=predict.NNET, na.rm=TRUE, factors=factors,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("Artificial Neural Network (package: nnet) prediction failed"))},
            silent=F)
        if (is.null(pnnet) == F) {
            results2 <- NNET.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pnnet, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "NNET"
                pnnet <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pnnet <- trunc(1000*pnnet)
            raster::writeRaster(x=pnnet, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pnnet, p)/1000
                abs1 <- raster::extract(pnnet, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["NNET"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["NNET"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pnnet, pt)/1000
                abs1 <- raster::extract(pnnet, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: ANN prediction (nnet package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["NNET"] <- -1 
        } 
    }
    if (output.weights["FDA"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Flexible Discriminant Analysis (package: mda)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "WARRNING: not certain whether the correct factor levels will be used", "\n", sep=""))
        } 
        results <- FDA.OLD
        pfda <- NULL
        fullname <- paste("models/", RASTER.species.name, "_FDA", sep="")
        tryCatch(pfda <- raster::predict(object=xn, model=results, na.rm=TRUE, type="posterior", index=2, factors=factors,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("FDA prediction failed"))},
            silent=F)
        if (is.null(pfda) == F) {
            results2 <- FDA.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pfda, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "FDA"
                pfda <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pfda <- trunc(1000*pfda)
            raster::writeRaster(x=pfda, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pfda, p)/1000
                abs1 <- raster::extract(pfda, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["FDA"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["FDA"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pfda, pt)/1000
                abs1 <- raster::extract(pfda, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: FDA prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["FDA"] <- -1 
        } 
    }
    if (output.weights["SVM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Support Vector Machines (package: kernlab)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "NOTE: SVM model with factors may require explicit dummy variables", "\n", sep=""))
        } 
        results <- SVM.OLD
        psvm <- NULL
        fullname <- paste("models/", RASTER.species.name, "_SVM", sep="")
        predict.svm2 <- as.function(kernlab::predict)
        tryCatch(psvm <- raster::predict(object=xn, model=results, fun=predict.svm2, na.rm=TRUE, type="probabilities", index=2, factors=factors,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("Support Vector Machines (package: kernlab) prediction failed"))},
            silent=F)
        if (is.null(psvm) == F) {
            results2 <- SVM.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=psvm, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "SVM"
                psvm <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            psvm <- trunc(1000*psvm)
            raster::writeRaster(x=psvm, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(psvm, p)/1000
                abs1 <- raster::extract(psvm, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["SVM"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["SVM"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(psvm, pt)/1000
                abs1 <- raster::extract(psvm, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: SVM prediction (kernlab package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["SVM"] <- -1 
        } 
    }
    if (output.weights["SVME"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Support Vector Machines (package: e1071)\n", sep=""))
        results <- SVME.OLD
        psvme <- NULL
        fullname <- paste("models/", RASTER.species.name, "_SVME", sep="")
        tryCatch(psvme <- raster::predict(object=xn, model=results, fun=predict.SVME, na.rm=TRUE, factors=factors,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("SVM prediction (e1071 package) failed"))},
            warning= function(war) {print(paste("SVM prediction (e1071 package) failed"))},
            silent=F)
        if (is.null(psvme) == F) {
            results2 <- SVME.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=psvme, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "SVME"
                psvme <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            psvme <- trunc(1000*psvme)
            raster::writeRaster(x=psvme, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(psvme, p)/1000
                abs1 <- raster::extract(psvme, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["SVME"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["SVME"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(psvme, pt)/1000
                abs1 <- raster::extract(psvme, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: SVM prediction (e1071 package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["SVME"] <- -1 
        }
    }
    if (output.weights["GLMNET"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". GLM with lasso or elasticnet regularization (package: glmnet)\n", sep=""))
        if (is.null(factors) == F) {
            cat(paste("\n", "NOTE: factors not considered (maybe consider dummy variables)", "\n", sep=""))
        }
        results <- GLMNET.OLD
        pglmnet <- NULL
        fullname <- paste("models/", RASTER.species.name, "_GLMNET", sep="")
        tryCatch(pglmnet <- raster::predict(object=xn, model=results, fun=predict.GLMNET, na.rm=TRUE, GLMNET.class=GLMNET.class,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("GLMNET prediction (glmnet package) failed"))},
            warning= function(war) {print(paste("GLMNET prediction (glmnet package) failed"))},
            silent=F)
        if (is.null(pglmnet) == F) {
            results2 <- GLMNET.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pglmnet, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "GLMNET"
                pglmnet <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pglmnet <- trunc(1000*pglmnet)
            raster::writeRaster(x=pglmnet, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pglmnet, p)/1000
                abs1 <- raster::extract(pglmnet, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["GLMNET"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["GLMNET"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pglmnet, pt)/1000
                abs1 <- raster::extract(pglmnet, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: GLMNET prediction (glmnet package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["GLMNET"] <- -1 
        }
    }
    if (output.weights["BIOCLIM.O"] > 0) {  
        mc <- mc+1
        cat(paste("\n", mc, ". original BIOCLIM algorithm (package: BiodiversityR)\n", sep=""))
        results <- BIOCLIM.O.OLD
        pbioO <- NULL
        fullname <- paste("models/", RASTER.species.name, "_BIOCLIMO", sep="")
        tryCatch(pbioO <- raster::predict(object=xn, model=results, fun=predict.BIOCLIM.O, na.rm=TRUE, 
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("original BIOCLIM prediction failed"))},
            silent=F)
        if (is.null(pbioO) == F) {
            results2 <- BIOCLIM.O.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pbioO, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "BIOCLIM.O"
                pbio <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pbioO <- trunc(1000*pbioO)
            raster::writeRaster(x=pbioO, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pbioO, p)/1000
                abs1 <- raster::extract(pbioO, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["BIOCLIM.O"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["BIOCLIM.O"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pbioO, pt)/1000
                abs1 <- raster::extract(pbioO, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: original BIOCLIM prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["BIOCLIM.O"] <- -1 
        }
    }
    if (output.weights["BIOCLIM"] > 0 || output.weights["DOMAIN"] > 0 || output.weights["MAHAL"] > 0 || output.weights["MAHAL01"] > 0) {
        if(is.null(factors) == F) {
            xn <- raster::dropLayer(xn, which(names(xn) %in% factors)) 
            xn <- raster::stack(xn)              
        }
    }
    if (output.weights["BIOCLIM"] > 0) {  
        mc <- mc+1
        cat(paste("\n", mc, ". BIOCLIM algorithm (package: dismo)\n", sep=""))
        results <- BIOCLIM.OLD
        pbio <- NULL
        fullname <- paste("models/", RASTER.species.name, "_BIOCLIM", sep="")
        tryCatch(pbio <- dismo::predict(object=results, x=xn, na.rm=TRUE, 
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("BIOCLIM prediction failed"))},
            silent=F)
        if (is.null(pbio) == F) {
            results2 <- BIOCLIM.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pbio, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "BIOCLIM"
                pbio <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pbio <- trunc(1000*pbio)
            raster::writeRaster(x=pbio, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pbio, p)/1000
                abs1 <- raster::extract(pbio, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["BIOCLIM"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["BIOCLIM"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pbio, pt)/1000
                abs1 <- raster::extract(pbio, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: BIOCLIM prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["BIOCLIM"] <- -1 
        }
    }
    if (output.weights["DOMAIN"] > 0) {
        if(is.null(factors) == F) {
            xn <- raster::dropLayer(xn, which(names(xn) %in% dummy.vars.noDOMAIN)) 
            xn <- raster::stack(xn)              
        }
    }
    if (output.weights["DOMAIN"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". DOMAIN algorithm (package: dismo)\n", sep=""))
        if(is.null(models.list$dummy.vars.noDOMAIN) == F) {
            xn <- raster::dropLayer(xn, which(names(xn) %in% models.list$dummy.vars.noDOMAIN))
            xn <- raster::stack(xn)          
        }
        results <- DOMAIN.OLD
        pdom <- NULL
        fullname <- paste("models/", RASTER.species.name, "_DOMAIN", sep="")
        tryCatch(pdom <- dismo::predict(object=results, x=xn, na.rm=TRUE, 
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("DOMAIN prediction failed"))},
            silent=F)
        if (is.null(pdom) == F) {
            results2 <- DOMAIN.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pdom, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "DOMAIN"
                pdom <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pdom <- trunc(1000*pdom)
            raster::writeRaster(x=pdom, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pdom, p)/1000
                abs1 <- raster::extract(pdom, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["DOMAIN"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["DOMAIN"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pdom, pt)/1000
                abs1 <- raster::extract(pdom, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: DOMAIN prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["DOMAIN"] <- -1 
        }
    }
    if (output.weights["MAHAL"] > 0 || output.weights["MAHAL01"] > 0) {  
        if(is.null(dummy.vars) == F) {
            xn <- raster::dropLayer(xn, which(names(xn) %in% dummy.vars))
            xn <- raster::stack(xn)            
        }
    }
    if (output.weights["MAHAL"] > 0) {  
        mc <- mc+1
        cat(paste("\n", mc, ". Mahalanobis algorithm (package: dismo)\n", sep=""))
        results <- MAHAL.OLD
        pmahal <- NULL
        fullname <- paste("models/", RASTER.species.name, "_MAHAL", sep="")
# not possible to use the predict.mahal function as raster::predict automatically reverts to dismo::predict for 'DistModel' objects
        results2 <- MAHAL.PROBIT.OLD
        if (is.null(results2) == F) {
            tryCatch(pmahal <- raster::predict(object=xn, model=results, fun=predict.MAHAL, na.rm=TRUE, PROBIT=FALSE,
                    filename=fullname, progress='text', overwrite=TRUE),
                error= function(err) {print(paste("Mahalanobis prediction failed"))},
                silent=F)
        }else{
            tryCatch(pmahal <- raster::predict(object=xn, model=results, fun=predict.MAHAL, na.rm=TRUE, PROBIT=TRUE,
                    filename=fullname, progress='text', overwrite=TRUE),
                error= function(err) {print(paste("Mahalanobis prediction failed"))},
                silent=F)
        }
        if (is.null(pmahal) == F) {
            results2 <- MAHAL.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pmahal, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "MAHAL"
                pmahal <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pmahal <- trunc(1000*pmahal)
            raster::writeRaster(x=pmahal, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pmahal, p)/1000
                abs1 <- raster::extract(pmahal, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["MAHAL"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["MAHAL"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pmahal, pt)/1000
                abs1 <- raster::extract(pmahal, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: Mahalanobis prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["MAHAL"] <- -1 
        }
    }
    if (output.weights["MAHAL01"] > 0) {  
        mc <- mc+1
        cat(paste("\n", mc, ". Mahalanobis algorithm (transformed within 0 to 1 interval)", "\n", sep=""))
        if(is.null(dummy.vars) == F) {
            xn <- raster::dropLayer(xn, which(names(xn) %in% dummy.vars))
            xn <- raster::stack(xn)            
        }
        results <- MAHAL01.OLD
        pmahal <- NULL
        fullname <- paste("models/", RASTER.species.name, "_MAHAL01", sep="")
# not possible to use the predict.mahal function as raster::predict automatically reverts to dismo::predict for 'DistModel' objects
        tryCatch(pmahal01 <- raster::predict(object=xn, model=results, fun=predict.MAHAL01, na.rm=TRUE, MAHAL.shape=MAHAL.shape,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("transformed Mahalanobis prediction failed"))},
            silent=F)
        if (is.null(pmahal01) == F) {
#            pmahal <- pmahal - 1 - MAHAL.shape
#            pmahal <- abs(pmahal)
#            pmahal <- MAHAL.shape / pmahal
            results2 <- MAHAL01.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pmahal01, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "MAHAL01"
                pmahal <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pmahal01 <- trunc(1000*pmahal01)
            raster::writeRaster(x=pmahal01, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pmahal01, p)/1000
                abs1 <- raster::extract(pmahal01, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["MAHAL01"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["MAHAL01"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pmahal01, pt)/1000
                abs1 <- raster::extract(pmahal01, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: transformed Mahalanobis prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["MAHAL01"] <- -1 
        }
    }
#
    if (prediction.failures == T) {
        cat(paste("\n", "WARNING: some predictions failed", sep = ""))
        cat(paste("\n", "actual weights that were used were (-1 indicates failed predictions):", "\n", sep = ""))
        print(output.weights)
        cat(paste("\n", "Because ensemble suitability would therefore underestimated", sep = ""))
        cat(paste("\n", "the ensemble will not be calibrated", "\n", sep = ""))
    }
#
# create ensembles if no prediction failures

    if (prediction.failures == F) {

    if (evaluate == T) {thresholds <- thresholds.raster}
    cat(paste("\n", "submodel thresholds for absence-presence used: ", "\n", sep = ""))
    if (evaluate == T) {
        cat(paste("(thresholds recalculated from raster layers)", "\n", sep = ""))
        thresholds <- thresholds.raster
    }
    print(thresholds)
#
    mc <- mc+1
    cat(paste("\n\n", mc, ". Ensemble algorithm\n", sep=""))
    ensemble.statistics["n.models"] <- sum(as.numeric(output.weights > 0))
    ensemble <- xn[[1]] == raster::NAvalue(xn[[1]])
    raster::setMinMax(ensemble)
    names(ensemble) <- raster.title
    raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
    enscount <- ensemble
    raster::setMinMax(enscount)
    names(enscount) <- paste(raster.title, "_count", sep="")
    raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    enspresence <- ensemble
    raster::setMinMax(enspresence)
    names(enspresence) <- paste(raster.title, "_presence", sep="")
    raster::writeRaster(x=enspresence, filename=rasterpresence, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    if (output.weights["MAXENT"] > 0) {
        ensemble <- ensemble + output.weights["MAXENT"] * pmaxent
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)       
        pmaxent <- pmaxent >= 1000 * thresholds["MAXENT"]
        enscount <- enscount + pmaxent
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (output.weights["MAXLIKE"] > 0) {
        ensemble <- ensemble + output.weights["MAXLIKE"] * pmaxlike
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)       
        pmaxent <- pmaxlike >= 1000 * thresholds["MAXLIKE"]
        enscount <- enscount + pmaxlike
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (output.weights["GBM"] > 0) {
        ensemble <- ensemble + output.weights["GBM"] * pgbm
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)      
        pgbm <- pgbm >= 1000 * thresholds["GBM"]
        enscount <- enscount + pgbm
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (output.weights["GBMSTEP"] > 0) {
        ensemble <- ensemble + output.weights["GBMSTEP"] * pgbms
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)      
        pgbms <- pgbms >= 1000 * thresholds["GBMSTEP"]
        enscount <- enscount + pgbms
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (output.weights["RF"] > 0) {
        ensemble <- ensemble + output.weights["RF"] * prf
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)      
        prf <- prf >= 1000 * thresholds["RF"]
        enscount <- enscount + prf
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (output.weights["GLM"] > 0) {
        ensemble <- ensemble + output.weights["GLM"] * pglm
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)    
        pglm <- pglm >= 1000 * thresholds["GLM"]
        enscount <- enscount + pglm
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (output.weights["GLMSTEP"] > 0) {
        ensemble <- ensemble + output.weights["GLMSTEP"] * pglms
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)      
        pglms <- pglms >= 1000 * thresholds["GLMSTEP"]
        enscount <- enscount + pglms
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (output.weights["GAM"] > 0) {
        ensemble <- ensemble + output.weights["GAM"] * pgam
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pgam <- pgam >= 1000 * thresholds["GAM"]
        enscount <- enscount + pgam
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (output.weights["GAMSTEP"] > 0) {
        ensemble <- ensemble + output.weights["GAMSTEP"] * pgams
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pgams <- pgams >= 1000 * thresholds["GAMSTEP"]
        enscount <- enscount + pgams
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (output.weights["MGCV"] > 0) {
        ensemble <- ensemble + output.weights["MGCV"] * pmgcv
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pmgcv <- pmgcv >= 1000 * thresholds["MGCV"]
        enscount <- enscount + pmgcv
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (output.weights["MGCVFIX"] > 0) {
        ensemble <- ensemble + output.weights["MGCVFIX"] * pmgcvf
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pmgcvf <- pmgcvf >= 1000 * thresholds["MGCVFIX"]
        enscount <- enscount + pmgcvf
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (output.weights["EARTH"] > 0) {
        ensemble <- ensemble + output.weights["EARTH"] * pearth
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pearth <- pearth >= 1000 * thresholds["EARTH"]
        enscount <- enscount + pearth
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (output.weights["RPART"] > 0) {
        ensemble <- ensemble + output.weights["RPART"] * prpart
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        prpart <- prpart >= 1000 * thresholds["RPART"]
        enscount <- enscount + prpart
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (output.weights["NNET"] > 0) {
        ensemble <- ensemble + output.weights["NNET"] * pnnet
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pnnet <- pnnet >= 1000 * thresholds["NNET"]
        enscount <- enscount + pnnet
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (output.weights["FDA"] > 0) {
        ensemble <- ensemble + output.weights["FDA"] * pfda
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pfda <- pfda >= 1000 * thresholds["FDA"]
        enscount <- enscount + pfda
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (output.weights["SVM"] > 0) {
        ensemble <- ensemble + output.weights["SVM"] * psvm
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        psvm <- psvm >= 1000 * thresholds["SVM"]
        enscount <- enscount + psvm
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (output.weights["SVME"] > 0) {
        ensemble <- ensemble + output.weights["SVME"] * psvme
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        psvme <- psvme >= 1000 * thresholds["SVME"]
        enscount <- enscount + psvme
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (output.weights["GLMNET"] > 0) {
        ensemble <- ensemble + output.weights["GLMNET"] * pglmnet
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pglmnet <- pglmnet >= 1000 * thresholds["GLMNET"]
        enscount <- enscount + pglmnet
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (output.weights["BIOCLIM.O"] > 0) {
        ensemble <- ensemble + output.weights["BIOCLIM.O"] * pbioO
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pbioO <- pbioO >= 1000 * thresholds["BIOCLIM.O"]
        enscount <- enscount + pbioO
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (output.weights["BIOCLIM"] > 0) {
        ensemble <- ensemble + output.weights["BIOCLIM"] * pbio
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pbio <- pbio >= 1000 * thresholds["BIOCLIM"]
        enscount <- enscount + pbio
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (output.weights["DOMAIN"] > 0) {
        ensemble <- ensemble + output.weights["DOMAIN"] * pdom
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pdom <- pdom >= 1000 * thresholds["DOMAIN"]
        enscount <- enscount + pdom
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (output.weights["MAHAL"] > 0) {
        ensemble <- ensemble + output.weights["MAHAL"] * pmahal
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pmahal <- pmahal >= 1000 * thresholds["MAHAL"]
        enscount <- enscount + pmahal
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (output.weights["MAHAL01"] > 0) {
        ensemble <- ensemble + output.weights["MAHAL01"] * pmahal01
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pmahal01 <- pmahal01 >= 1000 * thresholds["MAHAL01"]
        enscount <- enscount + pmahal01
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
#
# note that submodels had already been multiplied by 1000
    ensemble <- trunc(ensemble)
    raster::setMinMax(ensemble)
    ensemble.statistics["ensemble.min"] <- raster::minValue(ensemble)
    ensemble.statistics["ensemble.max"] <- raster::maxValue(ensemble)
#    names(ensemble) <- raster.title
#    raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
#  avoid possible problems with saving of names of the raster layers
    raster::writeRaster(ensemble, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- raster.title
    raster::writeRaster(working.raster, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
#
    raster::setMinMax(enscount)
    ensemble.statistics["count.min"] <- raster::minValue(enscount)
    ensemble.statistics["count.max"] <- raster::maxValue(enscount)
#    names(enscount) <- paste(raster.title, "_count", sep="")
#    raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
#  avoid possible problems with saving of names of the raster layers
    raster::writeRaster(enscount, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- paste(raster.title, "_count", sep="")
    raster::writeRaster(working.raster, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
#
    if (KML.out == T) {
        raster::writeRaster(enscount, filename="KMLworking.grd", overwrite=T)
        KMLworking.raster <- raster::raster("KMLworking.grd")
        names(KMLworking.raster) <- paste(raster.title, "_count", sep="")
        nmax <- sum(as.numeric(output.weights > 0))
        if (nmax > 3) {
            raster::KML(KMLworking.raster, filename=kmlcount, col=c("grey", "black", grDevices::rainbow(n=(nmax-2), start=0, end=1/3), "blue"),
                colNA=0, blur=10, overwrite=T, breaks=seq(from=-1, to=nmax, by=1))
        }else{
            raster::KML(KMLworking.raster, filename=kmlcount, col=c("grey", grDevices::rainbow(n=nmax, start=0, end=1/3)),
                colNA=0, blur=10, overwrite=TRUE, breaks=seq(from=-1, to=nmax, by=1))
        }
    }
#
    if(evaluate == T) {
        eval1 <- NULL
        cat(paste("\n", "Evaluation of created ensemble raster layer (", rasterfull, ") at locations p and a", "\n\n", sep = ""))
        pres_consensus <- raster::extract(ensemble, p)/1000
        abs_consensus <- raster::extract(ensemble, a)/1000
        eval1 <- evaluate(p=pres_consensus, a=abs_consensus)
        print(eval1)
        thresholds["ENSEMBLE"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
            threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres_consensus, Abs=abs_consensus)
        cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
        print(as.numeric(thresholds["ENSEMBLE"]))
    }
    if(retest == T) {
        eval1 <- NULL
        cat(paste("\n", "Evaluation of created ensemble raster layer (", rasterfull, ") at locations pt and at", "\n\n", sep = ""))
        pres_consensus <- raster::extract(ensemble, pt)/1000
        abs_consensus <- raster::extract(ensemble, at)/1000
        eval1 <- evaluate(p=pres_consensus, a=abs_consensus)
        print(eval1)
    }
    ensemble.statistics["ensemble.threshold"] <- thresholds["ENSEMBLE"]
#
    if (KML.out == T) {
        raster::writeRaster(ensemble, filename="KMLworking.grd", overwrite=T)
        KMLworking.raster <- raster::raster("KMLworking.grd")
        raster::setMinMax(KMLworking.raster)
        names(KMLworking.raster) <- raster.title
        raster::writeRaster(KMLworking.raster, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        thresholdx <- 1000 * as.numeric(thresholds["ENSEMBLE"])
        raster.min <- raster::minValue(KMLworking.raster)
        raster.max <- raster::maxValue(KMLworking.raster)
        abs.breaks <- 8
        pres.breaks <- 8
        seq1 <- round(seq(from=raster.min, to=thresholdx, length.out=abs.breaks), 4)
        seq1 <- seq1[1:(abs.breaks-1)]
        seq1[-abs.breaks]
        seq1 <- unique(seq1)
        seq2 <- round(seq(from = thresholdx, to = raster.max, length.out=pres.breaks), 4)
        seq2 <- unique(seq2)
        raster::KML(KMLworking.raster, filename=kmlfull, breaks = c(seq1, seq2), col = c(grDevices::rainbow(n=length(seq1), start=0, end =1/6), grDevices::rainbow(n=length(seq2)-1, start=3/6, end=4/6)), colNA = 0, 
            blur=KML.blur, maxpixels=KML.maxpixels, overwrite=TRUE)
    }
    enspresence <- ensemble >= 1000 * thresholds["ENSEMBLE"]
    raster::setMinMax(enspresence)
#    names(enspresence) <- paste(raster.title, "_presence", sep="")
#    raster::writeRaster(x=enspresence, filename=rasterpresence, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
#  avoid possible problems with saving of names of the raster layers
    raster::writeRaster(enspresence, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- paste(raster.title, "_presence", sep="")
    raster::writeRaster(working.raster, filename=rasterpresence, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
#
    if (KML.out == T) {
        raster::writeRaster(enspresence, filename="KMLworking.grd", overwrite=T)
        KMLworking.raster <- raster::raster("KMLworking.grd")
        names(KMLworking.raster) <- paste(raster.title, "_presence", sep="")
        raster::KML(KMLworking.raster, filename=kmlpresence, col=c("grey", "green"),
            colNA=0, blur=KML.blur, maxpixels=KML.maxpixels, overwrite=T)
    }
#
    cat(paste("\n", "End of modelling for organism: ", RASTER.species.orig, "\n\n", sep = ""))
    cat(paste("Predictions were made for RasterStack: ", stack.title, "\n\n", sep = ""))

    out <- list(ensemble.statistics=ensemble.statistics, output.weights=output.weights, thresholds=thresholds, call=match.call() )
    if (SINK==T && OLD.SINK==F) {sink(file=NULL, append=T)}
    return(out)

# end of prediction failures loop
}else{
    out  <- list(warning="prediction failure for some algorithms", output.weights=output.weights, call=match.call() )
    if (SINK==T && OLD.SINK==F) {sink(file=NULL, append=T)}
    return(out)
}

}

