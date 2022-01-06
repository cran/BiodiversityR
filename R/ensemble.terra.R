`ensemble.terra` <- function(
    xn=NULL, 
    models.list=NULL, 
    input.weights=models.list$output.weights,
    thresholds=models.list$thresholds,
    RASTER.species.name=models.list$species.name, 
    RASTER.stack.name="xnTitle", 
    RASTER.filetype="GTiff", RASTER.datatype="INT2S", RASTER.NAflag=-32767,
    RASTER.models.overwrite=TRUE,
    evaluate=FALSE, SINK=FALSE,
    p=models.list$p, a=models.list$a,
    pt=models.list$pt, at=models.list$at,
    CATCH.OFF=FALSE
)
{
    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if (is.null(xn) == T) {stop("value for parameter xn is missing (SpatRaster object)")}
    if(inherits(xn, "SpatRaster") == F) {stop("xn is not a SpatRaster object")}
    if (is.null(models.list) == T) {stop("provide 'models.list' as models will not be recalibrated and retested")}
    if (is.null(input.weights) == T) {input.weights <- models.list$output.weights}
    if (is.null(thresholds) == T) {stop("provide 'thresholds' as models will not be recalibrated and retested")}
    thresholds.raster <- thresholds
    if(is.null(p) == F) {
        p <- data.frame(p)
        names(p) <- c("x", "y")
    }
    if(is.null(a) == F) {
        a <- data.frame(a)
        names(a) <- c("x", "y")
    }
    if(is.null(pt) == F) {
        pt <- data.frame(pt)
        names(pt) <- c("x", "y")
    }
    if(is.null(at) == F) {
        at <- data.frame(at)
        names(at) <- c("x", "y")
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
            xn <- terra::subset(xn, which(names(xn) == vars.xn[i]))
            xn <- terra::rast(xn)
        }
    }
#
# declare categorical layers for xn
    factors <- models.list$factors
    if(is.null(factors) == F) {
        for (i in 1:length(factors)) {
            j <- which(names(xn) == factors[i])
            xn[[j]] <- terra::as.factor(xn[[j]])
        }
    }
    if(length(factors) == 0) {factors <- NULL}
    factlevels <- models.list$factlevels
    dummy.vars <- models.list$dummy.vars
    dummy.vars.noDOMAIN <- models.list$dummy.vars.noDOMAIN 
#
    if (is.null(input.weights) == F) {
        MAXENT <- max(c(input.weights["MAXENT"], -1), na.rm=T)
        MAXNET <- max(c(input.weights["MAXNET"], -1), na.rm=T)
        MAXLIKE <- max(c(input.weights["MAXLIKE"], -1), na.rm=T)
        GBM <- max(c(input.weights["GBM"], -1), na.rm=T)
        GBMSTEP <- max(c(input.weights["GBMSTEP"], -1), na.rm=T)
        RF <- max(c(input.weights["RF"], -1), na.rm=T)
        CF <- max(c(input.weights["CF"], -1), na.rm=T)
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
    MAXENT.OLD <- MAXNET.OLD <- MAXLIKE.OLD <- GBM.OLD <- GBMSTEP.OLD <- RF.OLD <- CF.OLD <- GLM.OLD <- GLMSTEP.OLD <- GAM.OLD <- GAMSTEP.OLD <- MGCV.OLD <- NULL
    MGCVFIX.OLD <- EARTH.OLD <- RPART.OLD <- NNET.OLD <- FDA.OLD <- SVM.OLD <- SVME.OLD <- GLMNET.OLD <- BIOCLIM.O.OLD <- BIOCLIM.OLD <- DOMAIN.OLD <- MAHAL.OLD <- MAHAL01.OLD <- NULL
# probit models, NULL if no probit model fitted
    MAXENT.PROBIT.OLD <- MAXNET.PROBIT.OLD <- MAXLIKE.PROBIT.OLD <- GBM.PROBIT.OLD <- GBMSTEP.PROBIT.OLD <- RF.PROBIT.OLD <- CF.PROBIT.OLD <- GLM.PROBIT.OLD <- GLMSTEP.PROBIT.OLD <- GAM.PROBIT.OLD <- GAMSTEP.PROBIT.OLD <- MGCV.PROBIT.OLD <- NULL
    MGCVFIX.PROBIT.OLD <- EARTH.PROBIT.OLD <- RPART.PROBIT.OLD <- NNET.PROBIT.OLD <- FDA.PROBIT.OLD <- SVM.PROBIT.OLD <- SVME.PROBIT.OLD <- GLMNET.PROBIT.OLD <- BIOCLIM.O.PROBIT.OLD <- BIOCLIM.PROBIT.OLD <- DOMAIN.PROBIT.OLD <- MAHAL.PROBIT.OLD <- MAHAL01.PROBIT.OLD <- NULL
    if (is.null(models.list) == F) {
        if (is.null(models.list$MAXENT) == F) {MAXENT.OLD <- models.list$MAXENT}
        if (is.null(models.list$MAXNET) == F) {
            MAXNET.OLD <- models.list$MAXNET
            MAXNET.clamp <- models.list$formulae$MAXNET.clamp
            MAXNET.type <- models.list$formulae$MAXNET.type
        }
        if (is.null(models.list$MAXLIKE) == F) {
            MAXLIKE.OLD <- models.list$MAXLIKE
            MAXLIKE.formula <- models.list$formulae$MAXLIKE.formula
        }
        if (is.null(models.list$GBM) == F) {GBM.OLD <- models.list$GBM}
        if (is.null(models.list$GBMSTEP) == F) {GBMSTEP.OLD <- models.list$GBMSTEP}
        if (is.null(models.list$RF) == F) {RF.OLD <- models.list$RF}
        if (is.null(models.list$CF) == F) {CF.OLD <- models.list$CF}
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
        if (is.null(models.list$GLMNET) == F) {
            GLMNET.OLD <- models.list$GLMNET
            GLMNET.class <- models.list$formulae$GLMNET.class
        }
        if (is.null(models.list$BIOCLIM.O) == F) {BIOCLIM.O.OLD <- models.list$BIOCLIM.O}
        if (is.null(models.list$BIOCLIM) == F) {BIOCLIM.OLD <- models.list$BIOCLIM}
        if (is.null(models.list$DOMAIN) == F) {DOMAIN.OLD <- models.list$DOMAIN}
        if (is.null(models.list$MAHAL) == F) {MAHAL.OLD <- models.list$MAHAL}
        if (is.null(models.list$MAHAL01) == F) {
            MAHAL01.OLD <- models.list$MAHAL01
            MAHAL.shape <- models.list$formulae$MAHAL.shape
        }
# probit models
        if (is.null(models.list$MAXENT.PROBIT) == F) {MAXENT.PROBIT.OLD <- models.list$MAXENT.PROBIT}
        if (is.null(models.list$MAXNET.PROBIT) == F) {MAXNET.PROBIT.OLD <- models.list$MAXNET.PROBIT}
        if (is.null(models.list$MAXLIKE.PROBIT) == F) {MAXLIKE.PROBIT.OLD <- models.list$MAXLIKE.PROBIT}
        if (is.null(models.list$GBM.PROBIT) == F) {GBM.PROBIT.OLD <- models.list$GBM.PROBIT}
        if (is.null(models.list$GBMSTEP.PROBIT) == F) {GBMSTEP.PROBIT.OLD <- models.list$GBMSTEP.PROBIT}
        if (is.null(models.list$RF.PROBIT) == F) {RF.PROBIT.OLD <- models.list$RF.PROBIT}
        if (is.null(models.list$CF.PROBIT) == F) {CF.PROBIT.OLD <- models.list$CF.PROBIT}
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
    if (MAXNET > 0) {
        if (! requireNamespace("maxnet")) {stop("Please install the maxnet package")}
            predict.maxnet2 <- function(object, newdata, clamp=F, type=c("cloglog")) {
                p <- predict(object=object, newdata=newdata, clamp=clamp, type=type)
                return(as.numeric(p))
            }
    }
    if (MAXLIKE > 0) {
        if (! requireNamespace("maxlike")) {stop("Please install the maxlike package")}
#        MAXLIKE.formula <- ensemble.formulae(xn, factors=factors)$MAXLIKE.formula
#        environment(MAXLIKE.formula) <- .BiodiversityR
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
    if (CF > 0) {
#  get the probabilities from RF
#  ensure that cases with missing values are removed
        if (! requireNamespace("party")) {stop("Please install the party package")}
        predict.CF <- function(object, newdata) {
# avoid problems with single variables, especially with terra::predict
            for (i in 1:ncol(newdata)) {
                if (is.integer(newdata[, i])) {newdata[, i] <- as.numeric(newdata[, i])}
            }
            p1 <- predict(object=object, newdata=newdata, type="prob")
            p <- numeric(length(p1))
            for (i in 1:length(p1)) {p[i] <- p1[[i]][2]}
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
#
    stack.title <- RASTER.stack.name
    if (gsub(".", "_", stack.title, fixed=T) != stack.title) {cat(paste("\n", "WARNING: title of stack (", stack.title, ") contains '.'", "\n\n", sep = ""))}
#
    raster.title <- paste(RASTER.species.name, "_", stack.title , sep="")
    rasterfull <- paste("ensembles//suitability//", raster.title , ".tif", sep="")
    rastercount <- paste("ensembles//count//", raster.title , ".tif", sep="")
    rasterpresence <- paste("ensembles//presence//", raster.title, ".tif", sep="")
#
    RASTER.species.orig <- RASTER.species.name
    if (RASTER.models.overwrite==T) {
        RASTER.species.name <- "working"
    }else{
        RASTER.species.name <- paste(RASTER.species.name, "_", stack.title, sep="")
    }
#
    cat(paste("\n", "Start of predictions for organism: ", RASTER.species.orig, "\n", sep = ""))
    cat(paste("Predictions for SpatRaster: ", stack.title, "\n", sep = ""))
    ensemble.statistics <- NULL
    cat(paste("ensemble raster layers will be saved in folder ", getwd(), "//ensembles", "\n\n", sep = ""))
    statistics.names <- c("n.models", "ensemble.threshold", "ensemble.min", "ensemble.max", "count.min", "count.max") 
    ensemble.statistics <- numeric(6)
    names(ensemble.statistics) <- statistics.names
#
# sometimes still error warnings for minimum and maximum values of the layers
# set minimum and maximum values for xn
    for (i in 1:terra::nlyr(xn)) {
        xn[[i]] <- terra::setMinMax(xn[[i]])
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
        fullname <- paste("models/", RASTER.species.name, "_MAXENT.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_MAXENT_step1.tif", sep="")
        if (CATCH.OFF == F) {
            tryCatch(pmaxent <- terra::predict(object=results, x=xn, na.rm=TRUE, 
                    filename=fullname, overwrite=TRUE),
                error= function(err) {print(paste("MAXENT prediction failed"))},
                silent=F)
            }else{
                pmaxent <- terra::predict(object=results, x=xn, na.rm=TRUE, 
                    filename=fullname, overwrite=TRUE)
            }
        if (is.null(pmaxent) == F) {
            results2 <- MAXENT.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste(fullname, "_step1", sep="")
                terra::writeRaster(x=pmaxent, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "MAXENT"
                pmaxent <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            pmaxent <- trunc(1000*pmaxent)
            terra::writeRaster(x=pmaxent, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(pmaxent, y=p)[,2]/1000
                abs1 <- terra::extract(pmaxent, y=a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["MAXENT"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["MAXENT"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(pmaxent, pt)[,2]/1000
                abs1 <- terra::extract(pmaxent, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: MAXENT prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["MAXENT"] <- -1
        }
    }
    if (output.weights["MAXNET"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Maximum entropy algorithm (package: maxnet)\n", sep=""))
        results <- MAXNET.OLD
        pmaxnet <- NULL
        fullname <- paste("models/", RASTER.species.name, "_MAXNET.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_MAXNET_step1.tif", sep="")
        if (CATCH.OFF == F) {
            tryCatch(pmaxnet <- terra::predict(object=xn, model=results, fun=predict.maxnet2, na.rm=TRUE, clamp=MAXNET.clamp, type=MAXNET.type,
                    filename=fullname, overwrite=TRUE),
                error= function(err) {print(paste("MAXNET prediction failed"))},
                silent=F)
            }else{
                pmaxnet <- terra::predict(object=xn, model=results, fun=predict.maxnet2, na.rm=TRUE, clamp=MAXNET.clamp, type=MAXNET.type,
                    filename=fullname, overwrite=TRUE)
            }
        if (is.null(pmaxnet) == F) {
            results2 <- MAXNET.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste(fullname, "_step1.tif", sep="")
                terra::writeRaster(x=pmaxnet, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "MAXNET"
                pmaxnet <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            pmaxnet <- trunc(1000*pmaxnet)
            terra::writeRaster(x=pmaxnet, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(pmaxnet, y=p)[, 2]/1000
                abs1 <- terra::extract(pmaxnet, y=a)[, 2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["MAXNET"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["MAXNET"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- (terra::extract(pmaxnet, pt)[,2])/1000
                abs1 <- (terra::extract(pmaxnet, at)[,2])/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: MAXNET prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["MAXNET"] <- -1
        }
    }
    if (output.weights["MAXLIKE"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Maxlike algorithm (package: maxlike)\n", sep=""))
        results <- MAXLIKE.OLD
        pmaxlike <- NULL
        fullname <- paste("models/", RASTER.species.name, "_MAXLIKE.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_MAXLIKE_step1.tif", sep="")
        xn.num <- terra::subset(xn, subset=models.list$num.vars)
        if (CATCH.OFF == F) {
            tryCatch(pmaxlike <- terra::predict(object=xn.num, model=results, na.rm=TRUE, 
                    filename=fullname, overwrite=TRUE),
                error= function(err) {print(paste("MAXLIKE prediction failed"))},
                silent=F)
        }else{
            pmaxlike <- terra::predict(object=xn.num, model=results, na.rm=TRUE, 
                filename=fullname, overwrite=TRUE)
        }
        if (is.null(pmaxlike) == F) {
            results2 <- MAXLIKE.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste("models//", "MAXLIKE_step1", sep="")
#                fullname2 <- paste(fullname, "_step1.tif", sep="")
                terra::writeRaster(x=pmaxlike, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "MAXLIKE"
                pmaxlike <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            pmaxlike <- trunc(1000*pmaxlike)
            terra::writeRaster(x=pmaxlike, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(pmaxlike, p)[,2]/1000
                abs1 <- terra::extract(pmaxlike, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["MAXLIKE"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["MAXLIKE"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(pmaxlike, pt)[,2]/1000
                abs1 <- terra::extract(pmaxlike, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: MAXLIKE prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["MAXLIKE"] <- -1
        }
    }
    if (output.weights["GBM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Generalized boosted regression modeling (package: gbm) \n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "WARNING: not certain whether the correct factor levels will be used", "\n", sep=""))
        } 
        results <- GBM.OLD
        pgbm <- NULL
        fullname <- paste("models/", RASTER.species.name, "_GBM.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_GBM_step1.tif", sep="")
        if (CATCH.OFF == F) {
            tryCatch(pgbm <- terra::predict(object=xn, model=results, na.rm=TRUE, factors=factlevels,
                    n.trees=results$n.trees, type="response", filename=fullname, overwrite=TRUE),
                error= function(err) {print(paste("GBM prediction failed"))},
                silent=F)
        }else{
            pgbm <- terra::predict(object=xn, model=results, na.rm=TRUE, factors=factlevels,
                n.trees=results$n.trees, type="response", filename=fullname, overwrite=TRUE)
        }
        if (is.null(pgbm) == F) {
            results2 <- GBM.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste(fullname, "_step1.tif", sep="")
                terra::writeRaster(x=pgbm, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "GBM"
                pgbm <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            pgbm <- trunc(1000*pgbm)
            terra::writeRaster(x=pgbm, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(pgbm, p)[,2]/1000
                abs1 <- terra::extract(pgbm, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["GBM"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["GBM"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(pgbm, pt)[,2]/1000
                abs1 <- terra::extract(pgbm, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
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
        fullname <- paste("models/", RASTER.species.name, "_GBMSTEP.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_GBMSTEP_step1.tif", sep="")
        if (CATCH.OFF == F) {
            tryCatch(pgbms <- terra::predict(object=xn, model=results, fun=gbm::predict.gbm, na.rm=TRUE, factors=factlevels,
                n.trees=results$n.trees, type="response", filename=fullname, overwrite=TRUE),
            error= function(err) {print(paste("stepwise GBM prediction failed"))},
            silent=F)
        }else{
            pgbms <- terra::predict(object=xn, model=results, fun=gbm::predict.gbm, na.rm=TRUE, factors=factlevels,
                n.trees=results$n.trees, type="response", filename=fullname, overwrite=TRUE)
        }
        if (is.null(pgbms) == F) {
            results2 <- GBMSTEP.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste(fullname, "_step1", sep="")
                terra::writeRaster(x=pgbms, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "GBMSTEP"
                pgbms <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            pgbms <- trunc(1000*pgbms)
# corrected writing in new format (August 2020)
            terra::writeRaster(x=pgbms, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(pgbms, p)[,2]/1000
                abs1 <- terra::extract(pgbms, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["GBMSTEP"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["GBMSTEP"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(pgbms, pt)[,2]/1000
                abs1 <- terra::extract(pgbms, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
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
        fullname <- paste("models/", RASTER.species.name, "_RF.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_RF_step1.tif", sep="")
        if (CATCH.OFF == F) {
            tryCatch(prf <- terra::predict(object=xn, model=results, fun=predict.RF, na.rm=TRUE, factors=factlevels,
                filename=fullname, overwrite=TRUE),
            error= function(err) {print(paste("random forest prediction failed"))},
            silent=F)
        }else{
            prf <- terra::predict(object=xn, model=results, fun=predict.RF, na.rm=TRUE, factors=factlevels,
                filename=fullname, overwrite=TRUE)
        }
        if (is.null(prf) == F) {
            results2 <- RF.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste(fullname, "_step1", sep="")
                terra::writeRaster(x=prf, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "RF"
                prf <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            prf <- trunc(1000*prf)
            terra::writeRaster(x=prf, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(prf, p)[,2]/1000
                abs1 <- terra::extract(prf, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["RF"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["RF"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(prf, pt)[,2]/1000
                abs1 <- terra::extract(prf, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: random forest prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["RF"] <- -1 
        }
    } 
    if (output.weights["CF"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Random forest algorithm (package: party)\n", sep=""))
        results <- CF.OLD
        pcf <- NULL
        fullname <- paste("models/", RASTER.species.name, "_CF.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_CF_step1.tif", sep="")
        if (CATCH.OFF == F) {
            tryCatch(pcf <- terra::predict(object=xn, model=results, fun=predict.CF, na.rm=TRUE, factors=factlevels,
                filename=fullname, overwrite=TRUE),
            error= function(err) {print(paste("random forest prediction failed"))},
            silent=F)
        }else{
            pcf <- terra::predict(object=xn, model=results, fun=predict.CF, na.rm=TRUE, factors=factlevels,
                filename=fullname, overwrite=TRUE)
        }
        if (is.null(pcf) == F) {
            results2 <- CF.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste(fullname, "_step1", sep="")
                terra::writeRaster(x=pcf, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "CF"
                pcf <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            pcf <- trunc(1000*pcf)
            terra::writeRaster(x=pcf, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(pcf, p)[,2]/1000
                abs1 <- terra::extract(pcf, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["CF"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["CF"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(pcf, pt)[,2]/1000
                abs1 <- terra::extract(pcf, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: random forest prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            output.weights["CF"] <- -1 
        }
    } 
    if (output.weights["GLM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Generalized Linear Model \n", sep=""))
        results <- GLM.OLD
        pglm <- NULL
        fullname <- paste("models/", RASTER.species.name, "_GLM.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_GLM_step1.tif", sep="")
        if (CATCH.OFF == T) {
            tryCatch(pglm <- terra::predict(object=xn, model=results, na.rm=TRUE, type="response", factors=factlevels,
                filename=fullname, overwrite=TRUE),
            error= function(err) {print(paste("GLM prediction failed"))},
            silent=F)
        }else{
            pglm <- terra::predict(object=xn, model=results, na.rm=TRUE, type="response", factors=factlevels,
                filename=fullname, overwrite=TRUE)
        }
        if (is.null(pglm) == F) {
            results2 <- GLM.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste(fullname, "_step1", sep="")
                terra::writeRaster(x=pglm, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "GLM"
                pglm <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            pglm <- trunc(1000*pglm)
            terra::writeRaster(x=pglm, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(pglm, p)[,2]/1000
                abs1 <- terra::extract(pglm, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["GLM"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["GLM"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(pglm, pt)[,2]/1000
                abs1 <- terra::extract(pglm, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
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
        fullname <- paste("models/", RASTER.species.name, "_GLMSTEP.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_GLMSTEP_step1.tif", sep="")
        if (CATCH.OFF == F) {
            tryCatch(pglms <- terra::predict(object=xn, model=results, na.rm=TRUE, type="response", factors=factlevels,
                filename=fullname, overwrite=TRUE),
            error= function(err) {print(paste("stepwise GLM prediction failed"))},
            silent=F)
        }else{
            pglms <- terra::predict(object=xn, model=results, na.rm=TRUE, type="response", factors=factlevels,
                filename=fullname, overwrite=TRUE)
        }
        if (is.null(pglms) == F) {
            results2 <- GLMSTEP.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste(fullname, "_step1", sep="")
                terra::writeRaster(x=pglms, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "GLMSTEP"
                pglms <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            pglms <- trunc(1000*pglms)
            terra::writeRaster(x=pglms, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(pglms, p)[,2]/1000
                abs1 <- terra::extract(pglms, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["GLMSTEP"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["GLMSTEP"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(pglms, pt)[,2]/1000
                abs1 <- terra::extract(pglms, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
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
        fullname <- paste("models/", RASTER.species.name, "_GAM.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_GAM_step1.tif", sep="")
        if (CATCH.OFF == F) {
            tryCatch(pgam <- terra::predict(object=xn, model=results, fun=gam::predict.Gam, na.rm=TRUE, type="response", factors=factlevels,
                filename=fullname, overwrite=TRUE),
            error= function(err) {print(paste("GAM (package: gam) prediction failed"))},
            silent=F)
        }else{
            pgam <- terra::predict(object=xn, model=results, fun=gam::predict.Gam, na.rm=TRUE, type="response", factors=factlevels,
                filename=fullname, overwrite=TRUE)
        }
        if (is.null(pgam) == F) {
            results2 <- GAM.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste(fullname, "_step1", sep="")
                terra::writeRaster(x=pgam, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "GAM"
                pgam <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            pgam <- trunc(1000*pgam)
            terra::writeRaster(x=pgam, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(pgam, p)[,2]/1000
                abs1 <- terra::extract(pgam, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["GAM"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["GAM"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(pgam, pt)[,2]/1000
                abs1 <- terra::extract(pgam, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
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
        fullname <- paste("models/", RASTER.species.name, "_GAMSTEP.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_GAMSTEP_step1.tif", sep="")
        if (CATCH.OFF == F) {
            tryCatch(pgams <- terra::predict(object=xn, model=results, fun=gam::predict.Gam, type="response", na.rm=TRUE, factors=factlevels,
                filename=fullname, overwrite=TRUE),
            error= function(err) {print(paste("stepwise GAM (package: gam) prediction failed"))},
            silent=F)
        }else{
            pgams <- terra::predict(object=xn, model=results, fun=gam::predict.Gam, type="response", na.rm=TRUE, factors=factlevels,
                filename=fullname, overwrite=TRUE)
        }
        if (is.null(pgams) == F) {
            results2 <- GAMSTEP.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste(fullname, "_step1", sep="")
                terra::writeRaster(x=pgams, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "GAMSTEP"
                pgams <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            pgams <- trunc(1000*pgams)
            terra::writeRaster(x=pgams, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(pgams, p)[,2]/1000
                abs1 <- terra::extract(pgams, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["GAMSTEP"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["GAMSTEP"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(pgams, pt)[,2]/1000
                abs1 <- terra::extract(pgams, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
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
        fullname <- paste("models/", RASTER.species.name, "_MGCV.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_MGCV_step1.tif", sep="")
        if (CATCH.OFF == F) {
            tryCatch(pmgcv <- terra::predict(object=xn, model=results, fun=predict.MGCV, na.rm=TRUE, type="response", factors=factlevels,
                filename=fullname, overwrite=TRUE),
            error= function(err) {print(paste("GAM (package: mgcv) prediction failed"))},
            silent=F)
        }else{
            pmgcv <- terra::predict(object=xn, model=results, fun=predict.MGCV, na.rm=TRUE, type="response", factors=factlevels,
                filename=fullname, overwrite=TRUE)
        }
        if (is.null(pmgcv) == F) {
            results2 <- MGCV.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste(fullname, "_step1", sep="")
                terra::writeRaster(x=pmgcv, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "MGCV"
                pmgcv <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            pmgcv <- trunc(1000*pmgcv)
            terra::writeRaster(x=pmgcv, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(pmgcv, p)[,2]/1000
                abs1 <- terra::extract(pmgcv, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["MGCV"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["MGCV"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(pmgcv, pt)[,2]/1000
                abs1 <- terra::extract(pmgcv, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
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
        fullname <- paste("models/", RASTER.species.name, "_MGCVFIX.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_MGCVFIX_step1.tif", sep="")
        if (CATCH.OFF == F) {
            tryCatch(pmgcvf <- terra::predict(object=xn, model=results, fun=predict.MGCV, na.rm=TRUE, type="response", factors=factlevels,
                filename=fullname, overwrite=TRUE),
            error= function(err) {print(paste("GAM with fixed d.f. regression splines (package: mgcv) prediction failed"))},
            silent=F)
        }else{
            pmgcvf <- terra::predict(object=xn, model=results, fun=predict.MGCV, na.rm=TRUE, type="response", factors=factlevels,
                filename=fullname, overwrite=TRUE)
        }
        if (is.null(pmgcvf) == F) {
            results2 <- MGCVFIX.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste(fullname, "_step1", sep="")
                terra::writeRaster(x=pmgcvf, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "MGCVFIX"
                pmgcvf <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            pmgcvf <- trunc(1000*pmgcvf)
            terra::writeRaster(x=pmgcvf, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(pmgcvf, p)[,2]/1000
                abs1 <- terra::extract(pmgcvf, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["MGCVFIX"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["MGCVFIX"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(pmgcvf, pt)[,2]/1000
                abs1 <- terra::extract(pmgcvf, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
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
        fullname <- paste("models/", RASTER.species.name, "_EARTH.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_EARTH_step1.tif", sep="")
        if (CATCH.OFF == F) {
            tryCatch(pearth <- terra::predict(object=xn, model=results, fun=predict.EARTH, na.rm=TRUE, type="response", factors=factlevels,
                filename=fullname, overwrite=TRUE),
            error= function(err) {print(paste("MARS (package: earth) prediction failed"))},
            silent=F)
        }else{
            pearth <- terra::predict(object=xn, model=results, fun=predict.EARTH, na.rm=TRUE, type="response", factors=factlevels,
                filename=fullname, overwrite=TRUE)
        }
        if (is.null(pearth) == F) {
            results2 <- EARTH.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste(fullname, "_step1", sep="")
                terra::writeRaster(x=pearth, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "EARTH"
                pearth <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            pearth <- trunc(1000*pearth)
            terra::writeRaster(x=pearth, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(pearth, p)[,2]/1000
                abs1 <- terra::extract(pearth, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["EARTH"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["EARTH"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(pearth, pt)[,2]/1000
                abs1 <- terra::extract(pearth, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
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
        fullname <- paste("models/", RASTER.species.name, "_RPART.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_RPART_step1.tif", sep="")
        if (CATCH.OFF == F) {
            tryCatch(prpart <- terra::predict(object=xn, model=results, na.rm=TRUE, type="prob", index=2, factors=factlevels,
                filename=fullname, overwrite=TRUE),
            error= function(err) {print(paste("RPART prediction failed"))},
            silent=F)
        }else{
            prpart <- terra::predict(object=xn, model=results, na.rm=TRUE, type="prob", index=2, factors=factlevels,
                filename=fullname, overwrite=TRUE)
        }
        if (is.null(prpart) == F) {
            results2 <- RPART.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste(fullname, "_step1", sep="")
                terra::writeRaster(x=prpart, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "RPART"
                prpart <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            prpart <- trunc(1000*prpart)
            terra::writeRaster(x=prpart, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(prpart, p)[,2]/1000
                abs1 <- terra::extract(prpart, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["RPART"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["RPART"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(prpart, pt)[,2]/1000
                abs1 <- terra::extract(prpart, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
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
        fullname <- paste("models/", RASTER.species.name, "_NNET.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_NNET_step1.tif", sep="")
        if (CATCH.OFF == F){
            tryCatch(pnnet <- terra::predict(object=xn, model=results, fun=predict.NNET, na.rm=TRUE, factors=factlevels,
                filename=fullname, overwrite=TRUE),
            error= function(err) {print(paste("Artificial Neural Network (package: nnet) prediction failed"))},
            silent=F)
        }else{
            pnnet <- terra::predict(object=xn, model=results, fun=predict.NNET, na.rm=TRUE, factors=factlevels,
                filename=fullname, overwrite=TRUE)
        }
        if (is.null(pnnet) == F) {
            results2 <- NNET.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste(fullname, "_step1", sep="")
                terra::writeRaster(x=pnnet, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "NNET"
                pnnet <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            pnnet <- trunc(1000*pnnet)
            terra::writeRaster(x=pnnet, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(pnnet, p)[,2]/1000
                abs1 <- terra::extract(pnnet, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["NNET"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["NNET"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(pnnet, pt)[,2]/1000
                abs1 <- terra::extract(pnnet, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
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
        fullname <- paste("models/", RASTER.species.name, "_FDA.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_FDA_step1.tif", sep="")
        if (CATCH.OFF == F) {
            tryCatch(pfda <- terra::predict(object=xn, model=results, na.rm=TRUE, type="posterior", index=2, factors=factlevels,
                filename=fullname, overwrite=TRUE),
            error= function(err) {print(paste("FDA prediction failed"))},
            silent=F)
        }else{
            pfda <- terra::predict(object=xn, model=results, na.rm=TRUE, type="posterior", index=2, factors=factlevels,
                filename=fullname, overwrite=TRUE)
        }
        if (is.null(pfda) == F) {
            results2 <- FDA.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste(fullname, "_step1", sep="")
                terra::writeRaster(x=pfda, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "FDA"
                pfda <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            pfda <- trunc(1000*pfda)
            terra::writeRaster(x=pfda, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(pfda, p)[,2]/1000
                abs1 <- terra::extract(pfda, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["FDA"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["FDA"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(pfda, pt)[,2]/1000
                abs1 <- terra::extract(pfda, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
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
        fullname <- paste("models/", RASTER.species.name, "_SVM.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_SVM_step1.tif", sep="")
        predict.svm2 <- as.function(kernlab::predict)
        if (CATCH.OFF == F) {
            tryCatch(psvm <- terra::predict(object=xn, model=results, fun=predict.svm2, na.rm=TRUE, type="probabilities", index=2, factors=factlevels,
                filename=fullname, overwrite=TRUE),
            error= function(err) {print(paste("Support Vector Machines (package: kernlab) prediction failed"))},
            silent=F)
        }else{
            psvm <- terra::predict(object=xn, model=results, fun=predict.svm2, na.rm=TRUE, type="probabilities", index=2, factors=factlevels,
                filename=fullname, overwrite=TRUE)
        }
        if (is.null(psvm) == F) {
            results2 <- SVM.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste(fullname, "_step1", sep="")
                terra::writeRaster(x=psvm, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "SVM"
                psvm <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            psvm <- trunc(1000*psvm)
            terra::writeRaster(x=psvm, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(psvm, p)[,2]/1000
                abs1 <- terra::extract(psvm, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["SVM"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["SVM"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(psvm, pt)[,2]/1000
                abs1 <- terra::extract(psvm, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
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
        fullname <- paste("models/", RASTER.species.name, "_SVME.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_SVME_step1.tif", sep="")
        if (CATCH.OFF == F) {
            tryCatch(psvme <- terra::predict(object=xn, model=results, fun=predict.SVME, na.rm=TRUE, factors=factlevels,
                filename=fullname, overwrite=TRUE),
            error= function(err) {print(paste("SVM prediction (e1071 package) failed"))},
            warning= function(war) {print(paste("SVM prediction (e1071 package) failed"))},
            silent=F)
        }else{
            psvme <- terra::predict(object=xn, model=results, fun=predict.SVME, na.rm=TRUE, factors=factlevels,
                filename=fullname, overwrite=TRUE)
        }
        if (is.null(psvme) == F) {
            results2 <- SVME.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste(fullname, "_step1", sep="")
                terra::writeRaster(x=psvme, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "SVME"
                psvme <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            psvme <- trunc(1000*psvme)
            terra::writeRaster(x=psvme, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(psvme, p)[,2]/1000
                abs1 <- terra::extract(psvme, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["SVME"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["SVME"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(psvme, pt)[,2]/1000
                abs1 <- terra::extract(psvme, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
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
        fullname <- paste("models/", RASTER.species.name, "_GLMNET.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_GLMNET_step1.tif", sep="")
        xn.num <- terra::subset(xn, subset=models.list$num.vars)
        if (CATCH.OFF == F) {
            tryCatch(pglmnet <- terra::predict(object=xn.num, model=results, fun=predict.GLMNET, na.rm=TRUE, GLMNET.class=GLMNET.class,
                filename=fullname, overwrite=TRUE),
            error= function(err) {print(paste("GLMNET prediction (glmnet package) failed"))},
            warning= function(war) {print(paste("GLMNET prediction (glmnet package) failed"))},
            silent=F)
        }else{
            pglmnet <- terra::predict(object=xn.num, model=results, fun=predict.GLMNET, na.rm=TRUE, GLMNET.class=GLMNET.class,
                filename=fullname, overwrite=TRUE)
        }
        if (is.null(pglmnet) == F) {
            results2 <- GLMNET.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste(fullname, "_step1", sep="")
                terra::writeRaster(x=pglmnet, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "GLMNET"
                pglmnet <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            pglmnet <- trunc(1000*pglmnet)
            terra::writeRaster(x=pglmnet, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(pglmnet, p)[,2]/1000
                abs1 <- terra::extract(pglmnet, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["GLMNET"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["GLMNET"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(pglmnet, pt)[,2]/1000
                abs1 <- terra::extract(pglmnet, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
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
        fullname <- paste("models/", RASTER.species.name, "_BIOCLIMO.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_BIOCLIMO_step1.tif", sep="")
        if (CATCH.OFF == F){
            tryCatch(pbioO <- terra::predict(object=xn, model=results, fun=predict.BIOCLIM.O, na.rm=TRUE, 
                filename=fullname, overwrite=TRUE),
            error= function(err) {print(paste("original BIOCLIM prediction failed"))},
            silent=F)
        }else{
            pbioO <- terra::predict(object=xn, model=results, fun=predict.BIOCLIM.O, na.rm=TRUE, 
                filename=fullname, overwrite=TRUE)
        }
        if (is.null(pbioO) == F) {
            results2 <- BIOCLIM.O.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste(fullname, "_step1", sep="")
                terra::writeRaster(x=pbioO, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "BIOCLIM.O"
                pbioO <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            pbioO <- trunc(1000*pbioO)
            terra::writeRaster(x=pbioO, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(pbioO, p)[,2]/1000
                abs1 <- terra::extract(pbioO, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["BIOCLIM.O"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["BIOCLIM.O"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(pbioO, pt)[,2]/1000
                abs1 <- terra::extract(pbioO, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
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
            xn <- terra::subset(xn, which((names(xn) %in% factors) == FALSE))
#            xn <- terra::rast(xn)              
        }
    }
    if (output.weights["BIOCLIM"] > 0) {  
        mc <- mc+1
        cat(paste("\n", mc, ". BIOCLIM algorithm (package: dismo)\n", sep=""))
        results <- BIOCLIM.OLD
        pbio <- NULL
        fullname <- paste("models/", RASTER.species.name, "_BIOCLIM.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_BIOCLIM_step1.tif", sep="")
        if (CATCH.OFF == F) {
#            tryCatch(pbio <- dismo::predict(object=results, x=xn, na.rm=TRUE, 
            tryCatch(pbio <- terra::predict(object=xn.num, model=results, na.rm=TRUE,
                filename=fullname, overwrite=TRUE),
            error= function(err) {print(paste("BIOCLIM prediction failed"))},
            silent=F)
        }else{
#            pbio <- dismo::predict(object=results, x=xn, na.rm=TRUE, 
#                filename=fullname, overwrite=TRUE)
            
            pbio <- terra::predict(object=xn.num, model=results, na.rm=TRUE,
                                      filename=fullname, overwrite=TRUE)            
            
        }
        if (is.null(pbio) == F) {
            results2 <- BIOCLIM.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste(fullname, "_step1", sep="")
                terra::writeRaster(x=pbio, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "BIOCLIM"
                pbio <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            pbio <- trunc(1000*pbio)
            terra::writeRaster(x=pbio, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(pbio, p)[,2]/1000
                abs1 <- terra::extract(pbio, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["BIOCLIM"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["BIOCLIM"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(pbio, pt)[,2]/1000
                abs1 <- terra::extract(pbio, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
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
            xn <- terra::subset(xn, which((names(xn) %in% dummy.vars.noDOMAIN) == FALSE)) 
#            xn <- terra::rast(xn)              
        }
    }
    if (output.weights["DOMAIN"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". DOMAIN algorithm (package: dismo)\n", sep=""))
        if(is.null(models.list$dummy.vars.noDOMAIN) == F) {
            xn <- terra::subset(xn, which((names(xn) %in% models.list$dummy.vars.noDOMAIN) == FALSE))
#            xn <- terra::rast(xn)          
        }
        results <- DOMAIN.OLD
        pdom <- NULL
        fullname <- paste("models/", RASTER.species.name, "_DOMAIN.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_DOMAIN_step1.tif", sep="")
        if (CATCH.OFF == F) {
#            tryCatch(pdom <- dismo::predict(object=results, x=xn, na.rm=TRUE, 
            tryCatch(pdom <- terra::predict(object=xn, model=results, na.rm=TRUE, 
                filename=fullname, overwrite=TRUE),
            error= function(err) {print(paste("DOMAIN prediction failed"))},
            silent=F)
        }else{
#            pdom <- dismo::predict(object=results, x=xn, na.rm=TRUE, 
#                filename=fullname, overwrite=TRUE)
            
            pdom <- terra::predict(object=xn, model=results, na.rm=TRUE,
                                   filename=fullname, overwrite=TRUE) 
            
        }
        if (is.null(pdom) == F) {
            results2 <- DOMAIN.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
#                fullname2 <- paste(fullname, "_step1", sep="")
                terra::writeRaster(x=pdom, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "DOMAIN"
                pdom <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            pdom <- trunc(1000*pdom)
            terra::writeRaster(x=pdom, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(pdom, p)[,2]/1000
                abs1 <- terra::extract(pdom, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["DOMAIN"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["DOMAIN"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(pdom, pt)[,2]/1000
                abs1 <- terra::extract(pdom, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
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
            xn <- terra::subset(xn, which((names(xn) %in% dummy.vars) == FALSE))
#            xn <- terra::rast(xn)            
        }
    }
    if (output.weights["MAHAL"] > 0) {  
        mc <- mc+1
        cat(paste("\n", mc, ". Mahalanobis algorithm (package: dismo)\n", sep=""))
        results <- MAHAL.OLD
        pmahal <- NULL
        fullname <- paste("models/", RASTER.species.name, "_MAHAL.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_MAHAL_step1.tif", sep="")
# not possible to use the predict.mahal function as terra::predict automatically reverts to dismo::predict for 'DistModel' objects
        results2 <- MAHAL.PROBIT.OLD
# PROBIT FALSE        
        if (is.null(results2) == F) {
            if (CATCH.OFF == F) {
                tryCatch(pmahal <- terra::predict(object=xn, model=results, fun=predict.MAHAL, na.rm=TRUE, PROBIT=FALSE,
                    filename=fullname, overwrite=TRUE),
                error= function(err) {print(paste("Mahalanobis prediction failed"))},
                silent=F)
            }else{       
                pmahal <- terra::predict(object=xn, model=results, fun=predict.MAHAL, na.rm=TRUE, PROBIT=FALSE,
                    filename=fullname, overwrite=TRUE)
            }
# PROBIT TRUE
        }else{
            if (CATCH.OFF == F) {
                tryCatch(pmahal <- terra::predict(object=xn, model=results, fun=predict.MAHAL, na.rm=TRUE, PROBIT=TRUE,
                    filename=fullname, overwrite=TRUE),
                error= function(err) {print(paste("Mahalanobis prediction failed"))},
                silent=F)
            }else{
                pmahal <- terra::predict(object=xn, model=results, fun=predict.MAHAL, na.rm=TRUE, PROBIT=TRUE,
                    filename=fullname, overwrite=TRUE)
            }
        }

        if (is.null(pmahal) == F) {
#            results2 <- MAHAL.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
 #               fullname2 <- paste(fullname, "_step1", sep="")
                terra::writeRaster(x=pmahal, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "MAHAL"
                pmahal <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            pmahal <- trunc(1000*pmahal)
            terra::writeRaster(x=pmahal, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(pmahal, p)[,2]/1000
                abs1 <- terra::extract(pmahal, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["MAHAL"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["MAHAL"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(pmahal, pt)[,2]/1000
                abs1 <- terra::extract(pmahal, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
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
            xn <- terra::subset(xn, which((names(xn) %in% dummy.vars) == FALSE))
#            xn <- terra::rast(xn)            
        }
        results <- MAHAL01.OLD
#        pmahal <- NULL
        fullname <- paste("models/", RASTER.species.name, "_MAHAL01.tif", sep="")
        fullname2 <- paste("models/", RASTER.species.name, "_MAHAL01_step1.tif", sep="")
# not possible to use the predict.mahal function as terra::predict automatically reverts to dismo::predict for 'DistModel' objects
        if (CATCH.OFF == F) {
            tryCatch(pmahal01 <- terra::predict(object=xn, model=results, fun=predict.MAHAL01, na.rm=TRUE, MAHAL.shape=MAHAL.shape,
                filename=fullname, overwrite=TRUE),
            error= function(err) {print(paste("transformed Mahalanobis prediction failed"))},
            silent=F)
        }else{
            pmahal01 <- terra::predict(object=xn, model=results, fun=predict.MAHAL01, na.rm=TRUE, MAHAL.shape=MAHAL.shape,
                filename=fullname, overwrite=TRUE)
        }
        if (is.null(pmahal01) == F) {
#            pmahal <- pmahal - 1 - MAHAL.shape
#            pmahal <- abs(pmahal)
#            pmahal <- MAHAL.shape / pmahal
            results2 <- MAHAL01.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
 #               fullname2 <- paste(fullname, "_step1", sep="")
                terra::writeRaster(x=pmahal01, filename=fullname2, overwrite=TRUE)
                explan.stack <- terra::rast(fullname2)
                names(explan.stack) <- "MAHAL01"
                pmahal01 <- terra::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, overwrite=TRUE)                
            }
            pmahal01 <- trunc(1000*pmahal01)
            terra::writeRaster(x=pmahal01, filename=fullname, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- terra::extract(pmahal01, p)[,2]/1000
                abs1 <- terra::extract(pmahal01, a)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
                print(eval1)
                thresholds.raster["MAHAL01"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
                    threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres1, Abs=abs1)
                cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
                print(as.numeric(thresholds.raster["MAHAL01"]))
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- terra::extract(pmahal01, pt)[,2]/1000
                abs1 <- terra::extract(pmahal01, at)[,2]/1000
                eval1 <- dismo::evaluate(p=pres1, a=abs1)
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
#    ensemble <- xn[[1]] == terra::NAflag(xn[[1]])
# avoid problems with SRS not matching
    ensemble <- terra::rast(fullname)
    ensemble <- ensemble > Inf
    terra::setMinMax(ensemble)
    names(ensemble) <- raster.title
    terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
    enscount <- ensemble 
    terra::setMinMax(enscount)
    names(enscount) <- paste(raster.title, "_count", sep="")
    terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    enspresence <- ensemble
    terra::setMinMax(enspresence)
    names(enspresence) <- paste(raster.title, "_presence", sep="")
    terra::writeRaster(x=enspresence, filename=rasterpresence, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    if (output.weights["MAXENT"] > 0) {
        ensemble <- ensemble + output.weights["MAXENT"] * pmaxent
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)       
        pmaxent <- pmaxent >= 1000 * thresholds["MAXENT"]
        enscount <- enscount + pmaxent
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["MAXNET"] > 0) {
        ensemble <- ensemble + output.weights["MAXNET"] * pmaxnet
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)       
        pmaxnet <- pmaxnet >= 1000 * thresholds["MAXNET"]
        enscount <- enscount + pmaxnet
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["MAXLIKE"] > 0) {
        ensemble <- ensemble + output.weights["MAXLIKE"] * pmaxlike
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)       
        pmaxlike <- pmaxlike >= 1000 * thresholds["MAXLIKE"]
        enscount <- enscount + pmaxlike
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["GBM"] > 0) {
        ensemble <- ensemble + output.weights["GBM"] * pgbm
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)      
        pgbm <- pgbm >= 1000 * thresholds["GBM"]
        enscount <- enscount + pgbm
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["GBMSTEP"] > 0) {
        ensemble <- ensemble + output.weights["GBMSTEP"] * pgbms
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)      
        pgbms <- pgbms >= 1000 * thresholds["GBMSTEP"]
        enscount <- enscount + pgbms
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["RF"] > 0) {
        ensemble <- ensemble + output.weights["RF"] * prf
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)      
        prf <- prf >= 1000 * thresholds["RF"]
        enscount <- enscount + prf
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["CF"] > 0) {
        ensemble <- ensemble + output.weights["CF"] * pcf
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)      
        pcf <- pcf >= 1000 * thresholds["CF"]
        enscount <- enscount + pcf
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["GLM"] > 0) {
        ensemble <- ensemble + output.weights["GLM"] * pglm
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)    
        pglm <- pglm >= 1000 * thresholds["GLM"]
        enscount <- enscount + pglm
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["GLMSTEP"] > 0) {
        ensemble <- ensemble + output.weights["GLMSTEP"] * pglms
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)      
        pglms <- pglms >= 1000 * thresholds["GLMSTEP"]
        enscount <- enscount + pglms
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["GAM"] > 0) {
        ensemble <- ensemble + output.weights["GAM"] * pgam
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pgam <- pgam >= 1000 * thresholds["GAM"]
        enscount <- enscount + pgam
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["GAMSTEP"] > 0) {
        ensemble <- ensemble + output.weights["GAMSTEP"] * pgams
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pgams <- pgams >= 1000 * thresholds["GAMSTEP"]
        enscount <- enscount + pgams
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["MGCV"] > 0) {
        ensemble <- ensemble + output.weights["MGCV"] * pmgcv
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pmgcv <- pmgcv >= 1000 * thresholds["MGCV"]
        enscount <- enscount + pmgcv
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["MGCVFIX"] > 0) {
        ensemble <- ensemble + output.weights["MGCVFIX"] * pmgcvf
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pmgcvf <- pmgcvf >= 1000 * thresholds["MGCVFIX"]
        enscount <- enscount + pmgcvf
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["EARTH"] > 0) {
        ensemble <- ensemble + output.weights["EARTH"] * pearth
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pearth <- pearth >= 1000 * thresholds["EARTH"]
        enscount <- enscount + pearth
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["RPART"] > 0) {
        ensemble <- ensemble + output.weights["RPART"] * prpart
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        prpart <- prpart >= 1000 * thresholds["RPART"]
        enscount <- enscount + prpart
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["NNET"] > 0) {
        ensemble <- ensemble + output.weights["NNET"] * pnnet
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pnnet <- pnnet >= 1000 * thresholds["NNET"]
        enscount <- enscount + pnnet
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["FDA"] > 0) {
        ensemble <- ensemble + output.weights["FDA"] * pfda
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pfda <- pfda >= 1000 * thresholds["FDA"]
        enscount <- enscount + pfda
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["SVM"] > 0) {
        ensemble <- ensemble + output.weights["SVM"] * psvm
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        psvm <- psvm >= 1000 * thresholds["SVM"]
        enscount <- enscount + psvm
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["SVME"] > 0) {
        ensemble <- ensemble + output.weights["SVME"] * psvme
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        psvme <- psvme >= 1000 * thresholds["SVME"]
        enscount <- enscount + psvme
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["GLMNET"] > 0) {
        ensemble <- ensemble + output.weights["GLMNET"] * pglmnet
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pglmnet <- pglmnet >= 1000 * thresholds["GLMNET"]
        enscount <- enscount + pglmnet
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["BIOCLIM.O"] > 0) {
        ensemble <- ensemble + output.weights["BIOCLIM.O"] * pbioO
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pbioO <- pbioO >= 1000 * thresholds["BIOCLIM.O"]
        enscount <- enscount + pbioO
 #       terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["BIOCLIM"] > 0) {
        ensemble <- ensemble + output.weights["BIOCLIM"] * pbio
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pbio <- pbio >= 1000 * thresholds["BIOCLIM"]
        enscount <- enscount + pbio
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["DOMAIN"] > 0) {
        ensemble <- ensemble + output.weights["DOMAIN"] * pdom
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pdom <- pdom >= 1000 * thresholds["DOMAIN"]
        enscount <- enscount + pdom
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["MAHAL"] > 0) {
        ensemble <- ensemble + output.weights["MAHAL"] * pmahal
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pmahal <- pmahal >= 1000 * thresholds["MAHAL"]
        enscount <- enscount + pmahal
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
    if (output.weights["MAHAL01"] > 0) {
        ensemble <- ensemble + output.weights["MAHAL01"] * pmahal01
#        terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pmahal01 <- pmahal01 >= 1000 * thresholds["MAHAL01"]
        enscount <- enscount + pmahal01
#        terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    }
#
terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
    
# note that submodels had already been multiplied by 1000
    ensemble <- trunc(ensemble)
    terra::setMinMax(ensemble)
    ensemble.statistics[c("ensemble.min", "ensemble.max")] <- as.numeric(terra::minmax(ensemble))
#    names(ensemble) <- raster.title
#    terra::writeRaster(x=ensemble, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
#  avoid possible problems with saving of names of the raster layers
    terra::writeRaster(ensemble, filename="working.tif", overwrite=T)
    working.raster <- terra::rast("working.tif")
    names(working.raster) <- raster.title
    terra::writeRaster(working.raster, filename=rasterfull, overwrite=TRUE, filetype=RASTER.filetype, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
#
    terra::setMinMax(enscount)
    ensemble.statistics[c("count.min", "count.max")] <- as.numeric(terra::minmax(enscount))
#    names(enscount) <- paste(raster.title, "_count", sep="")
#    terra::writeRaster(x=enscount, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
#  avoid possible problems with saving of names of the raster layers
    terra::writeRaster(enscount, filename="working.tif", overwrite=T)
    working.raster <- terra::rast("working.tif")
    names(working.raster) <- paste(raster.title, "_count", sep="")
    terra::writeRaster(working.raster, filename=rastercount, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
#
    if(evaluate == T) {
        eval1 <- NULL
        cat(paste("\n", "Evaluation of created ensemble raster layer (", rasterfull, ") at locations p and a", "\n\n", sep = ""))
        pres_consensus <- terra::extract(ensemble, p)[,2]/1000
        abs_consensus <- terra::extract(ensemble, a)[,2]/1000
        eval1 <- dismo::evaluate(p=pres_consensus, a=abs_consensus)
        print(eval1)
        thresholds["ENSEMBLE"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
            threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres_consensus, Abs=abs_consensus)
        cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
        print(as.numeric(thresholds["ENSEMBLE"]))
    }
    if(retest == T) {
        eval1 <- NULL
        cat(paste("\n", "Evaluation of created ensemble raster layer (", rasterfull, ") at locations pt and at", "\n\n", sep = ""))
        pres_consensus <- terra::extract(ensemble, pt)[,2]/1000
        abs_consensus <- terra::extract(ensemble, at)[,2]/1000
        eval1 <- dismo::evaluate(p=pres_consensus, a=abs_consensus)
        print(eval1)
    }
    ensemble.statistics["ensemble.threshold"] <- thresholds["ENSEMBLE"]
#
    enspresence <- ensemble >= 1000 * thresholds["ENSEMBLE"]
    terra::setMinMax(enspresence)
#    names(enspresence) <- paste(raster.title, "_presence", sep="")
#    terra::writeRaster(x=enspresence, filename=rasterpresence, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
#  avoid possible problems with saving of names of the raster layers
    terra::writeRaster(enspresence, filename="working.tif", overwrite=T)
    working.raster <- terra::rast("working.tif")
    names(working.raster) <- paste(raster.title, "_presence", sep="")
    terra::writeRaster(working.raster, filename=rasterpresence, overwrite=TRUE, filetype=RASTER.filetype, datatype="INT1U", NAflag=255)
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

