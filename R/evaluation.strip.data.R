`evaluation.strip.data` <- function(
    xn=NULL, ext=NULL,
    models.list=NULL, 
    input.weights=models.list$output.weights,
    steps=200, CATCH.OFF=FALSE
)
{
    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if (is.null(xn) == T) {stop("value for parameter xn is missing (RasterStack object)")}
    if (is.null(models.list) == T) {stop("provide 'models.list' as models will not be recalibrated and retested")}
    if (is.null(input.weights) == T) {input.weights <- models.list$output.weights}
    if (is.null(ext) == F) {
        if(length(xn@title) == 0) {xn@title <- "stack1"}
        title.old <- xn@title
        xn <- raster::crop(xn, y=ext, snap="in")
        xn <- raster::stack(xn)
        xn@title <- title.old
    }
#
# check if all variables are present
    vars <- models.list$vars
    vars.xn <- names(xn)
    nv <- length(vars) 
    for (i in 1:nv) {
        if (any(vars.xn==vars[i]) == F) {stop("explanatory variable '", vars[i], "' not among grid layers of RasterStack xn", "\n", sep = "")}
    }
    nv <- length(vars.xn) 
    for (i in 1:nv) {
        if (any(vars==vars.xn[i]) == F) {
            cat(paste("\n", "NOTE: RasterStack layer '", vars.xn[i], "' was not calibrated as explanatory variable", "\n", sep = ""))
            xn <- raster::dropLayer(xn, which(names(xn) %in% c(vars.xn[i]) ))
            xn <- raster::stack(xn)
        }
    }
    factors <- models.list$factors
    dummy.vars <- models.list$dummy.vars
    dummy.vars.noDOMAIN <- models.list$dummy.vars.noDOMAIN

#
# set minimum and maximum values for xn
    for (i in 1:raster::nlayers(xn)) {
        xn[[i]] <- raster::setMinMax(xn[[i]])
    }
# declare categorical layers for xn
#    factors <- models.list$factors
    if(is.null(factors) == F) {
        for (i in 1:length(factors)) {
            j <- which(names(xn) == factors[i])
            xn[[j]] <- raster::as.factor(xn[[j]])
        }
    }
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
        if (is.null(models.list$MAXLIKE) == F) {MAXLIKE.OLD <- models.list$MAXLIKE}
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
    if (FDA > 0) {
        if (! requireNamespace("mda")) {stop("Please install the mda package")}
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
    ws <- input.weights

#
# prepare data set
    vars.xn <- names(xn)
    nvars <- length(vars)
    nnum <- nvars - length(factors)
    nrows <- nnum * steps
    plot.data <- array(dim=c(nnum*steps, nvars+2), NA)
    dimnames(plot.data)[[2]] <- c("focal.var", "categorical", vars)
# for categorical variables first
    fixedlevel <- array(dim=c(nvars))
    for (i in 1:nvars) {
        if(any(vars[i] == factors) == T) {
            il <- which(vars.xn == vars[i])
            tabulation <- data.frame(raster::freq(xn[[il]]))
            NA.index <- !is.na(tabulation[,"value"])
            tabulation <- tabulation[NA.index,]           
            plot.data2 <- array(dim=c(nrow(tabulation), nvars+2), NA)
            plot.data2[, 1] <- rep(i, nrow(tabulation))
            plot.data2[, 2] <- rep(1, nrow(tabulation))
            plot.data2[, i+2] <- tabulation[,1]
            plot.data <- rbind(plot.data, plot.data2)
            fixedlevel[i] <- tabulation[which.max(tabulation[,2]),1]
        }
    }
    for (i in 1:nvars) {
        if(any(vars[i] == factors) == T) {
            index <- is.na(plot.data[,i+2])
            plot.data[index,i+2] <- fixedlevel[i]
        }
    }
    nrows <- nrow(plot.data)
# for numerical variables next
    for (i in 1:nvars) {
        if(any(vars[i] == factors) == F) {
            il <- which(vars.xn == vars[i]) 
            plot.data[,i+2] <- rep(raster::cellStats(xn[[il]], stat="mean"), nrows)
        }
    }
    j <- 0
    for (i in 1:nvars) {
        if(any(vars[i] == factors) == F) {
            il <- which(vars.xn == vars[i]) 
            j <- j+1
            startpos <- (j-1)*steps+1
            endpos <- (j-1)*steps+steps
            plot.data[startpos:endpos,1] <- rep(i, steps)
            plot.data[startpos:endpos,2] <- rep(0, steps)
            minv <- raster::minValue(xn[[il]])
            maxv <- raster::maxValue(xn[[il]])
            plot.data[startpos:endpos,i+2] <- seq(from=minv,to=maxv, length.out=steps)
         }
    }
# declare factor variables
    plot.data.vars <- data.frame(plot.data)
    plot.data.vars <- plot.data.vars[, which(names(plot.data.vars) %in% vars), drop=F]
    plot.data.numvars <- plot.data.vars

    for (i in 1:nvars) {
        if(any(vars[i] == factors) == T) {
            plot.data.vars[,i+2] <- factor(plot.data.vars[,i+2], levels=models.list$categories[[vars[i]]])
            plot.data.numvars <- plot.data.numvars[, which(names(plot.data.numvars) != vars[i]), drop=F]
        }
    }

    plot.data.domain <- plot.data.numvars
    for (i in 1:nvars) {
        if(any(vars[i] == dummy.vars) == T) {
            plot.data.domain <- plot.data.domain[, which(names(plot.data.domain) != dummy.vars[i]), drop=F]
        }
    }

    plot.data.mahal <- plot.data.numvars
    for (i in 1:nvars) {
        if(any(vars[i] == dummy.vars) == T) {
            plot.data.mahal <- plot.data.mahal[, which(names(plot.data.mahal) != dummy.vars[i]), drop=F]
        }
    }

    assign("plot.data.vars", plot.data.vars, envir=.BiodiversityR)
    assign("plot.data.numvars", plot.data.numvars, envir=.BiodiversityR)
    assign("plot.data.domain", plot.data.domain, envir=.BiodiversityR)
    assign("plot.data.mahal", plot.data.mahal, envir=.BiodiversityR)

    modelnames <- c("MAXENT", "MAXLIKE", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", 
        "MGCVFIX", "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "GLMNET", 
        "BIOCLIM.O", "BIOCLIM", "DOMAIN", "MAHAL", "MAHAL01", "ENSEMBLE")
    nmodels <- length(modelnames)
    modelout <- array(dim=c(nrows, nmodels), 0.0)
    dimnames(modelout)[[2]] <- modelnames
    plot.data <- cbind(plot.data, modelout)
    plot.data <- data.frame(plot.data)
    for (i in 1:nvars) {
        if(any(vars[i] == factors) == T) {
            plot.data[,i+2] <- factor(plot.data[,i+2], levels=models.list$categories[[vars[i]]])
        }
    }
#
# sometimes still error warnings for minimum and maximum values of the layers
# set minimum and maximum values for xn
    for (i in 1:raster::nlayers(xn)) {
        xn[[i]] <- raster::setMinMax(xn[[i]])
    }
#
# Different modelling algorithms
#
    if (MAXENT > 0) {
        results <- MAXENT.OLD
        if (CATCH.OFF == F) {
            tryCatch(plot.data[,"MAXENT"] <- dismo::predict(object=results, x=plot.data.vars),
                error= function(err) {print(paste("MAXENT prediction failed"))},
                silent=F)
        }else{
            plot.data[,"MAXENT"] <- dismo::predict(object=results, x=plot.data.vars)
        }
        results2 <- MAXENT.PROBIT.OLD
        if (is.null(results2) == F) {
            plot.data[,"MAXENT.step1"] <- plot.data[, "MAXENT"]
            plot.data[,"MAXENT"] <- predict.glm(object=results2, newdata=plot.data, type="response")
        }
    }
    if (MAXLIKE > 0) {
        results <- MAXLIKE.OLD
        if (CATCH.OFF == F) {
            tryCatch(plot.data[,"MAXLIKE"] <- predict(object=results, newdata=plot.data.vars),
                error= function(err) {print(paste("MAXLIKE prediction failed"))},
                silent=F)
        }else{
            plot.data[,"MAXLIKE"] <- predict(object=results, newdata=plot.data.vars)
        }
        results2 <- MAXLIKE.PROBIT.OLD
        if (is.null(results2) == F) {
            plot.data[, "MAXLIKE.step1"] <- plot.data[, "MAXLIKE"]
            plot.data[,"MAXLIKE"] <- predict.glm(object=results2, newdata=plot.data, type="response")
        }
    }
    if (GBM > 0) {
        results <- GBM.OLD
        if (CATCH.OFF == F) {
            tryCatch(plot.data[,"GBM"] <- gbm::predict.gbm(object=results, newdata=plot.data.vars, n.trees=results$n.trees, type="response"),
                error= function(err) {print(paste("GBM prediction failed"))},
                silent=F)
        }else{
            plot.data[,"GBM"] <- gbm::predict.gbm(object=results, newdata=plot.data.vars, n.trees=results$n.trees, type="response")
        }
        results2 <- GBM.PROBIT.OLD
        if (is.null(results2) == F) {
            plot.data[, "GBM.step1"] <- plot.data[, "GBM"]
            plot.data[,"GBM"] <- predict.glm(object=results2, newdata=plot.data, type="response")
        }
    }
    if (GBMSTEP > 0) {
        results <- GBMSTEP.OLD
        if (CATCH.OFF == F) {
            tryCatch(plot.data[,"GBMSTEP"] <- gbm::predict.gbm(object=results, newdata=plot.data.vars, n.trees=results$n.trees, type="response"),
                error= function(err) {print(paste("stepwise GBM prediction failed"))},
                silent=F)
        }else{
            plot.data[,"GBMSTEP"] <- gbm::predict.gbm(object=results, newdata=plot.data.vars, n.trees=results$n.trees, type="response")
        }
        results2 <- GBMSTEP.PROBIT.OLD
        if (is.null(results2) == F) {
            plot.data[,"GBMSTEP.step1"] <- plot.data[, "GBMSTEP"]
            plot.data[,"GBMSTEP"] <- predict.glm(object=results2, newdata=plot.data, type="response")
        }
    }
    if (RF > 0) {
        results <- RF.OLD
        if (CATCH.OFF == F) {
            tryCatch(plot.data[,"RF"] <- predict.RF(object=results, newdata=plot.data.vars),
                error= function(err) {print(paste("RF prediction failed"))},
                silent=F)
        }else{
            plot.data[,"RF"] <- predict.RF(object=results, newdata=plot.data.vars)
        }
        results2 <- RF.PROBIT.OLD
        if (is.null(results2) == F) {
            plot.data[,"RF.step1"] <- plot.data[,"RF"]
            plot.data[,"RF"] <- predict.glm(object=results2, newdata=plot.data, type="response")
        }
    }
    if (GLM > 0) {
        results <- GLM.OLD
        if (CATCH.OFF == F) {
            tryCatch(plot.data[,"GLM"] <- predict.glm(object=results, newdata=plot.data.vars, type="response"),
                error= function(err) {print(paste("GLM prediction failed"))},
                silent=F)
        }else{
            plot.data[,"GLM"] <- predict.glm(object=results, newdata=plot.data.vars, type="response")
        }
        results2 <- GLM.PROBIT.OLD
        if (is.null(results2) == F) {
            plot.data[,"GLM.step1"] <- plot.data[, "GLM"]
            plot.data[,"GLM"] <- predict.glm(object=results2, newdata=plot.data, type="response")
        }
    }
    if (GLMSTEP > 0) {
        results <- GLMSTEP.OLD
        if (CATCH.OFF == F) {
            tryCatch(plot.data[,"GLMSTEP"] <- predict.glm(object=results, newdata=plot.data.vars, type="response"),
                error= function(err) {print(paste("stepwise GLM prediction failed"))},
                silent=F)
        }else{
            plot.data[,"GLMSTEP"] <- predict.glm(object=results, newdata=plot.data.vars, type="response")
        }
        results2 <- GLMSTEP.PROBIT.OLD
        if (is.null(results2) == F) {
            plot.data[,"GLMSTEP.step1"] <- plot.data[, "GLMSTEP"]
            plot.data[,"GLMSTEP"] <- predict.glm(object=results2, newdata=plot.data, type="response")
        }
    }
    if (GAM > 0) {
        results <- GAM.OLD
        if (CATCH.OFF == F) {
            tryCatch(plot.data[,"GAM"] <- gam::predict.Gam(object=results, newdata=plot.data.vars, type="response"),
                error= function(err) {print(paste("GAM (package: gam) prediction failed"))},
                silent=F)
        }else{
            plot.data[,"GAM"] <- gam::predict.Gam(object=results, newdata=plot.data.vars, type="response")
        }
        results2 <- GAM.PROBIT.OLD
        if (is.null(results2) == F) {
            plot.data[,"GAM.step1"] <- plot.data[, "GAM"]
            plot.data[,"GAM"] <- predict.glm(object=results2, newdata=plot.data, type="response")
        }
    }
    if (GAMSTEP > 0) {
        results <- GAMSTEP.OLD
        if (CATCH.OFF == F) {
            tryCatch(plot.data[,"GAMSTEP"] <- gam::predict.Gam(object=results, newdata=plot.data.vars, type="response"),
                error= function(err) {print(paste("stepwise GAM prediction (gam package) failed"))},
                silent=F)
        }else{
            plot.data[,"GAMSTEP"] <- gam::predict.Gam(object=results, newdata=plot.data.vars, type="response")
        }
        results2 <- GAMSTEP.PROBIT.OLD
        if (is.null(results2) == F) {
            plot.data[, "GAMSTEP.step1"] <- plot.data[, "GAMSTEP"]
            plot.data[,"GAMSTEP"] <- predict.glm(object=results2, newdata=plot.data, type="response")
        }
    }
    if (MGCV > 0) {
        results <- MGCV.OLD
        if (CATCH.OFF == F) {
            tryCatch(plot.data[,"MGCV"] <- predict.MGCV(object=results, newdata=plot.data.vars),
                error= function(err) {print(paste("GAM prediction (mgcv package) failed"))},
                silent=F)
        }else{
            plot.data[,"MGCV"] <- predict.MGCV(object=results, newdata=plot.data.vars)
        }
        results2 <- MGCV.PROBIT.OLD
        if (is.null(results2) == F) {
            plot.data[,"MGCV.step1"] <- plot.data[,"MGCV"]
            plot.data[,"MGCV"] <- predict.glm(object=results2, newdata=plot.data, type="response")
        }
    }
    if (MGCVFIX > 0) {
        results <- MGCVFIX.OLD
        if (CATCH.OFF == F) {
            tryCatch(plot.data[,"MGCVFIX"] <- predict.MGCV(object=results, newdata=plot.data.vars),
                error= function(err) {print(paste("MGCVFIX prediction (mgcv package) failed"))},
                silent=F)
        }else{
            plot.data[,"MGCVFIX"] <- predict.MGCV(object=results, newdata=plot.data.vars)
        }
        results2 <- MGCVFIX.PROBIT.OLD
        if (is.null(results2) == F) {plot.data[,"MGCVFIX"] <- predict.glm(object=results2, newdata=plot.data, type="response")}
    }
    if (EARTH > 0) {
        results <- EARTH.OLD
        if (CATCH.OFF == F) {
            tryCatch(plot.data[,"EARTH"] <- predict.EARTH(object=results, newdata=plot.data.numvars),
                error= function(err) {print(paste("MARS prediction (earth package) failed"))},
                silent=F)
        }else{
            plot.data[,"EARTH"] <- predict.EARTH(object=results, newdata=plot.data.numvars)
        }
        results2 <- EARTH.PROBIT.OLD
        if (is.null(results2) == F) {
            plot.data[,"EARTH.step1"] <- plot.data[,"EARTH"]
            plot.data[,"EARTH"] <- predict.glm(object=results2, newdata=plot.data, type="response")
        }
    }
    if (RPART > 0) {
        results <- RPART.OLD
        if (CATCH.OFF == F) {
            tryCatch(plot.data[,"RPART"] <- predict(object=results, newdata=plot.data.vars, type="prob")[,2],
                error= function(err) {print(paste("RPART prediction failed"))},
                silent=F)
        }else{
            plot.data[,"RPART"] <- predict(object=results, newdata=plot.data.vars, type="prob")[,2]
        }
        results2 <- RPART.PROBIT.OLD
        if (is.null(results2) == F) {
            plot.data[,"RPART.step1"] <- plot.data[,"RPART"]
            plot.data[,"RPART"] <- predict.glm(object=results2, newdata=plot.data, type="response")
        }
    }
    if (NNET > 0) {
        results <- NNET.OLD
        if (CATCH.OFF == F) {
            tryCatch(plot.data[,"NNET"] <- predict.NNET(object=results, newdata=plot.data.vars),
                error= function(err) {print(paste("ANN prediction (nnet package) failed"))},
                silent=F)
        }else{
            plot.data[,"NNET"] <- predict.NNET(object=results, newdata=plot.data.vars)
        }
        results2 <- NNET.PROBIT.OLD
        if (is.null(results2) == F) {
            plot.data[,"NNET.step1"] <- plot.data[,"NNET"]
            plot.data[,"NNET"] <- predict.glm(object=results2, newdata=plot.data, type="response")
        }
    }
    if (FDA > 0) {
        results <- FDA.OLD
        if (CATCH.OFF == F) {
            tryCatch(plot.data[,"FDA"] <- predict(object=results, newdata=plot.data.vars, type="posterior")[,2],
                error= function(err) {print(paste("FDA prediction failed"))},
                silent=F)
        }else{
            plot.data[,"FDA"] <- predict(object=results, newdata=plot.data.vars, type="posterior")[,2]
        }
        results2 <- FDA.PROBIT.OLD
        if (is.null(results2) == F) {
            plot.data[,"FDA.step1"] <- plot.data[,"FDA"]
            plot.data[,"FDA"] <- predict.glm(object=results2, newdata=plot.data, type="response")
        }
    }
    if (SVM > 0) {
        results <- SVM.OLD
        if (CATCH.OFF == F) {
            tryCatch(plot.data[,"SVM"] <- kernlab::predict(object=results, newdata=plot.data.vars, type="probabilities")[,2],
                error= function(err) {print(paste("SVM prediction (kernlab package) failed"))},
                silent=F)
        }else{
            plot.data[,"SVM"] <- kernlab::predict(object=results, newdata=plot.data.vars, type="probabilities")[,2]
        }
        results2 <- SVM.PROBIT.OLD
        if (is.null(results2) == F) {
            plot.data[,"SVM.step1"] <- plot.data[,"SVM"]
            plot.data[,"SVM"] <- predict.glm(object=results2, newdata=plot.data, type="response")
        }
    }
    if (SVME > 0) {
        results <- SVME.OLD
        if (CATCH.OFF == F) {
            tryCatch(plot.data[,"SVME"] <- predict.SVME(model=results, newdata=plot.data.vars),
                error= function(err) {print(paste("SVM prediction (e1071 package) failed"))},
                warning= function(war) {print(paste("SVM prediction (e1071 package) failed"))},
                silent=F)
        }else{
            plot.data[,"SVME"] <- predict.SVME(model=results, newdata=plot.data.vars)
        }
        results2 <- SVME.PROBIT.OLD
        if (is.null(results2) == F) {
            plot.data[,"SVME.step1"] <- plot.data[,"SVME"]
            plot.data[,"SVME"] <- predict.glm(object=results2, newdata=plot.data, type="response")
        }
    }
    if (GLMNET > 0) {
        results <- GLMNET.OLD
        if (CATCH.OFF == F) {
            tryCatch(plot.data[,"GLMNET"] <- predict.GLMNET(model=results, newdata=plot.data.numvars, GLMNET.class=GLMNET.class),
                error= function(err) {print(paste("GLMNET prediction (glmnet package) failed"))},
                warning= function(war) {print(paste("GLMNET prediction (glmnet package) failed"))},
                silent=F)
        }else{
            plot.data[,"GLMNET"] <- predict.GLMNET(model=results, newdata=plot.data.numvars, GLMNET.class=GLMNET.class)
        }
        results2 <- GLMNET.PROBIT.OLD
        if (is.null(results2) == F) {
            plot.data[,"GLMNET.step1"] <- plot.data[,"GLMNET"]
            plot.data[,"GLMNET"] <- predict.glm(object=results2, newdata=plot.data, type="response")
        }
    }
    if (BIOCLIM.O > 0) {
        results <- BIOCLIM.O.OLD
        if (CATCH.OFF == F) {
            tryCatch(plot.data[,"BIOCLIM.O"] <- predict.BIOCLIM.O(object=results, newdata=plot.data.vars),
                error= function(err) {print(paste("original BIOCLIM prediction failed"))},
                silent=F)
        }else{
            plot.data[,"BIOCLIM.O"] <- predict.BIOCLIM.O(object=results, newdata=plot.data.vars)
        }
        results2 <- BIOCLIM.O.PROBIT.OLD
        if (is.null(results2) == F) {
            plot.data[,"BIOCLIM.O.step1"] <- plot.data[,"BIOCLIM.O"]
            plot.data[,"BIOCLIM.O"] <- predict.glm(object=results2, newdata=plot.data, type="response")
        }
    }
    if (BIOCLIM > 0 || DOMAIN > 0 || MAHAL > 0) {
        if(is.null(factors)==F) {
            for (i in 1:length(factors)) {
                plot.data.vars <- plot.data.vars[, which(names(plot.data.vars) != factors[i]), drop=F]
            }
        }
    }
    if (BIOCLIM > 0) {
        results <- BIOCLIM.OLD
        if (CATCH.OFF == F) {
            tryCatch(plot.data[,"BIOCLIM"] <- dismo::predict(object=results, x=plot.data.vars),
                error= function(err) {print(paste("BIOCLIM prediction failed"))},
                silent=F)
        }else{
            plot.data[,"BIOCLIM"] <- dismo::predict(object=results, x=plot.data.vars)
        }
        results2 <- BIOCLIM.PROBIT.OLD
        if (is.null(results2) == F) {
            plot.data[,"BIOCLIM.step1"] <- plot.data[,"BIOCLIM"]
            plot.data[,"BIOCLIM"] <- predict.glm(object=results2, newdata=plot.data, type="response")
        }
    }
    if (DOMAIN > 0) {
        if (CATCH.OFF == F) {
            tryCatch(plot.data[,"DOMAIN"] <- dismo::predict(object=results, x=plot.data.domain),
                error= function(err) {print(paste("DOMAIN prediction failed"))},
                silent=F)
        }else{
            plot.data[,"DOMAIN"] <- dismo::predict(object=results, x=plot.data.domain)
        }
        results2 <- DOMAIN.PROBIT.OLD
        if (is.null(results2) == F) {
            plot.data[,"DOMAIN.step1"] <- plot.data[,"DOMAIN"]
            plot.data[,"DOMAIN"] <- predict.glm(object=results2, newdata=plot.data, type="response")
        }
    }
    if (MAHAL > 0) {
        results <- MAHAL.OLD
        results2 <- MAHAL.PROBIT.OLD
        if (is.null(results2) == T) {
            if (CATCH.OFF == F) {
                tryCatch(plot.data[,"MAHAL"] <- predict.MAHAL(model=results, newdata=plot.data.mahal, PROBIT=FALSE),
                    error= function(err) {print(paste("Mahalanobis prediction failed"))},
                    silent=F)
            }else{
                plot.data[,"MAHAL"] <- predict.MAHAL(model=results, newdata=plot.data.mahal, PROBIT=FALSE)
            }
        }else{
            if (CATCH.OFF == F) {
                tryCatch(plot.data[,"MAHAL"] <- predict.MAHAL(model=results, newdata=plot.data.mahal, PROBIT=TRUE),
                    error= function(err) {print(paste("Mahalanobis prediction failed"))},
                    silent=F)
            }else{
                plot.data[,"MAHAL"] <- predict.MAHAL(model=results, newdata=plot.data.mahal, PROBIT=TRUE)
            }
            plot.data[,"MAHAL.step1"] <- plot.data[,"MAHAL"]
            plot.data[,"MAHAL"] <- predict.glm(object=results2, newdata=plot.data, type="response")
        }
    }
    if (MAHAL01 > 0) {
        results <- MAHAL01.OLD
        if (CATCH.OFF == F) {
            tryCatch(plot.data[,"MAHAL01"] <- predict.MAHAL01(model=results, newdata=plot.data.mahal, MAHAL.shape=MAHAL.shape),
                error= function(err) {print(paste("transformed Mahalanobis prediction failed"))},
                silent=F)
        }else{
            plot.data[,"MAHAL01"] <- predict.MAHAL01(model=results, newdata=plot.data.mahal, MAHAL.shape=MAHAL.shape)
        }
        results2 <- MAHAL.PROBIT.OLD
        if (is.null(results2) == F) {
            plot.data[,"MAHAL01.step1"] <- plot.data[,"MAHAL01"]
            plot.data[,"MAHAL01"] <- predict.glm(object=results2, newdata=plot.data, type="response")
        }
    }
#
    plot.data[,"ENSEMBLE"] <- ws["MAXENT"]*plot.data[,"MAXENT"] + ws["GBM"]*plot.data[,"GBM"] +
        ws["GBMSTEP"]*plot.data[,"GBMSTEP"] + ws["RF"]*plot.data[,"RF"] + ws["GLM"]*plot.data[,"GLM"] +
        ws["GLMSTEP"]*plot.data[,"GLMSTEP"] + ws["GAM"]*plot.data[,"GAM"] + ws["GAMSTEP"]*plot.data[,"GAMSTEP"] +
        ws["MGCV"]*plot.data[,"MGCV"] + ws["MGCVFIX"]*plot.data[,"MGCVFIX"] + ws["EARTH"]*plot.data[,"EARTH"] +
        ws["RPART"]*plot.data[,"RPART"] + ws["NNET"]*plot.data[,"NNET"] + ws["FDA"]*plot.data[,"FDA"] +
        ws["SVM"]*plot.data[,"SVM"] + ws["SVME"]*plot.data[,"SVME"] + ws["GLMNET"]*plot.data[,"GLMNET"]
        ws["BIOCLIM.O"]*plot.data[,"BIOCLIM.O"] + ws["BIOCLIM"]*plot.data[,"BIOCLIM"] +
        ws["DOMAIN"]*plot.data[,"DOMAIN"] + ws["MAHAL"]*plot.data[,"MAHAL"] + ws["MAHAL01"]*plot.data[,"MAHAL01"]

#
   for (i in 1:length(modelnames)) {
       if (sum(plot.data[, which(names(plot.data) == modelnames[i])]) == 0) {plot.data <- plot.data[, which(names(plot.data) != modelnames[i]), drop=F]} 
    }

    out <- list(plot.data=plot.data, TrainData=models.list$TrainData)

    remove(plot.data.vars, envir=.BiodiversityR)
    remove(plot.data.numvars, envir=.BiodiversityR)
    remove(plot.data.domain, envir=.BiodiversityR)
    remove(plot.data.mahal, envir=.BiodiversityR)

    return(out)

}

