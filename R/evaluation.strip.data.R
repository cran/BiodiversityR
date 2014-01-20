`evaluation.strip.data` <- function(
    xn=NULL, ext=NULL,
    models.list=NULL, 
    input.weights=models.list$output.weights,
    vars=models.list$vars, factors=models.list$factors, 
    dummy.vars=models.list$dummy.vars, 
    steps=50
)
{
    .BiodiversityR <- new.env()
    if (! require(dismo)) {stop("Please install the dismo package")}
    if (is.null(xn) == T) {stop("value for parameter xn is missing (RasterStack object)")}
    if (is.null(models.list) == T) {stop("provide 'models.list' as models will not be recalibrated and retested")}
    if (is.null(input.weights) == T) {input.weights <- models.list$output.weights}
    ws <- input.weights
    if (is.null(ext) == F) {
        if(length(xn@title) == 0) {xn@title <- "stack1"}
        title.old <- xn@title
        xn <- crop(xn, y=ext, snap="in")
        xn@title <- title.old
    }
#
# check if all variables are present
    if (is.null(vars) == T) {vars <- models.list$vars}
    vars.xn <- names(xn)
    nv <- length(vars) 
    for (i in 1:nv) {
        if (any(vars.xn==vars[i]) == F) {stop("explanatory variable '", vars[i], "' not among grid layers of RasterStack xn \n", sep = "")}
    }
    nv <- length(vars.xn) 
    for (i in 1:nv) {
        if (any(vars==vars.xn[i]) == F) {
            cat(paste("\n", "NOTE: RasterStack layer '", vars.xn[i], "' was not calibrated as explanatory variable", "\n", sep = ""))
            xn <- dropLayer(xn, which(names(xn) %in% c(vars.xn[i]) ))
        }
    }
    if (is.null(factors) == T) {factors <- models.list$factors}
    if (is.null(factors) == F) {
        factors <- as.character(factors)
        for (i in 1:length(factors)) {
            if(any(factors[i] == vars) == F) {
                cat(paste("\n", "Warning: ", factors[i], " not among variables", "\n", "\n", sep=""))
                factors <- factors[which(names(factors) != factors[i])]
            }
            if (length(factors) == 0) {factors <- NULL}
        }
    }
    if (is.null(dummy.vars) == T) {dummy.vars <- models.list$dummy.vars}
    if (is.null(dummy.vars) == F) {
        dummy.vars <- as.character(dummy.vars)
        for (i in 1:length(dummy.vars)) {
            if(any(dummy.vars[i] == vars) == F) {
                cat(paste("\n", "Warning: ", dummy.vars[i], " not among variables", "\n", "\n", sep=""))
                dummy.vars <- dummy.vars[which(names(dummy.vars) != dummy.vars[i])]
            }
            if (length(dummy.vars) == 0) {dummy.vars <- NULL}
        }
    }
#
# set minimum and maximum values for xn
    for (i in 1:nlayers(xn)) {
        xn[[i]] <- setMinMax(xn[[i]])
    }
    if(projection(xn)=="NA") {
        projection(xn) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    }
# declare categorical layers for xn
    factors <- models.list$factors
    if(is.null(factors) == F) {
        for (i in 1:length(factors)) {
            j <- which(names(xn) == factors[i])
            xn[[j]] <- raster::as.factor(xn[[j]])
        }
    }
#
    if (is.null(input.weights)==F) {
        MAXENT <- max(c(input.weights["MAXENT"], -1), na.rm=T)
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
        BIOCLIM <- max(c(input.weights["BIOCLIM"], -1), na.rm=T)
        DOMAIN <- max(c(input.weights["DOMAIN"], -1), na.rm=T)
        MAHAL<- max(c(input.weights["MAHAL"], -1), na.rm=T)
    }

#
    MAXENT.OLD <- GBM.OLD <- GBMSTEP.OLD <- RF.OLD <- GLM.OLD <- GLMSTEP.OLD <- GAM.OLD <- GAMSTEP.OLD <- MGCV.OLD <- NULL
    MGCVFIX.OLD <- EARTH.OLD <- RPART.OLD <- NNET.OLD <- FDA.OLD <- SVM.OLD <- SVME.OLD <- BIOCLIM.OLD <- DOMAIN.OLD <- MAHAL.OLD <- GEODIST.OLD <- NULL
    if (is.null(models.list) == F) {
        MAXENT.OLD <- models.list$MAXENT
        GBM.OLD <- models.list$GBM
        GBMSTEP.OLD <- models.list$GBMSTEP
        RF.OLD <- models.list$RF
        GLM.OLD <- models.list$GLM
        GLMSTEP.OLD <- models.list$GLMSTEP
        GAM.OLD <- models.list$GAM
        GAMSTEP.OLD <- models.list$GAMSTEP
        MGCV.OLD <- models.list$MGCV
        MGCVFIX.OLD <- models.list$MGCVFIX
        EARTH.OLD <- models.list$EARTH
        RPART.OLD <- models.list$RPART
        NNET.OLD <- models.list$NNET
        FDA.OLD <- models.list$FDA
        SVM.OLD <- models.list$SVM
        SVME.OLD <- models.list$SVME
        BIOCLIM.OLD <- models.list$BIOCLIM
        DOMAIN.OLD <- models.list$DOMAIN
        MAHAL.OLD <- models.list$MAHAL
        GEODIST.OLD <- models.list$GEODIST
    }
    if (MAXENT > 0) {
	jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
	if (!file.exists(jar)) {stop('maxent program is missing: ', jar, '\nPlease download it here: http://www.cs.princeton.edu/~schapire/maxent/')}
    }
    if (GBM > 0) {
        if (! require(gbm)) {stop("Please install the gbm package")}
    }
    if (GBMSTEP > 0) {
        if (! require(gbm)) {stop("Please install the gbm package")}
    }
    if (RF > 0) {
        if (! require(randomForest)) {stop("Please install the randomForest package")}
    }
    if (GLMSTEP > 0) {
        if (! require(MASS)) {stop("Please install the MASS package")}
    }
    if (GAM > 0  || GAMSTEP > 0) {
        cat(paste("\n"))
        try(detach(package:mgcv), silent=T)
        if (! require(gam)) {stop("Please install the gam package")}
    }
    if (MGCV > 0 || MGCVFIX > 0) {
        cat(paste("\n"))
        try(detach(package:gam), silent=T)
        cat(paste("\n"))
        if (! require(mgcv)) {stop("Please install the mgcv package")}
#         get the probabilities from MGCV
            predict.mgcv <- function(object, newdata, type="response") {
                p <- predict(object=object, newdata=newdata, type=type)
                return(as.numeric(p))
            } 
    }
    if (EARTH > 0) {
        if (! require(earth)) {stop("Please install the earth package")}
#         get the probabilities from earth
            predict.earth2 <- function(object, newdata, type="response") {
                p <- predict(object=object, newdata=newdata, type=type)
                return(as.numeric(p))
            }
    }
    if (RPART > 0) {
        if (! require(rpart)) {stop("Please install the rpart package")}
    }
    if (NNET > 0) {
        if (! require(nnet)) {stop("Please install the nnet package")}
#         get the probabilities from nnet
            predict.nnet2 <- function(object, newdata, type="raw") {
                p <- predict(object=object, newdata=newdata, type=type)
                return(as.numeric(p))
            }
    }
    if (FDA > 0) {
        if (! require(mda)) {stop("Please install the mda package")}
    }
    if (SVM > 0) {
        if (! require(kernlab)) {stop("Please install the kernlab package")}
    }
    if (SVME > 0) {
        if (! require(e1071)) {stop("Please install the e1071 package")}
#         get the probabilities from svm
            predict.svme <- function(model, newdata) {
                p <- predict(model, newdata, probability=T)
                return(attr(p, "probabilities")[,1])
             }
    }
    if (MAHAL > 0) {
        MAHAL.shape <- models.list$formulae$MAHAL.shape
#         get the probabilities from mahal
            predict.mahal <- function(model, newdata, MAHAL.shape) {
                p <- dismo::predict(object=model, x=newdata)
                p <- p - 1 - MAHAL.shape
                p <- abs(p)
                p <- MAHAL.shape / p
                return(p)
             }
    }
# 
    ws <- input.weights
    prediction.failures <- FALSE
#
# prepare data set
    nvars <- length(vars)
    nnum <- nvars - length(factors)
    nrows <- nnum * steps
    plot.data <- array(dim=c(nnum*steps, nvars+2), NA)
    dimnames(plot.data)[[2]] <- c("focal.var","categorical",vars)
# for categorical variables first
    fixedlevel <- array(dim=c(nvars))
    for (i in 1:nvars) {
        if(any(vars[i] == factors) == T) {
            il <- which(vars.xn == vars[i])       
            tabulation <- data.frame(freq(xn[[il]]))
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
            plot.data[,i+2] <- rep(cellStats(xn[[il]], stat="mean"), nrows)
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
            minv <- minValue(xn[[il]])
            maxv <- maxValue(xn[[il]])
            plot.data[startpos:endpos,i+2] <- seq(from=minv,to=maxv, length.out=steps)
         }
    }
# declare factor variables
    plot.data.vars <- data.frame(plot.data)
    for (i in 1:nvars) {
        if(any(vars[i] == factors) == T) {
            plot.data.vars[,i+2] <- factor(plot.data.vars[,i+2], levels=models.list$categories[[vars[i]]])
        }
    }
    modelnames <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL", "ENSEMBLE")
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
    for (i in 1:nlayers(xn)) {
        xn[[i]] <- setMinMax(xn[[i]])
    }
#
# Different modelling algorithms
#
    if (MAXENT > 0) {
        results <- MAXENT.OLD
        tryCatch(plot.data[,"MAXENT"] <- dismo::predict(object=results, x=plot.data.vars),
            error= function(err) {print(paste("MAXENT prediction failed"))},
            silent=F)
        plot.data[,"MAXENT"] <- trunc(1000*plot.data[,"MAXENT"])
    }
    if (GBM > 0) {
        results <- GBM.OLD
        tryCatch(plot.data[,"GBM"] <- predict(object=results, newdata=plot.data.vars, n.trees=results$n.trees, type="response"),
            error= function(err) {print(paste("GBM prediction failed"))},
            silent=F)
        plot.data[,"GBM"] <- trunc(1000*plot.data[,"GBM"])
    }
    if (GBMSTEP > 0) {
        results <- GBMSTEP.OLD
        tryCatch(plot.data[,"GBMSTEP"] <- predict(object=results, newdata=plot.data.vars, n.trees=results$n.trees, type="response"),
            error= function(err) {print(paste("stepwise GBM prediction failed"))},
            silent=F)
        plot.data[,"GBMSTEP"] <- trunc(1000*plot.data[,"GBMSTEP"])
    }
    if (RF > 0) {
        results <- RF.OLD
        tryCatch(plot.data[,"RF"] <- predict(object=results, newdata=plot.data.vars, type="response"),
            error= function(err) {print(paste("random forest prediction failed"))},
            silent=F)
        plot.data[,"RF"] <- trunc(1000*plot.data[,"RF"])
    }
    if (GLM > 0) {
        results <- GLM.OLD
        tryCatch(plot.data[,"GLM"] <- predict(object=results, newdata=plot.data.vars, type="response"),
            error= function(err) {print(paste("GLM prediction failed"))},
            silent=F)
        plot.data[,"GLM"] <- trunc(1000*plot.data[,"GLM"])
    }
    if (GLMSTEP > 0) {
        results <- GLMSTEP.OLD
        tryCatch(plot.data[,"GLMSTEP"] <- predict(object=results, newdata=plot.data.vars, type="response"),
            error= function(err) {print(paste("stepwise GLM prediction failed"))},
            silent=F)
        plot.data[,"GLMSTEP"] <- trunc(1000*plot.data[,"GLMSTEP"])
    }
    if (GAM > 0 || GAMSTEP > 0) {
        cat(paste("\n\n"))
        try(detach(package:mgcv), silent=T)
        require(gam, quietly=T)
    }
    if (GAM > 0) {
        results <- GAM.OLD
        tryCatch(plot.data[,"GAM"] <- predict(object=results, newdata=plot.data.vars, type="response"),
            error= function(err) {print(paste("GAM prediction (gam package) failed"))},
            silent=F)
        plot.data[,"GAM"] <- trunc(1000*plot.data[,"GAM"])
    }
    if (GAMSTEP > 0) {
        results <- GAMSTEP.OLD
        tryCatch(plot.data[,"GAMSTEP"] <- predict(object=results, newdata=plot.data.vars, type="response"),
            error= function(err) {print(paste("stepwise GAM prediction (gam package) failed"))},
            silent=F)
        plot.data[,"GAMSTEP"] <- trunc(1000*plot.data[,"GAMSTEP"])
    }
    if (MGCV > 0 || MGCVFIX > 0) {
        cat(paste("\n\n"))
        try(detach(package:gam), silent=T)
        require(mgcv, quietly=T)
    }
    if (MGCV > 0) {
        results <- MGCV.OLD
        tryCatch(plot.data[,"MGCV"] <- predict.mgcv(object=results, newdata=plot.data.vars),
            error= function(err) {print(paste("GAM prediction (mgcv package) failed"))},
            silent=F)
        plot.data[,"MGCV"] <- trunc(1000*plot.data[,"MGCV"])
    }
    if (MGCVFIX > 0) {
        results <- MGCVFIX.OLD
        tryCatch(plot.data[,"MGCVFIX"] <- predict.mgcv(object=results, newdata=plot.data.vars),
            error= function(err) {print(paste("MGCVFIX prediction (mgcv package) failed"))},
            silent=F)
        plot.data[,"MGCVFIX"] <- trunc(1000*plot.data[,"MGCVFIX"])
    }
    if (EARTH > 0) {
        results <- EARTH.OLD
        tryCatch(plot.data[,"EARTH"] <- predict.earth2(object=results, newdata=plot.data.vars),
            error= function(err) {print(paste("MARS prediction (earth package) failed"))},
            silent=F)
        plot.data[,"EARTH"] <- trunc(1000*plot.data[,"EARTH"])
    }
    if (RPART > 0) {
        results <- RPART.OLD
        tryCatch(plot.data[,"RPART"] <- predict(object=results, newdata=plot.data.vars, type="prob")[,2],
            error= function(err) {print(paste("RPART prediction failed"))},
            silent=F)
        plot.data[,"RPART"] <- trunc(1000*plot.data[,"RPART"])
    }
    if (NNET > 0) {
        results <- NNET.OLD
        tryCatch(plot.data[,"NNET"] <- predict.nnet2(object=results, newdata=plot.data.vars),
            error= function(err) {print(paste("ANN prediction (nnet package) failed"))},
            silent=F)
        plot.data[,"NNET"] <- trunc(1000*plot.data[,"NNET"])
    }
    if (FDA > 0) {
        results <- FDA.OLD
        tryCatch(plot.data[,"FDA"] <- predict(object=results, newdata=plot.data.vars, type="posterior")[,2],
            error= function(err) {print(paste("FDA prediction failed"))},
            silent=F)
        plot.data[,"FDA"] <- trunc(1000*plot.data[,"FDA"])
    }
    if (SVM > 0) {
        results <- SVM.OLD
        tryCatch(plot.data[,"SVM"] <- kernlab::predict(object=results, newdata=plot.data.vars, type="probabilities")[,2],
            error= function(err) {print(paste("SVM prediction (kernlab package) failed"))},
            silent=F)
        plot.data[,"SVM"] <- trunc(1000*plot.data[,"SVM"])
    }
    if (SVME > 0) {
        results <- SVME.OLD
        tryCatch(plot.data[,"SVME"] <- predict.svme(model=results, newdata=plot.data.vars),
            error= function(err) {print(paste("SVM prediction (e1071 package) failed"))},
            warning= function(war) {print(paste("SVM prediction (e1071 package) failed"))},
            silent=F)
        plot.data[,"SVME"] <- trunc(1000*plot.data[,"SVME"])
    }
    if (BIOCLIM > 0 || DOMAIN > 0 || MAHAL > 0) {
        if(is.null(factors)==F) {
            for (i in 1:length(factors)) {
                plot.data.vars <- plot.data.vars[, which(colnames(plot.data.vars) != factors[i])]
            }
        }
        if(is.null(dummy.vars)==F) {
            for (i in 1:length(dummy.vars)) {
                plot.data.vars <- plot.data.vars[, which(colnames(plot.data.vars) != dummy.vars[i])]
            }
        }
    }
    if (BIOCLIM > 0) {
        results <- BIOCLIM.OLD
        tryCatch(plot.data[,"BIOCLIM"] <- dismo::predict(object=results, x=plot.data.vars),
            error= function(err) {print(paste("BIOCLIM prediction failed"))},
            silent=F)
        plot.data[,"BIOCLIM"] <- trunc(1000*plot.data[,"BIOCLIM"])
    }
    if (DOMAIN > 0) {
        results <- DOMAIN.OLD
        tryCatch(plot.data[,"DOMAIN"] <- dismo::predict(object=results, x=plot.data.vars),
            error= function(err) {print(paste("DOMAIN prediction failed"))},
            silent=F)
        plot.data[,"DOMAIN"] <- trunc(1000*plot.data[,"DOMAIN"])
    }
    if (MAHAL > 0) {
        results <- MAHAL.OLD
        tryCatch(plot.data[,"MAHAL"] <- predict.mahal(model=results, newdata=plot.data.vars, MAHAL.shape=MAHAL.shape),
            error= function(err) {print(paste("Mahalanobis prediction failed"))},
            silent=F)
        plot.data[,"MAHAL"] <- trunc(1000*plot.data[,"MAHAL"])
    }
#
    if (prediction.failures == T) {
        cat(paste("\n", "WARNING: some predictions failed","\n", sep = ""))
        cat(paste("\n", "actual weights that were used were (-1 indicates failed predictions):","\n", sep = ""))
        print(ws)
        ws[which(ws==-1)] <- 0
    }
#
    plot.data[,"ENSEMBLE"] <- ws["MAXENT"]*plot.data[,"MAXENT"] + ws["GBM"]*plot.data[,"GBM"] +
        ws["GBMSTEP"]*plot.data[,"GBMSTEP"] + ws["RF"]*plot.data[,"RF"] + ws["GLM"]*plot.data[,"GLM"] +
        ws["GLMSTEP"]*plot.data[,"GLMSTEP"] + ws["GAM"]*plot.data[,"GAM"] + ws["GAMSTEP"]*plot.data[,"GAMSTEP"] +
        ws["MGCV"]*plot.data[,"MGCV"] + ws["MGCVFIX"]*plot.data[,"MGCVFIX"] + ws["EARTH"]*plot.data[,"EARTH"] +
        ws["RPART"]*plot.data[,"RPART"] + ws["NNET"]*plot.data[,"NNET"] + ws["FDA"]*plot.data[,"FDA"] +
        ws["SVM"]*plot.data[,"SVM"] + ws["SVME"]*plot.data[,"SVME"] + ws["BIOCLIM"]*plot.data[,"BIOCLIM"] +
        ws["DOMAIN"]*plot.data[,"DOMAIN"] + ws["MAHAL"]*plot.data[,"MAHAL"]
    plot.data[,"ENSEMBLE"] <- trunc(plot.data[,"ENSEMBLE"])

    return(plot.data)
}

