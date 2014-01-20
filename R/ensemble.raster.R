`ensemble.raster` <- function(
    xn=NULL, ext=NULL,
    models.list=NULL, 
    input.weights=models.list$output.weights,
    thresholds=models.list$thresholds,
    RASTER.species.name="Species001", RASTER.stack.name=xn@title, 
    RASTER.format="raster", RASTER.datatype="INT2S", RASTER.NAflag=-32767,
    RASTER.models.overwrite=TRUE,
    KML.out=FALSE, KML.maxpixels=100000, KML.blur=10,
    evaluate=FALSE, SINK=FALSE,
    p=models.list$p, a=models.list$a,
    pt=models.list$pt, at=models.list$at
)
{
    .BiodiversityR <- new.env()
    if (! require(dismo)) {stop("Please install the dismo package")}
    if (is.null(xn) == T) {stop("value for parameter xn is missing (RasterStack object)")}
    if (is.null(models.list) == T) {stop("provide 'models.list' as models will not be recalibrated and retested")}
    if (is.null(input.weights) == T) {input.weights <- models.list$output.weights}
    if (is.null(thresholds) == T) {stop("provide 'thresholds' as models will not be recalibrated and retested")}
    retest <- F
    if (evaluate == T) { 
        if (is.null(p)==T || is.null(a)==T) {
            cat(paste("\n", "NOTE: not possible to evaluate the models since locations p and a are not provided", "\n", sep = ""))
            evaluate <- F
        }
        if (is.null(pt)==F && is.null(at)==F) {
            if(identical(pt, p) == F || identical(at, a) == F)  {retest <- T}
        }
    }
    if (is.null(ext) == F) {
        if(length(xn@title) == 0) {xn@title <- "stack1"}
        title.old <- xn@title
        xn <- crop(xn, y=ext, snap="in")
        xn@title <- title.old
    }

# create output file
    if (RASTER.species.name == "Species001") {
        RASTER.species.name <- models.list$species.name
    }
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
        cat(paste("\n\n", "RESULTS (ensemble.raster function)", "\n", sep=""), file=paste.file, append=T)
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
        if (any(vars.xn==vars[i]) == F) {stop("explanatory variable '", vars[i], "' not among grid layers of RasterStack xn \n", sep = "")}
    }
    nv <- length(vars.xn) 
    for (i in 1:nv) {
        if (any(vars==vars.xn[i]) == F) {
            cat(paste("\n", "NOTE: RasterStack layer '", vars.xn[i], "' was not calibrated as explanatory variable", "\n", sep = ""))
            xn <- dropLayer(xn, which(names(xn) %in% c(vars.xn[i]) ))
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
        categories <- models.list$categories
    }
    dummy.vars <- models.list$dummy.vars 
#
    KML.blur <- trunc(KML.blur)
    if (KML.blur < 1) {KML.blur <- 1}
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
# prepare for raster output
    dir.create("models", showWarnings = F)
    dir.create("ensembles", showWarnings = F)
    dir.create("ensembles/count", showWarnings = F)
    dir.create("ensembles/presence", showWarnings = F)
    stack.title <- RASTER.stack.name
#    stack.title <- xn@title
    if(KML.out == T) {
        dir.create("kml", showWarnings = F)
        dir.create("kml/count", showWarnings = F)
        dir.create("kml/presence", showWarnings = F)
    }
    rasterfull <- paste("ensembles/", RASTER.species.name, "_", stack.title , sep="")
    kmlfull <- paste("kml/", RASTER.species.name, "_", stack.title , sep="")
    raster.title <- paste(RASTER.species.name, "_", stack.title , sep="")
    rastercount <- paste("ensembles/count/", RASTER.species.name, "_", stack.title , sep="")
    kmlcount <- paste("kml/count/", RASTER.species.name, "_", stack.title , sep="")
    rasterpresence <- paste("ensembles/presence/", RASTER.species.name, "_", stack.title, sep="")
    kmlpresence <- paste("kml/presence/", RASTER.species.name, "_", stack.title, sep="")
    RASTER.species.orig <- RASTER.species.name
    if (RASTER.models.overwrite==T) {
        RASTER.species.name <- "working"
    }else{
        RASTER.species.name <- paste(RASTER.species.name, "_", stack.title, sep="")
    }
#
#
    cat(paste("\n", "Start of modelling for organism: ", RASTER.species.orig, "\n", sep = ""))
    cat(paste("Predictions for RasterStack: ", stack.title, "\n", sep = ""))
    ensemble.statistics <- NULL
    cat(paste("ensemble raster layers will be saved in folder ", getwd(), "/ensembles", "\n\n", sep = ""))
    statistics.names <- c("n.models", "ensemble.threshold", "ensemble.min", "ensemble.max", "count.min", "count.max") 
    ensemble.statistics <- numeric(6)
    names(ensemble.statistics) <- statistics.names
#
# sometimes still error warnings for minimum and maximum values of the layers
# set minimum and maximum values for xn
    for (i in 1:nlayers(xn)) {
        xn[[i]] <- setMinMax(xn[[i]])
    }
#
#
# count models
    mc <- 0
#
# start raster layer creations
    if (ws["MAXENT"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Maximum entropy algorithm (package: dismo)\n", sep=""))
        # Put the file 'maxent.jar' in the 'java' folder of dismo
        # the file 'maxent.jar' can be obtained from from http://www.cs.princeton.edu/~schapire/maxent/.
        jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
        results <- MAXENT.OLD
        pmaxent <- NULL
        fullname <- paste("models/", RASTER.species.name, "_MAXENT", sep="")
        tryCatch(pmaxent <- raster::predict(object=results, x=xn, na.rm=TRUE, 
                filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
            error= function(err) {print(paste("MAXENT prediction failed"))},
            silent=F)
        if (is.null(pmaxent) == F) {
            pmaxent <- trunc(1000*pmaxent)
            writeRaster(x=pmaxent, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- extract(pmaxent, p)
                abs1 <- extract(pmaxent, a)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- extract(pmaxent, pt)
                abs1 <- extract(pmaxent, at)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: MAXENT prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["MAXENT"] <- -1 
        }
    }
    if (ws["GBM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Generalized boosted regression modeling (package: gbm) \n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "WARRNING: not certain whether the correct factor levels will be used", "\n", sep=""))
        } 
        results <- GBM.OLD
        pgbm <- NULL
        fullname <- paste("models/", RASTER.species.name, "_GBM", sep="")
        tryCatch(pgbm <- raster::predict(object=xn, model=results, na.rm=TRUE, factors=categories,
                n.trees=results$n.trees, type="response", 
                filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
            error= function(err) {print(paste("GBM prediction failed"))},
            silent=F)
        if (is.null(pgbm) == F) {
            pgbm <- trunc(1000*pgbm)
            writeRaster(x=pgbm, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- extract(pgbm, p)
                abs1 <- extract(pgbm, a)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- extract(pgbm, pt)
                abs1 <- extract(pgbm, at)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: GBM prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["GBM"] <- -1 
        }
    }
    if (ws["GBMSTEP"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". gbm step algorithm (package: dismo)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "WARRNING: not certain whether the correct factor levels will be used", "\n", sep=""))
        } 
        results <- GBMSTEP.OLD
        pgbms <- NULL
        fullname <- paste("models/", RASTER.species.name, "_GBMSTEP", sep="")
        tryCatch(pgbms <- raster::predict(object=xn, model=results, na.rm=TRUE, factors=categories,
                n.trees=results$n.trees, type="response", 
                filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
            error= function(err) {print(paste("stepwise GBM prediction failed"))},
            silent=F)
        if (is.null(pgbms) == F) {
            pgbms <- trunc(1000*pgbms)
            writeRaster(x=pgbms, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- extract(pgbms, p)
                abs1 <- extract(pgbms, a)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- extract(pgbms, pt)
                abs1 <- extract(pgbms, at)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: stepwise GBM prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["GBMSTEP"] <- -1 
        }
    }
    if (ws["RF"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Random forest algorithm (package: randomForest)\n", sep=""))
        results <- RF.OLD
        prf <- NULL
        fullname <- paste("models/", RASTER.species.name, "_RF", sep="")
        tryCatch(prf <- raster::predict(object=xn, model=results, na.rm=TRUE, factors=categories,
                filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
            error= function(err) {print(paste("random forest prediction failed"))},
            silent=F)
        if (is.null(prf) == F) {
            prf <- trunc(1000*prf)
            writeRaster(x=prf, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- extract(prf, p)
                abs1 <- extract(prf, a)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- extract(prf, pt)
                abs1 <- extract(prf, at)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: random forest prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["RF"] <- -1 
        }
    } 
    if (ws["GLM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Generalized Linear Model \n", sep=""))
        results <- GLM.OLD
        pglm <- NULL
        fullname <- paste("models/", RASTER.species.name, "_GLM", sep="")
        tryCatch(pglm <- raster::predict(object=xn, model=results, na.rm=TRUE, type="response", factors=categories,
                filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
            error= function(err) {print(paste("GLM prediction failed"))},
            silent=F)
        if (is.null(pglm) == F) {
            pglm <- trunc(1000*pglm)
            writeRaster(x=pglm, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- extract(pglm, p)
                abs1 <- extract(pglm, a)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- extract(pglm, pt)
                abs1 <- extract(pglm, at)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: GLM prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["GLM"] <- -1 
        }
    }
    if (ws["GLMSTEP"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Stepwise Generalized Linear Model \n", sep=""))
        results <- GLMSTEP.OLD
        pglms <- NULL
        fullname <- paste("models/", RASTER.species.name, "_GLMSTEP", sep="")
        tryCatch(pglms <- raster::predict(object=xn, model=results, na.rm=TRUE, type="response", factors=categories,
                filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
            error= function(err) {print(paste("stepwise GLM prediction failed"))},
            silent=F)
        if (is.null(pglms) == F) {
            pglms <- trunc(1000*pglms)
            writeRaster(x=pglms, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- extract(pglms, p)
                abs1 <- extract(pglms, a)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- extract(pglms, pt)
                abs1 <- extract(pglms, at)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: stepwise GLM prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["GLMSTEP"] <- -1 
        } 
    }
    if (ws["GAM"] > 0 || ws["GAMSTEP"] > 0) {
        cat(paste("\n"))
        try(detach(package:mgcv), silent=T)
        require(gam, quietly=T)
    }
    if (ws["GAM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Generalized Additive Model (package: gam)\n", sep=""))
        results <- GAM.OLD
        pgam <- NULL
        fullname <- paste("models/", RASTER.species.name, "_GAM", sep="")
        tryCatch(pgam <- raster::predict(object=xn, model=results, na.rm=TRUE, type="response", factors=categories,
                filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
            error= function(err) {print(paste("GAM prediction (gam package) failed"))},
            silent=F)
        if (is.null(pgam) == F) {
            pgam <- trunc(1000*pgam)
            writeRaster(x=pgam, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- extract(pgam, p)
                abs1 <- extract(pgam, a)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- extract(pgam, pt)
                abs1 <- extract(pgam, at)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: GAM prediction (gam package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["GAM"] <- -1 
        } 
    }
    if (ws["GAMSTEP"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Stepwise Generalized Additive Model (package: gam)\n", sep=""))
        results <- GAMSTEP.OLD
        pgams <- NULL
        fullname <- paste("models/", RASTER.species.name, "_GAMSTEP", sep="")
        tryCatch(pgams <- raster::predict(object=xn, model=results, type="response", na.rm=TRUE, factors=categories,
                filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
            error= function(err) {print(paste("stepwise GAM prediction (gam package) failed"))},
            silent=F)
        if (is.null(pgams) == F) {
            pgams <- trunc(1000*pgams)
            writeRaster(x=pgams, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- extract(pgams, p)
                abs1 <- extract(pgams, a)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- extract(pgams, pt)
                abs1 <- extract(pgams, at)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: stepwise GAM prediction (gam package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["GAMSTEP"] <- -1 
        } 
    }
    if (ws["MGCV"] > 0 || ws["MGCVFIX"] > 0) {
        cat(paste("\n"))
        try(detach(package:gam), silent=T)
        require(mgcv, quietly=T)
    }
    if (ws["MGCV"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Generalized Additive Model (package: mgcv)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "WARRNING: not certain whether the correct factor levels will be used", "\n", sep=""))
        } 
        results <- MGCV.OLD
        pmgcv <- NULL
        fullname <- paste("models/", RASTER.species.name, "_MGCV", sep="")
        tryCatch(pmgcv <- raster::predict(object=xn, model=results, na.rm=TRUE, type="response", factors=categories,
                filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
            error= function(err) {print(paste("GAM prediction (mgcv package) failed"))},
            silent=F)
        if (is.null(pmgcv) == F) {
            pmgcv <- trunc(1000*pmgcv)
            writeRaster(x=pmgcv, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- extract(pmgcv, p)
                abs1 <- extract(pmgcv, a)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- extract(pmgcv, pt)
                abs1 <- extract(pmgcv, at)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: GAM prediction (mgcv package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["MGCV"] <- -1 
        } 
    }
    if (ws["MGCVFIX"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". GAM with fixed d.f. regression splines (package: mgcv)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "WARRNING: not certain whether the correct factor levels will be used", "\n", sep=""))
        } 
        results <- MGCVFIX.OLD
        pmgcvf <- NULL
        fullname <- paste("models/", RASTER.species.name, "_MGCVFIX", sep="")
        tryCatch(pmgcvf <- raster::predict(object=xn, model=results, na.rm=TRUE, type="response", factors=categories,
                filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
            error= function(err) {print(paste("MGCVFIX prediction (mgcv package) failed"))},
            silent=F)
        if (is.null(pmgcvf) == F) {
            pmgcvf <- trunc(1000*pmgcvf)
            writeRaster(x=pmgcvf, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- extract(pmgcvf, p)
                abs1 <- extract(pmgcvf, a)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- extract(pmgcvf, pt)
                abs1 <- extract(pmgcvf, at)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: GAM prediction (mgcv package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["MGCVFIX"] <- -1 
        } 
    }
    if (ws["EARTH"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Multivariate Adaptive Regression Splines (package: earth)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "NOTE: MARS (earth package) with factors may require explicit dummy variables", "\n", sep=""))
        } 
        results <- EARTH.OLD
        pearth <- NULL
        fullname <- paste("models/", RASTER.species.name, "_EARTH", sep="")
        tryCatch(pearth <- raster::predict(object=xn, model=results, fun=predict.earth2, na.rm=TRUE, type="response", factors=categories,
                filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
            error= function(err) {print(paste("MARS prediction (earth package) failed"))},
            silent=F)
        if (is.null(pearth) == F) {
            pearth <- trunc(1000*pearth)
            writeRaster(x=pearth, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- extract(pearth, p)
                abs1 <- extract(pearth, a)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- extract(pearth, pt)
                abs1 <- extract(pearth, at)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: MARS prediction (earth package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["EARTH"] <- -1 
        } 
    }
    if (ws["RPART"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Recursive Partitioning And Regression Trees (package: rpart)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "WARRNING: not certain whether the correct factor levels will be used", "\n", sep=""))
        } 
        results <- RPART.OLD
        prpart <- NULL
        fullname <- paste("models/", RASTER.species.name, "_RPART", sep="")
        tryCatch(prpart <- raster::predict(object=xn, model=results, na.rm=TRUE, type="prob", index=2, factors=categories,
                filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
            error= function(err) {print(paste("RPART prediction failed"))},
            silent=F)
        if (is.null(prpart) == F) {
            prpart <- trunc(1000*prpart)
            writeRaster(x=prpart, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- extract(prpart, p)
                abs1 <- extract(prpart, a)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- extract(prpart, pt)
                abs1 <- extract(prpart, at)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: RPART prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["RPART"] <- -1 
        } 
    }
    if (ws["NNET"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Artificial Neural Network (package: nnet)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "WARRNING: not certain whether the correct factor levels will be used", "\n", sep=""))
        } 
        results <- NNET.OLD
        pnnet <- NULL
        fullname <- paste("models/", RASTER.species.name, "_NNET", sep="")
        tryCatch(pnnet <- raster::predict(object=xn, model=results, fun=predict.nnet2, na.rm=TRUE, factors=categories,
                filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
            error= function(err) {print(paste("ANN prediction (nnet package) failed"))},
            silent=F)
        if (is.null(pnnet) == F) {
            pnnet <- trunc(1000*pnnet)
            writeRaster(x=pnnet, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- extract(pnnet, p)
                abs1 <- extract(pnnet, a)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- extract(pnnet, pt)
                abs1 <- extract(pnnet, at)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: ANN prediction (nnet package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["NNET"] <- -1 
        } 
    }
    if (ws["FDA"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Flexible Discriminant Analysis (package: mda)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "WARRNING: not certain whether the correct factor levels will be used", "\n", sep=""))
        } 
        results <- FDA.OLD
        pfda <- NULL
        fullname <- paste("models/", RASTER.species.name, "_FDA", sep="")
        tryCatch(pfda <- raster::predict(object=xn, model=results, na.rm=TRUE, type="posterior", index=2, factors=categories,
                filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
            error= function(err) {print(paste("FDA prediction failed"))},
            silent=F)
        if (is.null(pfda) == F) {
            pfda <- trunc(1000*pfda)
            writeRaster(x=pfda, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- extract(pfda, p)
                abs1 <- extract(pfda, a)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- extract(pfda, pt)
                abs1 <- extract(pfda, at)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: FDA prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["FDA"] <- -1 
        } 
    }

    if (ws["SVM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Support Vector Machines (package: kernlab)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "NOTE: SVM model with factors may require explicit dummy variables", "\n", sep=""))
        } 
        results <- SVM.OLD
        psvm <- NULL
        fullname <- paste("models/", RASTER.species.name, "_SVM", sep="")
        predfun <- as.function(kernlab::predict)
        tryCatch(psvm <- raster::predict(object=xn, model=results, fun=predfun, na.rm=TRUE, type="probabilities", index=2, factors=categories,
                filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
            error= function(err) {print(paste("SVM prediction (kernlab package) failed"))},
            silent=F)
        if (is.null(psvm) == F) {
            psvm <- trunc(1000*psvm)
            writeRaster(x=psvm, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- extract(psvm, p)
                abs1 <- extract(psvm, a)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- extract(psvm, pt)
                abs1 <- extract(psvm, at)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: SVM prediction (kernlab package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["SVM"] <- -1 
        } 
    }

    if (ws["SVME"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Support Vector Machines (package: e1071)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "NOTE: SVME model with factors may require explicit dummy variables", "\n", sep=""))
        }
        results <- SVME.OLD
        psvme <- NULL
        fullname <- paste("models/", RASTER.species.name, "_SVME", sep="")
        tryCatch(psvme <- raster::predict(object=xn, model=results, fun=predict.svme, na.rm=TRUE, factors=categories,
                filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
            error= function(err) {print(paste("SVM prediction (e1071 package) failed"))},
            warning= function(war) {print(paste("SVM prediction (e1071 package) failed"))},
            silent=F)
        if (is.null(psvme) == F) {
            psvme <- trunc(1000*psvme)
            writeRaster(x=psvme, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- extract(psvme, p)
                abs1 <- extract(psvme, a)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- extract(psvme, pt)
                abs1 <- extract(psvme, at)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: SVM prediction (e1071 package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["SVME"] <- -1 
        }
    }
    if (BIOCLIM > 0 || DOMAIN > 0 || MAHAL > 0) {
        if(is.null(factors) == F) {
            xn <- dropLayer(xn, which(names(xn) %in% factors))               
        }
        if(is.null(dummy.vars) == F) {
            xn <- dropLayer(xn, which(names(xn) %in% dummy.vars))               
        }
    }
    if (ws["BIOCLIM"] > 0) {  
        mc <- mc+1
        cat(paste("\n", mc, ". BIOCLIM algorithm (package: dismo)\n", sep=""))
        results <- BIOCLIM.OLD
        pbio <- NULL
        fullname <- paste("models/", RASTER.species.name, "_BIOCLIM", sep="")
        tryCatch(pbio <- dismo::predict(object=results, x=xn, na.rm=TRUE, 
                filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
            error= function(err) {print(paste("BIOCLIM prediction failed"))},
            silent=F)
        if (is.null(pbio) == F) {
            pbio <- trunc(1000*pbio)
            writeRaster(x=pbio, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- extract(pbio, p)
                abs1 <- extract(pbio, a)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- extract(pbio, pt)
                abs1 <- extract(pbio, at)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: BIOCLIM prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["BIOCLIM"] <- -1 
        }
    }
    if (ws["DOMAIN"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". DOMAIN algorithm (package: dismo)\n", sep=""))
        results <- DOMAIN.OLD
        pdom <- NULL
        fullname <- paste("models/", RASTER.species.name, "_DOMAIN", sep="")
        tryCatch(pdom <- dismo::predict(object=results, x=xn, na.rm=TRUE, 
                filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
            error= function(err) {print(paste("DOMAIN prediction failed"))},
            silent=F)
        if (is.null(pdom) == F) {
            pdom <- trunc(1000*pdom)
            writeRaster(x=pdom, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- extract(pdom, p)
                abs1 <- extract(pdom, a)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- extract(pdom, pt)
                abs1 <- extract(pdom, at)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: DOMAIN prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["DOMAIN"] <- -1 
        }
    }
    if (ws["MAHAL"] > 0) {  
        mc <- mc+1
        cat(paste("\n", mc, ". Mahalanobis algorithm (package: dismo)\n", sep=""))
        results <- MAHAL.OLD
        pmahal <- NULL
        fullname <- paste("models/", RASTER.species.name, "_MAHAL", sep="")
# not possible to use the predict.mahal function as raster::predict automatically reverts to dismo::predict for 'DistModel' objects
        tryCatch(pmahal <- dismo::predict(object=results, x=xn, na.rm=TRUE,
                filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
            error= function(err) {print(paste("Mahalanobis prediction failed"))},
            silent=F)
        if (is.null(pmahal) == F) {
            pmahal <- pmahal - 1 - MAHAL.shape
            pmahal <- abs(pmahal)
            pmahal <- MAHAL.shape / pmahal
            pmahal <- trunc(1000*pmahal)
            writeRaster(x=pmahal, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- extract(pmahal, p)
                abs1 <- extract(pmahal, a)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- extract(pmahal, pt)
                abs1 <- extract(pmahal, at)
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: Mahalanobis prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["MAHAL"] <- -1 
        }
    }
    if (prediction.failures == T) {
        cat(paste("\n", "WARNING: some predictions failed","\n", sep = ""))
        cat(paste("\n", "actual weights that were used were (-1 indicates failed predictions):","\n", sep = ""))
        print(ws)
        ws[which(ws==-1)] <- 0
    }
#
# create ensembles
    mc <- mc+1
    cat(paste("\n\n", mc, ". Ensemble algorithm\n", sep=""))
    ensemble.statistics["n.models"] <- sum(as.numeric(ws > 0))
    ensemble <- xn[[1]] == NAvalue(xn[[1]])
    setMinMax(ensemble)
    names(ensemble) <- raster.title
    writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
    enscount <- ensemble
    setMinMax(enscount)
    names(enscount) <- paste(raster.title, "_count", sep="")
    writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    enspresence <- ensemble
    setMinMax(enspresence)
    names(enspresence) <- paste(raster.title, "_presence", sep="")
    writeRaster(x=enspresence, filename=rasterpresence, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    if (ws["MAXENT"] > 0) {
        ensemble <- ensemble + ws["MAXENT"] * pmaxent
        writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)       
        pmaxent <- pmaxent > thresholds["MAXENT"]
        enscount <- enscount + pmaxent
        writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["GBM"] > 0) {
        ensemble <- ensemble + ws["GBM"] * pgbm
        writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)      
        pgbm <- pgbm > thresholds["GBM"]
        enscount <- enscount + pgbm
        writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["GBMSTEP"] > 0) {
        ensemble <- ensemble + ws["GBMSTEP"] * pgbms
        writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)      
        pgbms <- pgbms > thresholds["GBMSTEP"]
        enscount <- enscount + pgbms
        writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["RF"] > 0) {
        ensemble <- ensemble + ws["RF"] * prf
        writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)      
        prf <- prf > thresholds["RF"]
        enscount <- enscount + prf
        writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["GLM"] > 0) {
        ensemble <- ensemble + ws["GLM"] * pglm
        writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)    
        pglm <- pglm > thresholds["GLM"]
        enscount <- enscount + pglm
        writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["GLMSTEP"] > 0) {
        ensemble <- ensemble + ws["GLMSTEP"] * pglms
        writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)      
        pglms <- pglms > thresholds["GLMSTEP"]
        enscount <- enscount + pglms
        writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["GAM"] > 0) {
        ensemble <- ensemble + ws["GAM"] * pgam
        writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pgam <- pgam > thresholds["GAM"]
        enscount <- enscount + pgam
        writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["GAMSTEP"] > 0) {
        ensemble <- ensemble + ws["GAMSTEP"] * pgams
        writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pgams <- pgams > thresholds["GAMSTEP"]
        enscount <- enscount + pgams
        writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["MGCV"] > 0) {
        ensemble <- ensemble + ws["MGCV"] * pmgcv
        writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pmgcv <- pmgcv > thresholds["MGCV"]
        enscount <- enscount + pmgcv
        writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["MGCVFIX"] > 0) {
        ensemble <- ensemble + ws["MGCVFIX"] * pmgcvf
        writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pmgcvf <- pmgcvf > thresholds["MGCVFIX"]
        enscount <- enscount + pmgcvf
        writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["EARTH"] > 0) {
        ensemble <- ensemble + ws["EARTH"] * pearth
        writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pearth <- pearth > thresholds["EARTH"]
        enscount <- enscount + pearth
        writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["RPART"] > 0) {
        ensemble <- ensemble + ws["RPART"] * prpart
        writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        prpart <- prpart > thresholds["RPART"]
        enscount <- enscount + prpart
        writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["NNET"] > 0) {
        ensemble <- ensemble + ws["NNET"] * pnnet
        writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pnnet <- pnnet > thresholds["NNET"]
        enscount <- enscount + pnnet
        writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["FDA"] > 0) {
        ensemble <- ensemble + ws["FDA"] * pfda
        writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pfda <- pfda > thresholds["FDA"]
        enscount <- enscount + pfda
        writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["SVM"] > 0) {
        ensemble <- ensemble + ws["SVM"] * psvm
        writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        psvm <- psvm > thresholds["SVM"]
        enscount <- enscount + psvm
        writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["SVME"] > 0) {
        ensemble <- ensemble + ws["SVME"] * psvme
        writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        psvme <- psvme > thresholds["SVME"]
        enscount <- enscount + psvme
        writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["BIOCLIM"] > 0) {
        ensemble <- ensemble + ws["BIOCLIM"] * pbio
        writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pbio <- pbio > thresholds["BIOCLIM"]
        enscount <- enscount + pbio
        writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["DOMAIN"] > 0) {
        ensemble <- ensemble + ws["DOMAIN"] * pdom
        writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pdom <- pdom > thresholds["DOMAIN"]
        enscount <- enscount + pdom
        writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["MAHAL"] > 0) {
        ensemble <- ensemble + ws["MAHAL"] * pmahal
        writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pmahal <- pmahal > thresholds["MAHAL"]
        enscount <- enscount + pmahal
        writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    ensemble <- trunc(ensemble)
    setMinMax(ensemble)
    ensemble.statistics["ensemble.min"] <- minValue(ensemble)
    ensemble.statistics["ensemble.max"] <- maxValue(ensemble)
    names(ensemble) <- raster.title
    writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
    if (KML.out == T) {
        thresholdx <- thresholds["ENSEMBLE"]
        seq1 <- seq(from = 0, to = thresholdx, length.out = 10)
        seq2 <- seq(from = thresholdx, to = 1000, length.out = 11)
        KML(ensemble, filename=kmlfull, col = c(rainbow(n = 10, start = 0, end = 1/6), rainbow(n = 10, start = 3/6, end = 4/6)), colNA = 0, 
            blur=KML.blur, maxpixels=KML.maxpixels, overwrite=T, breaks = c(seq1, seq2))
    }
    setMinMax(enscount)
    ensemble.statistics["count.min"] <- minValue(enscount)
    ensemble.statistics["count.max"] <- maxValue(enscount)
    names(enscount) <- paste(raster.title, "_count", sep="")
    writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    if (KML.out == T) {
        nmax <- sum(as.numeric(ws > 0))
        if (nmax > 3) {
            KML(enscount, filename=kmlcount, col=c("grey", rainbow(n=(nmax-1), start=0, end=1/3), "blue"),
                colNA=0, blur=10, overwrite=T, breaks=seq(from=-1, to=nmax, by=1))
        }else{
            KML(enscount, filename=kmlcount, col=c("grey", rainbow(n=nmax, start=0, end=1/3)),
                colNA=0, blur=10, overwrite=TRUE, breaks=seq(from=-1, to=nmax, by=1))
        }
    }
    ensemble.statistics["ensemble.threshold"] <- thresholds["ENSEMBLE"]
    enspresence <- ensemble > thresholds["ENSEMBLE"]
    setMinMax(enspresence)
    names(enspresence) <- paste(raster.title, "_presence", sep="")
    writeRaster(x=enspresence, filename=rasterpresence, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    if (KML.out == T) {
        KML(enspresence, filename=kmlpresence, col=c("grey", "green"),
            colNA=0, blur=KML.blur, maxpixels=KML.maxpixels, overwrite=T)
    }
    if(evaluate == T) {
        eval1 <- NULL
        cat(paste("\n", "Evaluation of created ensemble raster layer (", rasterfull, ") at locations p and a", "\n\n", sep = ""))
        pres_consensus <- extract(ensemble, p)
        abs_consensus <- extract(ensemble, a)
        eval1 <- evaluate(p=pres_consensus, a=abs_consensus)
        print(eval1)
    }
    if(retest == T) {
        eval1 <- NULL
        cat(paste("\n", "Evaluation of created ensemble raster layer (", rasterfull, ") at locations pt and at", "\n\n", sep = ""))
        pres_consensus <- extract(ensemble, pt)
        abs_consensus <- extract(ensemble, at)
        eval1 <- evaluate(p=pres_consensus, a=abs_consensus)
        print(eval1)
    }
    cat(paste("\n", "End of modelling for organism: ", RASTER.species.orig, "\n", sep = ""))
    cat(paste("Predictions were made for RasterStack: ", stack.title, "\n\n", sep = ""))
    result <- list(ensemble.statistics=ensemble.statistics, call=match.call() )
    if (SINK==T && OLD.SINK==F) {sink(file=NULL, append=T)}
    return(result)
}

