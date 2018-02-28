`ensemble.calibrate.models` <- function(
    x=NULL, p=NULL,
    a=NULL, an=1000, excludep=FALSE, 
    k=0, pt=NULL, at=NULL, SSB.reduce=FALSE, CIRCLES.d=250000,
    TrainData=NULL, TestData=NULL,
    VIF=FALSE, COR=FALSE,
    SINK=FALSE, PLOTS=FALSE, CATCH.OFF=FALSE,
    threshold.method="spec_sens", threshold.sensitivity=0.9, threshold.PresenceAbsence=FALSE,
    evaluations.keep=FALSE, 
    models.list=NULL, models.keep=FALSE, 
    models.save=FALSE, species.name="Species001",
    ENSEMBLE.tune=FALSE, 
    ENSEMBLE.best=0, ENSEMBLE.min=0.7, ENSEMBLE.exponent=1.0, ENSEMBLE.weight.min=0.05,
    input.weights=NULL, 
    MAXENT=1, MAXLIKE=1, GBM=1, GBMSTEP=1, RF=1, GLM=1, GLMSTEP=1, 
    GAM=1, GAMSTEP=1, MGCV=1, MGCVFIX=0,
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, GLMNET=1,
    BIOCLIM.O=0, BIOCLIM=1, DOMAIN=1, MAHAL=1, MAHAL01=1,
    PROBIT=FALSE,
    Yweights="BIOMOD", 
    layer.drops=NULL, factors=NULL, dummy.vars=NULL,
    formulae.defaults=TRUE, maxit=100,
    MAXENT.a=NULL, MAXENT.an=10000, 
    MAXENT.path=paste(getwd(), "/models/maxent_", species.name,  sep=""),
    MAXLIKE.formula=NULL, MAXLIKE.method="BFGS",
    GBM.formula=NULL, GBM.n.trees=2001, 
    GBMSTEP.gbm.x=2:(ncol(TrainData.vars)+1), GBMSTEP.tree.complexity=5, GBMSTEP.learning.rate=0.005, 
    GBMSTEP.bag.fraction=0.5, GBMSTEP.step.size=100, 
    RF.formula=NULL, RF.ntree=751, RF.mtry=floor(sqrt(ncol(TrainData.vars))),
    GLM.formula=NULL, GLM.family=binomial(link="logit"),
    GLMSTEP.steps=1000, STEP.formula=NULL, GLMSTEP.scope=NULL, GLMSTEP.k=2,
    GAM.formula=NULL, GAM.family=binomial(link="logit"),
    GAMSTEP.steps=1000, GAMSTEP.scope=NULL, GAMSTEP.pos=1,
    MGCV.formula=NULL, MGCV.select=FALSE,
    MGCVFIX.formula=NULL, 
    EARTH.formula=NULL, EARTH.glm=list(family=binomial(link="logit"), maxit=maxit),
    RPART.formula=NULL, RPART.xval=50,
    NNET.formula=NULL, NNET.size=8, NNET.decay=0.01,
    FDA.formula=NULL,
    SVM.formula=NULL,
    SVME.formula=NULL, 
    GLMNET.nlambda=100, GLMNET.class=FALSE,
    BIOCLIM.O.fraction=0.9,
    MAHAL.shape=1
)
{
    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    k <- as.integer(k)
# check data
    if (is.null(TrainData) == T) {
        if(is.null(x) == T) {stop("value for parameter x is missing (RasterStack object)")}
        if(inherits(x,"RasterStack") == F) {stop("x is not a RasterStack object")}
#        if(raster::projection(x)=="NA") {
#            raster::projection(x) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
#        }
        if(is.null(p) == T) {stop("presence locations are missing (parameter p)")}
    }
    if(is.null(p) == F) {names(p) <- c("x", "y")}
    if(is.null(a) == F) {names(a) <- c("x", "y")}
    if(is.null(pt) == F) {names(pt) <- c("x", "y")}
    if(is.null(at) == F) {names(at) <- c("x", "y")}
    if(is.null(MAXENT.a) == F) {names(MAXENT.a) <- c("x", "y")}
#
    if(models.save==T) {
        models.keep <- TRUE
        dir.create("models", showWarnings = F)
    }

# create output file
    dir.create("outputs", showWarnings = F)
    paste.file <- paste(getwd(), "/outputs/", species.name, "_output.txt", sep="")
    OLD.SINK <- TRUE
    if (sink.number(type="output") == 0) {OLD.SINK <- F}
    if (SINK==T && OLD.SINK==F) {
        if (file.exists(paste.file) == F) {
            cat(paste("\n", "NOTE: results captured in file: ", paste.file, "\n", sep = ""))
        }else{
            cat(paste("\n", "NOTE: results appended in file: ", paste.file, "\n", sep = ""))
        }
        cat(paste(sep="\n\n", "RESULTS (ensemble.calibrate.models function)"), file=paste.file, sep="\n\n", append=T)
        sink(file=paste.file, append=T)
        cat(paste(date()), sep="\n")
        print(match.call())
        cat(paste("   "), sep="\n")
    }

# check TrainData
    if (is.null(TrainData) == F) {
        TrainData <- data.frame(TrainData)
        if (names(TrainData)[1] !="pb") {stop("first column for TrainData should be 'pb' containing presence (1) and absence (0) data")}
        if ((is.null(x) == F) && (raster::nlayers(x) != (ncol(TrainData)-1))) {
            cat(paste("\n", "WARNING: different number of explanatory variables in rasterStack and TrainData", sep = ""))
        }
    }

# modify list of variables
# if TrainData is provided, then this data set takes precedence over raster x in the selection of variables

# 
# modify TrainData if layer.drops
    if (is.null(TrainData) == F) {
        if (is.null(layer.drops) == F) {
            vars <- names(TrainData)
            layer.drops <- as.character(layer.drops)
            dummy.vars <- as.character(dummy.vars)
            nd <- length(layer.drops)
            for (i in 1:nd) {     
                if (any(vars==layer.drops[i]) == FALSE) {
                    cat(paste("\n", "WARNING: variable to exclude '", layer.drops[i], "' not among columns of TrainData", sep = ""))
                }else{
                    cat(paste("\n", "NOTE: variable '", layer.drops[i], "' will not be included as explanatory variable", sep = ""))
                    TrainData <- TrainData[, which(names(TrainData) != layer.drops[i]), drop=F]
                    if (is.null(TestData) == F) {TestData <- TestData[, which(names(TestData) != layer.drops[i]), drop=F]}
                    vars <- names(TrainData)
                    if (length(factors) > 0) {
                        factors <- factors[factors != layer.drops[i]]                        
                    }
                    if (length(dummy.vars) > 0) {
                        dummy.vars <- dummy.vars[dummy.vars != layer.drops[i]]
                    }
                }
            }
            if(length(layer.drops) == 0) {layer.drops <- NULL} 
            if(length(factors) == 0) {factors <- NULL}
            if(length(dummy.vars) == 0) {dummy.vars <- NULL}
        }
        if (is.null(factors) == F) {
            vars <- names(TrainData)
            factors <- as.character(factors)
            nf <- length(factors)
            old.factors <- factors
            for (i in 1:nf) {
                if (any(vars==old.factors[i]) == FALSE) {
                    cat(paste("\n", "WARNING: categorical variable '", old.factors[i], "' not among columns of TrainData", sep = ""))
                    factors <- factors[factors != old.factors[i]]
                }
            }
            if(length(factors) == 0) {factors <- NULL}
        }
        if (is.null(dummy.vars) == F) {
            vars <- names(TrainData)
            dummy.vars <- as.character(dummy.vars)
            nf <- length(dummy.vars)
            old.dummy.vars <- dummy.vars
            for (i in 1:nf) {
                if (any(vars==old.dummy.vars[i]) == FALSE) {
                    cat(paste("\n", "WARNING: dummy variable '", old.dummy.vars[i], "' not among columns of TrainData", sep = ""))
                    dummy.vars <- dummy.vars[dummy.vars != old.dummy.vars[i]]
                }
            }
            if(length(dummy.vars) == 0) {dummy.vars <- NULL}
        }
    }
# 

# modify RasterStack x only if this RasterStack was provided
    if (is.null(x) == F) {
# same variables as TrainData in the rasterstack
        if (is.null(TrainData) == F) {
            vars <- names(TrainData)
            vars <- vars[which(vars!="pb")]
            x <- raster::subset(x, subset=vars)
            x <- raster::stack(x)
        }
        if (is.null(TrainData) == T) {
            if (is.null(layer.drops) == F) {
                vars <- names(x)
                layer.drops <- as.character(layer.drops)
                factors <- as.character(factors)
                dummy.vars <- as.character(dummy.vars)
                nd <- length(layer.drops)
                for (i in 1:nd) {     
                    if (any(vars==layer.drops[i])==FALSE) {
                        cat(paste("\n", "WARNING: variable to exclude '", layer.drops[i], "' not among grid layers", sep = ""))
                    }else{
                        cat(paste("\n", "NOTE: variable '", layer.drops[i], "' will not be included as explanatory variable", sep = ""))
                        x <- raster::dropLayer(x, which(names(x) %in% c(layer.drops[i]) ))
                        x <- raster::stack(x)
                        vars <- names(x)
                        if (length(factors) > 0) {
                            factors <- factors[factors != layer.drops[i]]
                        }
                        if (length(dummy.vars) > 0) {
                            dummy.vars <- dummy.vars[dummy.vars != layer.drops[i]]
                        }
                    }
                }
                if(length(layer.drops) == 0) {layer.drops <- NULL}
                if(length(factors) == 0) {factors <- NULL}
                if(length(dummy.vars) == 0) {dummy.vars <- NULL}
            }
            if (is.null(factors) == F) {
                vars <- names(x)
                factors <- as.character(factors)
                nf <- length(factors)
                old.factors <- factors
                for (i in 1:nf) {
                    if (any(vars==old.factors[i])==FALSE) {
                        cat(paste("\n", "WARNING: categorical variable '", old.factors[i], "' not among grid layers", sep = ""))
                        factors <- factors[factors != old.factors[i]]
                    }
                }
                if(length(factors) == 0) {factors <- NULL}
            }
            if (is.null(dummy.vars) == F) {
                vars <- names(x)
                dummy.vars <- as.character(dummy.vars)
                nf <- length(dummy.vars)
                old.dummy.vars <- dummy.vars
                for (i in 1:nf) {
                    if (any(vars==old.dummy.vars[i]) == FALSE) {
                        cat(paste("\n", "WARNING: dummy variable '", old.dummy.vars[i], "' not among grid layers", sep = ""))
                        dummy.vars <- dummy.vars[dummy.vars != old.dummy.vars[i]]
                    }
                }
                if(length(dummy.vars) == 0) {dummy.vars <- NULL}
            }
        }
        # set minimum and maximum values
            for (i in 1:raster::nlayers(x)) {
                x[[i]] <- raster::setMinMax(x[[i]])
            }
        # declare factor layers
        if(is.null(factors)==F) {
            for (i in 1:length(factors)) {
                j <- which(names(x) == factors[i])
                x[[j]] <- raster::as.factor(x[[j]])
            }
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
    ws <- as.numeric(c(MAXENT, MAXLIKE, GBM, GBMSTEP, RF, GLM, GLMSTEP, GAM, GAMSTEP, MGCV, 
        MGCVFIX, EARTH, RPART, NNET, FDA, SVM, SVME, GLMNET, 
        BIOCLIM.O, BIOCLIM, DOMAIN, MAHAL, MAHAL01))
    names(ws) <- c("MAXENT", "MAXLIKE", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", 
        "MGCVFIX", "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "GLMNET", 
        "BIOCLIM.O", "BIOCLIM", "DOMAIN", "MAHAL", "MAHAL01")
    ws <- ensemble.weights(weights=ws, exponent=1, best=0, min.weight=0)
#
    thresholds <- c(ws, -1)
    names(thresholds) <- c(names(ws), "ENSEMBLE")
    AUC.calibration <- thresholds
    AUC.calibration[] <- NA
    AUC.testing <- AUC.calibration
#
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
        if (is.null(models.list$BIOCLIM.O) == F) {BIOCLIM.O.OLD <- models.list$BIOCLIM.O}
        if (is.null(models.list$BIOCLIM) == F) {BIOCLIM.OLD <- models.list$BIOCLIM}
        if (is.null(models.list$DOMAIN) == F) {DOMAIN.OLD <- models.list$DOMAIN}
        if (is.null(models.list$MAHAL) == F) {MAHAL.OLD <- models.list$MAHAL}
        if (is.null(models.list$MAHAL01) == F) {MAHAL01.OLD <- models.list$MAHAL01}
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

# check formulae and packages
    if (ws["MAXENT"] > 0) {
        jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
        if (!file.exists(jar)) {stop('maxent program is missing: ', jar, '\n', 'Please download it here: http://www.cs.princeton.edu/~schapire/maxent/')}
    }
    if (formulae.defaults == T) {
        if (is.null(TrainData) == T) {
            formulae <- ensemble.formulae(x, layer.drops=layer.drops, factors=factors, dummy.vars=dummy.vars, weights=ws)
        }else{
            formulae <- ensemble.formulae(TrainData, layer.drops=layer.drops, factors=factors, dummy.vars=dummy.vars, weights=ws)
        }
    }
    if (ws["MAXLIKE"] > 0) {
        if (! requireNamespace("maxlike")) {stop("Please install the maxlike package")}
        if (is.null(MAXLIKE.formula) == T && formulae.defaults == T) {MAXLIKE.formula <- formulae$MAXLIKE.formula}
        if (is.null(MAXLIKE.formula) == T) {
            cat(paste("\n", "NOTE: not possible to calibrate MAXLIKE as no explanatory variables available" ,"\n", sep = ""))
            ws["MAXLIKE"] <- 0
        }else{
            environment(MAXLIKE.formula) <- .BiodiversityR
        }
    }
    if (ws["GBM"] > 0) {
        if (! requireNamespace("gbm")) {stop("Please install the gbm package")}
	requireNamespace("splines")
        if (is.null(GBM.formula) == T && formulae.defaults == T) {GBM.formula <- formulae$GBM.formula}
        if (is.null(GBM.formula) == T) {
            cat(paste("\n", "NOTE: not possible to calibrate GBM as no explanatory variables available" ,"\n", sep = ""))
            ws["GBM"] <- 0
        }else{
            environment(GBM.formula) <- .BiodiversityR
        }
    }
    if (ws["GBMSTEP"] > 0) {
#        if (! require(gbm)) {stop("Please install the gbm package")}
    }
    if (ws["RF"] > 0) {
#        if (! require(randomForest)) {stop("Please install the randomForest package")}
        if (is.null(RF.formula) == T && formulae.defaults == T) {RF.formula <- formulae$RF.formula}
        if (is.null(RF.formula) == T) {
            cat(paste("\n", "NOTE: not possible to calibrate RF as no explanatory variables available" ,"\n", sep = ""))
            ws["RF"] <- 0
        }else{
            environment(RF.formula) <- .BiodiversityR
            if (identical(RF.ntree, trunc(RF.ntree/2)) == F) {RF.ntree <- RF.ntree + 1}
#  get the probabilities from RF
            predict.RF <- function(object, newdata) {
                p <- predict(object=object, newdata=newdata, type="response")
                return(as.numeric(p))
            }
        }
    }
    if (ws["GLM"] > 0) {
        if (is.null(GLM.formula) == T && formulae.defaults == T) {GLM.formula <- formulae$GLM.formula}
        if (is.null(GLM.formula) == T) {
            cat(paste("\n", "NOTE: not possible to calibrate GLM as no explanatory variables available" ,"\n", sep = ""))
            ws["GLM"] <- 0
        }else{
            environment(GLM.formula) <- .BiodiversityR
            assign("GLM.family", GLM.family, envir=.BiodiversityR)
        }
    }
    if (ws["GLMSTEP"] > 0) {
#        if (! require(MASS)) {stop("Please install the MASS package")}
        if (is.null(STEP.formula) == T && formulae.defaults == T) {STEP.formula <- formulae$STEP.formula}
        if (is.null(GLMSTEP.scope) == T && formulae.defaults == T) {GLMSTEP.scope <- formulae$GLMSTEP.scope}
        if (is.null(STEP.formula) == T) {
            cat(paste("\n", "NOTE: not possible to calibrate GLMSTEP as no explanatory variables available" ,"\n", sep = ""))
            ws["GLMSTEP"] <- 0
        }else{
            environment(STEP.formula) <- .BiodiversityR
            assign("GLM.family", GLM.family, envir=.BiodiversityR)
        }
    }
    if (ws["GAM"] > 0) {
        if (is.null(GAM.formula) == T && formulae.defaults == T) {GAM.formula <- formulae$GAM.formula}
        if (is.null(GAM.formula) == T) {
            cat(paste("\n", "NOTE: not possible to calibrate GAM as no explanatory variables available" ,"\n", sep = ""))
            ws["GAM"] <- 0
        }else{
            environment(GAM.formula) <- .BiodiversityR
            assign("GAM.family", GAM.family, envir=.BiodiversityR)
        }
    }
    if (ws["GAMSTEP"] > 0) {
        if (is.null(STEP.formula) == T && formulae.defaults == T) {STEP.formula <- formulae$STEP.formula}
        if (is.null(GAMSTEP.scope) == T && formulae.defaults == T) {GAMSTEP.scope <- formulae$GAMSTEP.scope}
        if (is.null(STEP.formula) == T) {
            cat(paste("\n", "NOTE: not possible to calibrate GAMSTEP as no explanatory variables available" ,"\n", sep = ""))
            ws["GAMSTEP"] <- 0
        }else{
            environment(STEP.formula) <- .BiodiversityR
            assign("GAM.family", GAM.family, envir=.BiodiversityR)
        }
    }
    if (ws["MGCV"] > 0 || ws["MGCVFIX"] > 0) {
        cat(paste("\n\n"))
#         get the probabilities from MGCV
            predict.MGCV <- function(object, newdata, type="response") {
                p <- mgcv::predict.gam(object=object, newdata=newdata, type=type)
                return(as.numeric(p))
            }
#        options(warn=0)
    }
    if (ws["MGCV"] > 0) {
        if (is.null(MGCV.formula) == T && formulae.defaults == T) {MGCV.formula <- formulae$MGCV.formula}
        if (is.null(MGCV.formula) == T) {
            cat(paste("\n", "NOTE: not possible to calibrate MGCV as no explanatory variables available" ,"\n", sep = ""))
            ws["MGCV"] <- 0
        }else{
            environment(MGCV.formula) <- .BiodiversityR
            assign("GAM.family", GAM.family, envir=.BiodiversityR)
        }
    }
    if (ws["MGCVFIX"] > 0) {
        if (is.null(MGCVFIX.formula) == T && formulae.defaults == T) {MGCVFIX.formula <- formulae$MGCVFIX.formula}
        if (is.null(MGCVFIX.formula) == T) {stop("Please provide the MGCVFIX.formula (hint: use ensemble.formulae function)")}

        if (is.null(MGCVFIX.formula) == T) {
            cat(paste("\n", "NOTE: not possible to calibrate MGCVFIX as no explanatory variables available" ,"\n", sep = ""))
            ws["MGCVFIX"] <- 0
        }else{
            environment(MGCVFIX.formula) <- .BiodiversityR
            assign("GAM.family", GAM.family, envir=.BiodiversityR)
        }
    }
    if (ws["EARTH"] > 0) {
#        if (! require(earth)) {stop("Please install the earth package")}
        if (is.null(EARTH.formula) == T && formulae.defaults == T) {EARTH.formula <- formulae$EARTH.formula}
        if (is.null(EARTH.formula) == T) {
            cat(paste("\n", "NOTE: not possible to calibrate EARTH as no explanatory variables available" ,"\n", sep = ""))
            ws["EARTH"] <- 0
        }else{
            environment(EARTH.formula) <- .BiodiversityR
#   get the probabilities from earth
            predict.EARTH <- function(object, newdata, type="response") {
                p <- predict(object=object, newdata=newdata, type=type)
                return(as.numeric(p))
            }
        }
    }
    if (ws["RPART"] > 0) {
#        if (! require(rpart)) {stop("Please install the rpart package")}
        if (is.null(RPART.formula) == T && formulae.defaults == T) {RPART.formula <- formulae$RPART.formula}
        if (is.null(RPART.formula) == T) {
            cat(paste("\n", "NOTE: not possible to calibrate RPART as no explanatory variables available" ,"\n", sep = ""))
            ws["RPART"] <- 0
        }else{
            environment(RPART.formula) <- .BiodiversityR
        }
    }
    if (ws["NNET"] > 0) {
#        if (! require(nnet)) {stop("Please install the nnet package")}
        if (is.null(NNET.formula) == T && formulae.defaults == T) {NNET.formula <- formulae$NNET.formula}
        if (is.null(NNET.formula) == T) {
            cat(paste("\n", "NOTE: not possible to calibrate NNET as no explanatory variables available" ,"\n", sep = ""))
            ws["NNET"] <- 0
        }else{
             environment(NNET.formula) <- .BiodiversityR
#   get the probabilities from nnet
            predict.NNET <- function(object, newdata, type="raw") {
                p <- predict(object=object, newdata=newdata, type=type)
                return(as.numeric(p))
            }
        }
    }
    if (ws["FDA"] > 0) {
#        if (! require(mda)) {stop("Please install the mda package")}
        if (is.null(FDA.formula) == T && formulae.defaults == T) {FDA.formula <- formulae$FDA.formula}
        if (is.null(FDA.formula) == T) {
            cat(paste("\n", "NOTE: not possible to calibrate FDA as no explanatory variables available" ,"\n", sep = ""))
            ws["FDA"] <- 0
        }else{
            environment(FDA.formula) <- .BiodiversityR
        }
    }
    if (ws["SVM"] > 0) {
#        if (! require(kernlab)) {stop("Please install the kernlab package")}
        if (is.null(SVM.formula) == T && formulae.defaults == T) {SVM.formula <- formulae$SVM.formula}
        if (is.null(SVM.formula) == T) {
            cat(paste("\n", "NOTE: not possible to calibrate SVM as no explanatory variables available" ,"\n", sep = ""))
            ws["SVM"] <- 0
        }else{
            environment(SVM.formula) <- .BiodiversityR
        }
    }
    if (ws["SVME"] > 0) {
#        if (! require(e1071)) {stop("Please install the e1071 package")}
        if (is.null(SVME.formula) == T && formulae.defaults == T) {SVME.formula <- formulae$SVME.formula}
        if (is.null(SVME.formula) == T) {stop("Please provide the SVME.formula (hint: use ensemble.formulae function)")}

        if (is.null(SVME.formula) == T) {
            cat(paste("\n", "NOTE: not possible to calibrate SVME as no explanatory variables available" ,"\n", sep = ""))
            ws["SVME"] <- 0
        }else{
            environment(SVME.formula) <- .BiodiversityR
#  get the probabilities from svm
            predict.SVME <- function(model, newdata) {
                p <- predict(model, newdata, probability=T)
                return(attr(p, "probabilities")[,1])
             }
        }
    }
    if (ws["GLMNET"] > 0) {
        if (! requireNamespace("glmnet")) {stop("Please install the glmnet package")}
#         get the mean probabilities from glmnet
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
    if (ws["BIOCLIM.O"] > 0) {
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
                datai <- newdata[i,]
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
    if (ws["MAHAL"] > 0) {
#         get the probabilities from mahal
            predict.MAHAL <- function(model, newdata, PROBIT) {
                p <- dismo::predict(object=model, x=newdata)
                if (PROBIT == F) {
                    p[p<0] <- 0
                    p[p>1] <- 1
                }
                return(as.numeric(p))
             }
    }
    if (ws["MAHAL01"] > 0) {
#         get the probabilities from transformed mahal
            predict.MAHAL01 <- function(model, newdata, MAHAL.shape) {
                p <- dismo::predict(object=model, x=newdata)
                p <- p - 1 - MAHAL.shape
                p <- abs(p)
                p <- MAHAL.shape / p
                return(p)
             }
    }

# create TrainData and TestData

#    if(length(ENSEMBLE.exponent) > 1 || length(ENSEMBLE.best) > 1 || length(ENSEMBLE.min) > 1) {ENSEMBLE.tune <- TRUE}
    no.tests <- FALSE
    if (is.null(pt)==T && is.null(at)==T && is.null(TestData)==T && k < 2) {
        no.tests <- TRUE
        if (ENSEMBLE.tune == T) {
            cat(paste("\n", "WARNING: not possible to tune (ENSEMBLE.tune=FALSE) as no Test Data available or will be created", sep = ""))
            ENSEMBLE.tune <- FALSE
        }
    }

    if (is.null(TrainData) == F) {
        if(any(is.na(TrainData))) {
            cat(paste("\n", "WARNING: sample units with missing data removed from calibration data","\n\n",sep = ""))
        }        
        TrainValid <- complete.cases(TrainData)
        TrainData <- TrainData[TrainValid,,drop=F]
        if (is.null(p) == F) {
            TrainValid <- complete.cases(TrainData[TrainData[,"pb"]==1,])
            p <- p[TrainValid,]
        }
        if (is.null(a) == F) {
            TrainValid <- complete.cases(TrainData[TrainData[,"pb"]==0,])            
            a <- a[TrainValid,]
        }
        if(is.null(TestData) == T) {
            TestData <- TrainData
            if (k > 1) {
                groupp <- dismo::kfold(TrainData, k=k, by=TrainData[,"pb"])
                TrainData.c <- TrainData[groupp != 1,]
                TestData <- TrainData[groupp == 1,]
                TrainData <- TrainData.c
            }
        } 
    }else{
        if (is.null(a)==T) {
            if (excludep == T) {
                a <- dismo::randomPoints(x[[1]], n=an, p=p, excludep=T)
            }else{
                a <- dismo::randomPoints(x[[1]], n=an, p=NULL, excludep=F)
            }        
        }
        if (is.null(pt)==T && is.null(TestData)) {pt <- p}
        if (k > 1 && identical(pt, p) == T) {
            groupp <- dismo::kfold(p, k=k)
            pc <- p[groupp != 1,]
            pt <- p[groupp == 1,]
            p <- pc
        }
        if (is.null(at)==T && is.null(TestData)) {at <- a}
        if (k > 1 && identical(at, a) == T) {
            groupa <- dismo::kfold(a, k=k)
            ac <- a[groupa != 1,]
            at <- a[groupa == 1,]
            a <- ac
        }
# check for spatial sorting bias (are the testing absences farther away than testing presences)
        if (is.null(p)==F && identical(pt, p)==F) {
            sb.bias <- dismo::ssb(p=pt, a=at, reference=p)
            sb.bias2 <- sb.bias[, 1]/sb.bias[, 2]
            cat(paste("\n", "Spatial sorting bias (dismo package, no bias=1, extreme bias=0): ",  sb.bias2, "\n", sep = ""))
        }
# attempt to reduce spatial bias by searching absence testing locations within circles around all known presences
        if (SSB.reduce == T) {
            if (identical(pt, p) == T) {
                cat(paste("\n", "No search for testing absences in circular neighbourhoods since no separate testing presences", sep = ""))
                SSB.reduce <- FALSE
            }else{
                d.km <- CIRCLES.d/1000
                cat(paste("\n", "Random selection of testing absences in circular neighbourhoods of ", d.km, " km", sep = ""))
                pres_all <- rbind(pt, p)
                circles.calibrate <- dismo::circles(p=pres_all, lonlat=raster::isLonLat(x[[1]]), d=CIRCLES.d)
                circles.predicted <- dismo::predict(circles.calibrate, x[[1]])
                at <- dismo::randomPoints(circles.predicted, n=nrow(at), p=pres_all, excludep=T)
                sb.bias <- dismo::ssb(p=pt, a=at, reference=p)
                sb.bias2 <- sb.bias[, 1]/sb.bias[, 2]
                cat(paste("\n", "Spatial sorting bias with new testing absences: ",  sb.bias2, "\n", sep = ""))
            }
        }
#
        if (length(names(x)) == 1) {
            xdouble <- raster::stack(x, x)
            TrainData <- dismo::prepareData(x=xdouble, p, b=a, factors=factors, xy=FALSE)
            TrainData <- TrainData[, -3]
            names(TrainData)[2] <- names(x)
        }else{
            TrainData <- dismo::prepareData(x, p, b=a, factors=factors, xy=FALSE)
        }
        if(any(is.na(TrainData[TrainData[,"pb"]==1,]))) {
            cat(paste("\n", "WARNING: presence locations with missing data removed from calibration data","\n\n",sep = ""))
        }
        TrainValid <- complete.cases(TrainData[TrainData[,"pb"]==1,])
        p <- p[TrainValid,]
        if(any(is.na(TrainData[TrainData[,"pb"]==0,]))) {
            cat(paste("\n", "WARNING: background locations with missing data removed from calibration data","\n\n",sep = ""))
        }
        TrainValid <- complete.cases(TrainData[TrainData[,"pb"]==0,])
        a <- a[TrainValid,]
        if (length(names(x)) == 1) {
            xdouble <- raster::stack(x, x)
            TrainData <- dismo::prepareData(x=xdouble, p, b=a, factors=factors, xy=FALSE)
            TrainData <- TrainData[, -3]
            names(TrainData)[2] <- names(x)
        }else{
            TrainData <- dismo::prepareData(x, p, b=a, factors=factors, xy=FALSE)
        }
    }
#
    if (no.tests == F) {
        if (is.null(TestData) == F) {
            TestData <- data.frame(TestData)
            if(any(is.na(TestData))) {
                cat(paste("\n", "WARNING: sample units with missing data removed from testing data","\n\n",sep = ""))
            }
            TestValid <- complete.cases(TestData)
            TestData <- TestData[TestValid,,drop=F]
            if (is.null(pt) == F) {
                TestValid <- complete.cases(TestData[TestData[,"pb"]==1,])
                pt <- pt[TestValid,]
            }
            if (is.null(at) == F) {
                TestValid <- complete.cases(TestData[TestData[,"pb"]==0,])
                at <- at[TestValid,]
            }
            if (all(names(TestData)!="pb") == T) {stop("one column needed of 'pb' with presence and absence for TestData")} 
        }else{
            if (length(names(x)) == 1) {
                xdouble <- raster::stack(x, x)
                TestData <- dismo::prepareData(x=xdouble, p=pt, b=at, factors=factors, xy=FALSE)
                TestData <- TestData[, -3]
                names(TestData)[2] <- names(x)
            }else{
                TestData <- dismo::prepareData(x, p=pt, b=at, factors=factors, xy=FALSE)
            }
            if(any(is.na(TestData[TestData[,"pb"]==1,]))) {
                cat(paste("\n", "WARNING: presence locations with missing data removed from evaluation data","\n\n",sep = ""))
            }
            TestValid <- complete.cases(TestData[TestData[,"pb"]==1,])
            pt <- pt[TestValid,]
            if(any(is.na(TestData[TestData[,"pb"]==0,]))) {
                cat(paste("\n", "WARNING: background locations with missing data removed from evaluation data","\n\n",sep = ""))
            }
            TestValid <- complete.cases(TestData[TestData[,"pb"]==0,])
            at <- at[TestValid,]
            TestData <- dismo::prepareData(x, p=pt, b=at, factors=factors, xy=FALSE)
            if (length(names(x)) == 1) {
                xdouble <- raster::stack(x, x)
                TestData <- dismo::prepareData(x=xdouble, p=pt, b=at, factors=factors, xy=FALSE)
                TestData <- TestData[, -3]
                names(TestData)[2] <- names(x)
            }else{
                TestData <- dismo::prepareData(x, p=pt, b=at, factors=factors, xy=FALSE)
            }
        }
    }
#
# check if TestData is different from TrainData
    if (identical(TrainData, TestData) == T) {no.tests <- TRUE}
#
# include all possible factor levels in TrainData (especially important if models are kept)
    if (is.null(factors)==F && is.null(x)==T && no.tests==F) {
        for (i in 1:length(factors)) {
            if (identical(levels(droplevels(TrainData[,factors[i]])), levels(droplevels(TestData[,factors[i]])))==F) {
                cat(paste("\n", "WARNING: differences in factor levels between calibration and evaluation data (variable ", factors[i], ")", "\n", sep = ""))
                cat(paste("Same levels set for both data sets to avoid problems with some evaluations", "\n", sep = ""))
                cat(paste("However, some predictions may still fail", "\n", sep = ""))
                uniquelevels <- unique(c(levels(droplevels(TrainData[,factors[i]])), levels(droplevels(TestData[,factors[i]]))))
                levels(TrainData[,factors[i]]) <- uniquelevels
                levels(TestData[,factors[i]]) <- uniquelevels
            }
        }
    }
    if(is.null(factors)==F && is.null(x)==F) {
        if(models.keep==T) {
            categories <- as.list(factors)
            names(categories) <- factors
        }
        for (i in 1:length(factors)) {
            all.categories <- raster::freq(x[[which(names(x) == factors[i])]])[,1]
            all.categories <- all.categories[is.na(all.categories) == F]
            all.categories <- as.character(all.categories)
            if(models.keep==T) {
                categories[[as.name(factors[i])]] <- all.categories
            }
            train.categories <- levels(droplevels(TrainData[,factors[i]]))
            new.categories <- c(all.categories[is.na(match(all.categories, train.categories))])
            if (length(new.categories) > 0) {
                cat(paste("\n", "The following levels were initially not captured by TrainData for factor '", factors[i], "'\n", sep = ""))
                print(new.categories)
                if (is.null(x)==F && is.null(p)==F && is.null(a)==F) {
# step 1: search if suitable presence locations in TestData
                    if (no.tests == F){
                        for (j in 1:length(new.categories)) {
                            if (any(TestData[TestData[,"pb"]==1, factors[i]] == new.categories[j])) {    
                                cat(paste("Warning: level '", new.categories[j], "' available for presence location in Test Data", "\n", sep = ""))
                            }
                        }
                    }
# step 2: stratified background sample
                    strat1 <- raster::sampleStratified(x[[which(names(x) == factors[i])]], size=1, exp=1, na.rm=TRUE, xy=FALSE)
                    strat1 <- strat1[which(strat1[,2] %in% new.categories), 1]
                    xy1 <- raster::xyFromCell(x[[which(names(x) == factors[i])]], cell=strat1, spatial=FALSE)
                    a <- rbind(a, xy1)
                    if (length(names(x)) == 1) {
                        xdouble <- raster::stack(x, x)
                        TrainData <- dismo::prepareData(x=xdouble, p, b=a, factors=factors, xy=FALSE)
                        TrainData <- TrainData[, -3]
                        names(TrainData)[2] <- names(x)
                    }else{
                        TrainData <- dismo::prepareData(x, p, b=a, factors=factors, xy=FALSE)
                    }
                    TrainValid <- complete.cases(TrainData[TrainData[,"pb"]==0,])
                    a <- a[TrainValid,]
                    if (length(names(x)) == 1) {
                        xdouble <- raster::stack(x, x)
                        TrainData <- dismo::prepareData(x=xdouble, p, b=a, factors=factors, xy=FALSE)
                        TrainData <- TrainData[, -3]
                        names(TrainData)[2] <- names(x)
                    }else{
                        TrainData <- dismo::prepareData(x, p, b=a, factors=factors, xy=FALSE)
                    }
                    train.categories <- levels(droplevels(TrainData[,factors[i]]))
                    new.categories <- all.categories[is.na(match(all.categories, train.categories))]
                    if (length(new.categories) == 0) {
                        cat(paste("All levels have now been included as background data for TrainData for factor '", factors[i], "'\n", sep = ""))
                    }else{
                       cat(paste("The following levels were not captured by TrainData for factor '", factors[i], "'\n", sep = ""))
                       print(new.categories)
                       cat(paste("\n", "Attempt to include these levels was complicated by missing values in other layers", "\n", sep = ""))
                    }
                }
            }
# step 3: also modify test data, but only if no circular neighbourhood
            if (no.tests == F) {
                test.categories <- levels(droplevels(TestData[,factors[i]]))
                new.categories <- c(all.categories[is.na(match(all.categories, test.categories))])
                if (length(new.categories)>0  && SSB.reduce==F) {
                    cat(paste("\n", "The following levels were initially not captured by TestData for factor '", factors[i], "'\n", sep = ""))
                    print(new.categories)
                    if (is.null(x)==F && is.null(pt)==F && is.null(at)==F) {
                        strat1 <- raster::sampleStratified(x[[which(names(x) == factors[i])]], size=1, exp=1, na.rm=TRUE, xy=FALSE)
                        strat1 <- strat1[which(strat1[,2] %in% new.categories), 1]
                        xy1 <- raster::xyFromCell(x[[which(names(x) == factors[i])]], cell=strat1, spatial=FALSE)
                        at <- rbind(at, xy1)
                        if (length(names(x)) == 1) {
                            xdouble <- raster::stack(x, x)
                            TestData <- dismo::prepareData(x=xdouble, p=pt, b=at, factors=factors, xy=FALSE)
                            TestData <- TestData[, -3]
                            names(TestData)[2] <- names(x)
                        }else{
                            TestData <- dismo::prepareData(x, p=pt, b=at, factors=factors, xy=FALSE)
                        }
                        TestValid <- complete.cases(TestData[TestData[,"pb"]==0,])
                        at <- at[TestValid,]
                        if (length(names(x)) == 1) {
                            xdouble <- raster::stack(x, x)
                            TestData <- dismo::prepareData(x=xdouble, p=pt, b=at, factors=factors, xy=FALSE)
                            TestData <- TestData[, -3]
                            names(TestData)[2] <- names(x)
                        }else{
                            TestData <- dismo::prepareData(x, p=pt, b=at, factors=factors, xy=FALSE)
                        }
                        test.categories <- levels(droplevels(TestData[,factors[i]]))
                        new.categories <- all.categories[is.na(match(all.categories, test.categories))]
                        if (length(new.categories) == 0) {
                            cat(paste("All levels have now been included as background data for TestData for factor '", factors[i], "'\n", sep = ""))
                        }else{
                            cat(paste("The following levels were not captured by TestData for factor '", factors[i], "'\n", sep = ""))
                            print(new.categories)
                            cat(paste("\n", "Attempt to include these levels was complicated by missing values in other layers", "\n", sep = ""))
                        }
                    }
                }
                if (length(new.categories)>0  && SSB.reduce==T) {
                    cat(paste("\n", "Note that the following levels were not captured in the circular neighbourhood by TestData for factor '", factors[i], "'\n", sep = ""))
                    print(new.categories)
                }
            }
        }
    }
#
    if (sum(ws, na.rm=T) > 0) {
        cat(paste("\n", "Summary of Training data set used for calibrations (rows: ", nrow(TrainData),  ")\n", sep = ""))
        print(summary(TrainData))
        if (no.tests == F) {
            cat(paste("\n", "Summary of Testing data set used for evaluations (rows: ", nrow(TestData),  ")\n", sep = ""))
            print(summary(TestData))
        }else{
            cat(paste("\n", "(no tests with separate data set)", "\n", sep = ""))
        }
    }
#
    dummy.vars.noDOMAIN <- NULL
#
    if(models.keep == T) {
        models <- list(MAXENT=NULL, MAXLIKE=NULL, GBM=NULL, GBMSTEP=NULL, RF=NULL, GLM=NULL, 
            GLMSTEP=NULL, GAM=NULL, GAMSTEP=NULL, MGCV=NULL, MGCVFIX=NULL, EARTH=NULL, RPART=NULL, 
            NNET=NULL, FDA=NULL, SVM=NULL, SVME=NULL, GLMNET=NULL, BIOCLIM.O=NULL, BIOCLIM=NULL, DOMAIN=NULL, MAHAL=NULL, MAHAL01=NULL,
            formulae=NULL, output.weights=NULL,
            TrainData=NULL, TestData=NULL, p=NULL, a=NULL, pt=NULL, at=NULL, 
            MAXENT.a=NULL, 
            vars=NULL, factors=NULL, categories=NULL, dummy.vars=NULL, dummy.vars.noDOMAIN=NULL, 
            thresholds=NULL, threshold.method=NULL, threshold.sensitivity=NULL, threshold.PresenceAbsence=NULL,
            species.name=NULL)
        models$TrainData <- TrainData
        models$p <- p
        models$a <- a
        models$MAXENT.a <- MAXENT.a
        if (no.tests==F) {models$pt <- pt}
        if (no.tests==F) {models$at <- at}
        vars <- names(TrainData)
        vars <- vars[which(vars!="pb")]
        models$vars <- vars
        models$factors <- factors
        if(is.null(factors)==F) {models$categories <- categories}
        models$dummy.vars <- dummy.vars
        models$threshold.method <- threshold.method
        models$threshold.sensitivity <- threshold.sensitivity
        if (threshold.PresenceAbsence == T) {
            models$threshold.PresenceAbsence <- TRUE
        }else{
            models$threshold.PresenceAbsence <- FALSE
        }
        models$species.name <- species.name
        if (no.tests == F) {models$TestData <- TestData}
    }else{
        models <- NULL
    }
#
# Data frames for distance-based methods, SVME and GLMNET
    TrainData.vars <- TrainData[,names(TrainData) != "pb", drop=F]
    assign("TrainData.vars", TrainData.vars, envir=.BiodiversityR)
    TrainData.numvars <- TrainData.vars
    if (is.null(factors) == F) {
        for (i in 1:length(factors)) {
            TrainData.numvars <- TrainData.numvars[, which(names(TrainData.numvars) != factors[i]), drop=F]
        }
    }
    assign("TrainData.numvars", TrainData.numvars, envir=.BiodiversityR)
    TrainData.pres <- TrainData[TrainData[,"pb"]==1,,drop=F]
    TrainData.pres <- TrainData.pres[,names(TrainData.pres) != "pb", drop=F]
    assign("TrainData.pres", TrainData.pres, envir=.BiodiversityR)
    var.names <- names(TrainData.pres)
    if (no.tests == F) {
        TestData.vars <- TestData[,names(TestData) != "pb", drop=F]
        assign("TestData.vars", TestData.vars, envir=.BiodiversityR)
        TestData.numvars <- TestData.vars
        if (is.null(factors) == F) {
            for (i in 1:length(factors)) {
                TestData.numvars <- TestData.numvars[, which(names(TestData.numvars) != factors[i]), drop=F]
            }
        }
        assign("TestData.numvars", TestData.numvars, envir=.BiodiversityR)
    }
#
# make MAXENT.TrainData
    if (ws["MAXENT"] > 0  || ws["MAXLIKE"] > 0) {
        if (is.null(MAXENT.a)==T) {
            MAXENT.a <- dismo::randomPoints(x[[1]], n=MAXENT.an, p=p, excludep=T)       
        }
        colnames(MAXENT.a) <- colnames(p)
        if (length(names(x)) == 1) {
            xdouble <- raster::stack(x, x)
            MAXENT.TrainData <- dismo::prepareData(x, p, b=MAXENT.a, factors=factors, xy=FALSE)
            MAXENT.TrainData <- MAXENT.TrainData[, -3]
            names(MAXENT.TrainData)[2] <- names(x)
        }else{
            MAXENT.TrainData <- dismo::prepareData(x, p, b=MAXENT.a, factors=factors, xy=FALSE)
        }
        if(any(is.na(MAXENT.TrainData[MAXENT.TrainData[,"pb"]==0,]))) {
            cat(paste("\n", "WARNING: background locations with missing data removed from MAXENT calibration data","\n\n",sep = ""))
        }
        TrainValid <- complete.cases(MAXENT.TrainData[MAXENT.TrainData[,"pb"]==0,])
        MAXENT.a <- MAXENT.a[TrainValid,]
        if (length(names(x)) == 1) {
            xdouble <- raster::stack(x, x)
            MAXENT.TrainData <- dismo::prepareData(x=xdouble, p, b=MAXENT.a, factors=factors, xy=FALSE)
            MAXENT.TrainData <- MAXENT.TrainData[, -3]
            names(MAXENT.TrainData)[2] <- names(x)
        }else{
            MAXENT.TrainData <- dismo::prepareData(x, p, b=MAXENT.a, factors=factors, xy=FALSE)
        }
        MAXENT.pa <- as.vector(MAXENT.TrainData[ , "pb"])
        MAXENT.TrainData <- MAXENT.TrainData[, which(names(MAXENT.TrainData) != "pb"), drop=F]
        cat(paste("\n", "Summary of Training data set used for calibration of MAXENT or MAXLIKE model (rows: ", nrow(MAXENT.TrainData),  ", presence locations: ", sum(MAXENT.pa), ")\n", sep = ""))
        print(summary(MAXENT.TrainData))
        assign("MAXENT.TrainData", MAXENT.TrainData, envir=.BiodiversityR)
        assign("MAXENT.pa", MAXENT.pa, envir=.BiodiversityR)
    }  
#
    newVIF <- NULL
    if (VIF == T  && length(names(TrainData.vars)) > 1 ) {
#        if (! require(car)) {stop("Please install the car package")}
# only use background data
        TrainDataNum <- TrainData[TrainData[,"pb"]==0,]
        LM.formula <- ensemble.formulae(TrainData, factors=factors)$RF.formula
# create possible response
        TrainDataNum[,"pb"] <- mean(as.numeric(TrainDataNum[,2]))
        vifresult <- NULL
        if (CATCH.OFF == F) {
            tryCatch(vifresult <- car::vif(lm(formula=LM.formula, data=TrainDataNum)),
                error= function(err) {print(paste("\n", "WARNING: VIF (package: car) evaluation failed", "\n", sep=""))},
                        silent=F)
        }else{
            vifresult <- car::vif(lm(formula=LM.formula, data=TrainDataNum))
        }
        if (is.null(vifresult) == F) {
            cat(paste("\n", "Variance inflation (package: car)", "\n", sep = ""))        
            print(vifresult)
        }
        cat(paste("\n", "VIF directly calculated from linear model with focal numeric variable as response", "\n", sep = ""))
        TrainDataNum <- TrainDataNum[,names(TrainDataNum)!="pb"]
        varnames <- names(TrainDataNum)
        newVIF <- numeric(length=length(varnames))
        newVIF[] <- NA
        names(newVIF) <- varnames
        for (i in 1:length(varnames)) {
            response.name <- varnames[i]
            explan.names <- varnames[-i]
            if ((response.name  %in% factors) == F) {
                LM.formula <- as.formula(paste(response.name, "~", paste(explan.names, collapse="+"), sep=""))
                newVIF[i] <- summary(lm(formula=LM.formula, data=TrainDataNum))$r.squared
            }
        }
        newVIF <- 1/(1-newVIF)
        newVIF <- sort(newVIF, decreasing=T, na.last=T)
        print(newVIF)
    }
    if (COR == T) {
        TrainDataNum <- TrainData[, names(TrainData) != "pb", drop=F]
        if(is.null(factors)) {
            for (i in 1:length(factors)) {
                TrainDataNum <- TrainDataNum[, names(TrainDataNum) != factors[i], drop=F]            
            }
        }
        corresult <- cor(TrainDataNum)
        corresult <- round(100*corresult, digits=2)
        cat(paste("\n", "Correlation between numeric variables (as percentage)", "\n", sep = ""))        
        print(corresult)
    }
#
    modelresults <- data.frame(array(dim=c(nrow(TrainData), length(ws)+1), 0))
    names(modelresults) <- c("MAXENT", "MAXLIKE", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", 
        "MGCVFIX", "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "GLMNET", 
        "BIOCLIM.O", "BIOCLIM", "DOMAIN", "MAHAL", "MAHAL01", "ENSEMBLE")
    TrainData <- cbind(TrainData, modelresults)
    assign("TrainData", TrainData, envir=.BiodiversityR)
    if (no.tests == F) {
        modelresults <- data.frame(array(dim=c(nrow(TestData), length(ws)+1), 0))
        names(modelresults) <- c("MAXENT", "MAXLIKE", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", 
            "MGCVFIX", "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "GLMNET", 
            "BIOCLIM.O", "BIOCLIM", "DOMAIN", "MAHAL", "MAHAL01", "ENSEMBLE")
        TestData <- cbind(TestData, modelresults)
        assign("TestData", TestData, envir=.BiodiversityR)
    }
    weights <- as.numeric(array(dim=length(ws), 0))
    names(weights) <- c("MAXENT", "MAXLIKE", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX", 
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "GLMNET", "BIOCLIM.O", "BIOCLIM", "DOMAIN", "MAHAL", "MAHAL01")
#
    if(evaluations.keep==T) {
        evaluations <- list(MAXENT.C=NULL, MAXENT.T=NULL, MAXLIKE.C=NULL, MAXLIKE.T=NULL,
            GBM.trees=NULL, GBM.C=NULL, GBM.T=NULL, GBMSTEP.trees=NULL, GBMSTEP.C=NULL, GBMSTEP.T=NULL, 
            RF.C=NULL, RF.T=NULL, GLM.C=NULL, GLM.T=NULL, GLMS.C=NULL, GLMS.T=NULL, 
            GAM.C=NULL, GAM.T=NULL, GAMS.C=NULL, GAMS.T=NULL, MGCV.C=NULL, MGCV.T=NULL, MGCVF.C=NULL, MGCVF.T=NULL,
            EARTH.C=NULL, EARTH.T=NULL, RPART.C=NULL, RPART.T=NULL,
            NNET.C=NULL, NNET.T=NULL, FDA.C=NULL, FDA.T=NULL, SVM.C=NULL, SVM.T=NULL, SVME.C=NULL, SVME.T=NULL, GLMNET.C=NULL, GLMNET.T=NULL,
            BIOCLIM.O.C=NULL, BIOCLIM.O.T=NULL, BIOCLIM.C=NULL, BIOCLIM.T=NULL, DOMAIN.C=NULL, DOMAIN.T=NULL, MAHAL.C=NULL, MAHAL.T=NULL, MAHAL01.C=NULL, MAHAL01.T=NULL,
            ENSEMBLE.C=NULL, ENSEMBLE.T=NULL, STRATEGY.weights=NULL,
            TrainData=NULL, TestData=NULL, MAXENT.a=NULL,
            factors=NULL, dummy.vars=NULL)
        evaluations$factors <- factors
        evaluations$dummy.vars <- dummy.vars
    }else{
        evaluations <- NULL
    }
#
    Yweights1 <- Yweights
    if (Yweights == "BIOMOD") {
        #have equal weight of presence vs. background
        Yweights1 <- numeric(length = nrow(TrainData))
        pressum <- sum(TrainData[,"pb"]==1)
        abssum <- sum(TrainData[,"pb"]==0)
        Yweights1[which(TrainData[, "pb"] == 1)] <- 1
        Yweights1[which(TrainData[, "pb"] == 0)] <- pressum/abssum
    }
    if (Yweights == "equal") {
        Yweights1 <- numeric(length = nrow(TrainData))
        Yweights1[] <- 1
    }
    assign("Yweights1", Yweights1, envir=.BiodiversityR)
#
# prepare for calculation of deviance
    obs1 <- TrainData[, "pb"]
#
# count models
    mc <- 0
#
# Different modelling algorithms
#
    if(ws["MAXENT"] > 0 && length(names(MAXENT.TrainData)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "MAXENT model can therefore not be calibrated", "\n", sep = ""))
        ws["MAXENT"] <- weights["MAXENT"] <- 0
    }
    if (ws["MAXENT"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Maximum entropy algorithm (package: dismo)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        # Put the file 'maxent.jar' in the 'java' folder of dismo
        # the file 'maxent.jar' can be obtained from from http://www.cs.princeton.edu/~schapire/maxent/.
        jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
        if(is.null(MAXENT.OLD) == T) {
            if (CATCH.OFF == F) {
                tryCatch(results <- dismo::maxent(x=MAXENT.TrainData, p=MAXENT.pa, factors=factors, path=MAXENT.path),
                    error= function(err) {print(paste("MAXENT calibration failed"))},
                    silent=F)
            }else{
                results <- dismo::maxent(x=MAXENT.TrainData, p=MAXENT.pa, factors=factors, path=MAXENT.path)
            }
        }else{ 
            results <- MAXENT.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data of other algorithms","\n\n", sep = ""))
            TrainData[,"MAXENT"] <- dismo::predict(object=results, x=TrainData.vars)
            if (PROBIT == T) {
                if(is.null(MAXENT.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ MAXENT"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- MAXENT.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n", sep = ""))
                TrainData[,"MAXENT.step1"] <- TrainData[,"MAXENT"]
                TrainData[,"MAXENT"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "MAXENT"] < 0), "MAXENT"] <- 0
                TrainData[which(TrainData[, "MAXENT"] > 1), "MAXENT"] <- 1               
            }
            pred1 <- TrainData[, "MAXENT"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"MAXENT"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"MAXENT"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
#            thresholds["MAXENT"] <- threshold(eval1, sensitivity=threshold.sensitivity)[[threshold.method]]
            thresholds["MAXENT"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["MAXENT"]))
            weights["MAXENT"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["MAXENT"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n", sep = ""))
                TestData[,"MAXENT"] <- dismo::predict(object=results, x=TestData.vars)
                if (PROBIT == T) {
                    TestData[,"MAXENT.step1"] <- TestData[,"MAXENT"]
                    TestData[,"MAXENT"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }else{
                    TestData[which(TestData[, "MAXENT"] < 0), "MAXENT"] <- 0
                    TestData[which(TestData[, "MAXENT"] > 1), "MAXENT"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1,"MAXENT"]
                TestAbs <- TestData[TestData[,"pb"]==0,"MAXENT"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["MAXENT"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["MAXENT"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: MAXENT evaluation failed","\n\n",sep = ""))
                    ws["MAXENT"] <- 0
                    weights["MAXENT"] <- 0
                    TrainData[,"MAXENT"] <- 0
                    TestData[,"MAXENT"] <- 0
                }
            }
            if(evaluations.keep ==T) {
                evaluations$MAXENT.C <- eval1
                evaluations$MAXENT.T <- eval2
            }
            if (models.keep==T) {
                models$MAXENT <- results
                models$MAXENT.PROBIT <- results2
            }
        }else{ 
            cat(paste("\n", "WARNING: MAXENT calibration failed", "\n", "\n"))
            ws["MAXENT"] <- weights["MAXENT"] <- 0
            TrainData[,"MAXENT"] <- 0
            if (no.tests == F) {TestData[,"MAXENT"] <- 0}
        }
    }
    if(ws["MAXLIKE"] > 0 && length(names(x)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "MAXLIKE model can therefore not be calibrated", "\n", sep = ""))
        ws["MAXLIKE"] <- weights["MAXLIKE"]<- 0
    }
    if (ws["MAXLIKE"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Maxlike algorithm (package: maxlike)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        MAXLIKE.x <- MAXENT.TrainData[MAXENT.pa == 1, ]
        MAXLIKE.z <- MAXENT.TrainData[MAXENT.pa == 0, ]
        if(is.null(MAXLIKE.OLD) == T) {
            if (CATCH.OFF == F) {
                tryCatch(results <- maxlike::maxlike(formula=MAXLIKE.formula, rasters=NULL, points=NULL, x=MAXLIKE.x, z=MAXLIKE.z, 
                    method=MAXLIKE.method, control=list(maxit=maxit)),
                    error= function(err) {print(paste("MAXLIKE calibration failed"))},
                    silent=F)
            }else{
                results <- maxlike::maxlike(formula=MAXLIKE.formula, rasters=NULL, points=NULL, x=MAXLIKE.x, z=MAXLIKE.z, 
                    method=MAXLIKE.method, control=list(maxit=maxit))
            }
        }else{ 
            results <- MAXLIKE.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data of other algorithms", "\n\n", sep = ""))
            TrainData[, "MAXLIKE"] <- predict(results, newdata=TrainData.numvars)
            if (PROBIT == T) {
                if(is.null(MAXLIKE.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ MAXLIKE"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- MAXLIKE.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n", sep = ""))
                TrainData[,"MAXLIKE.step1"] <- TrainData[,"MAXLIKE"]
                TrainData[,"MAXLIKE"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "MAXLIKE"] < 0), "MAXLIKE"] <- 0
                TrainData[which(TrainData[, "MAXLIKE"] > 1), "MAXLIKE"] <- 1               
            }   
            pred1 <- TrainData[, "MAXLIKE"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1, "MAXLIKE"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0, "MAXLIKE"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
#            thresholds["MAXLIKE"] <- threshold(eval1, sensitivity=threshold.sensitivity)[[threshold.method]]
            thresholds["MAXLIKE"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["MAXLIKE"]))
            weights["MAXLIKE"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["MAXLIKE"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n", sep = ""))
                TestData[, "MAXLIKE"] <- predict(results, newdata=TestData.numvars)
                if (PROBIT == T) {
                    TestData[,"MAXLIKE.step1"] <- TestData[,"MAXLIKE"]
                    TestData[,"MAXLIKE"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }else{
                    TestData[which(TestData[, "MAXLIKE"] < 0), "MAXLIKE"] <- 0
                    TestData[which(TestData[, "MAXLIKE"] > 1), "MAXLIKE"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1, "MAXLIKE"]
                TestAbs <- TestData[TestData[,"pb"]==0, "MAXLIKE"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["MAXLIKE"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["MAXLIKE"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: MAXLIKE evaluation failed","\n\n",sep = ""))
                    ws["MAXLIKE"] <- 0
                    weights["MAXLIKE"] <- 0
                    TrainData[,"MAXLIKE"] <- 0
                    TestData[,"MAXLIKE"] <- 0
                }
            }
            if(evaluations.keep ==T) {
                evaluations$MAXLIKE.C <- eval1
                evaluations$MAXLIKE.T <- eval2
            }
            if (models.keep==T) {
                models$MAXLIKE <- results
                models$MAXLIKE.PROBIT <- results2
                models$formulae$MAXLIKE.formula <- MAXLIKE.formula
            }
        }else{ 
            cat(paste("\n", "WARNING: MAXLIKE calibration failed", "\n", "\n"))
            ws["MAXLIKE"] <- weights["MAXLIKE"] <- 0
            TrainData[,"MAXLIKE"] <- 0
            if (no.tests == F) {TestData[,"MAXLIKE"] <- 0}
        }
    }
    if(ws["GBM"] > 0 && length(names(TrainData.vars)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "GBM model can therefore not be calibrated", "\n", sep = ""))
        ws["GBM"] <- weights["GBM"] <- 0
    }
    if (ws["GBM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Generalized boosted regression modeling (package: gbm) \n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(GBM.OLD) == T) {
            if (CATCH.OFF == F) {    
                tryCatch(results <- gbm::gbm(formula=GBM.formula, data=TrainData, weights=Yweights1, distribution="bernoulli", 
                    interaction.depth=7, shrinkage=0.001, bag.fraction=0.5, train.fraction=1, 
                    n.trees=GBM.n.trees, verbose=F, cv.folds=5),
                error= function(err) {print(paste("GBM calibration failed"))},
                silent=F)
            }else{
                results <- gbm::gbm(formula=GBM.formula, data=TrainData, weights=Yweights1, distribution="bernoulli", 
                    interaction.depth=7, shrinkage=0.001, bag.fraction=0.5, train.fraction=1, 
                    n.trees=GBM.n.trees, verbose=F, cv.folds=5)
            }
        }else{ 
            results <- GBM.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"GBM"] <- gbm::predict.gbm(object=results, newdata=TrainData.vars, n.trees=GBM.n.trees, type="response")
            if (PROBIT == T) {
                if(is.null(GBM.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ GBM"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- GBM.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)", "\n\n", sep = ""))
                TrainData[,"GBM.step1"] <- TrainData[,"GBM"]
                TrainData[,"GBM"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "GBM"] < 0), "GBM"] <- 0
                TrainData[which(TrainData[, "GBM"] > 1), "GBM"] <- 1               
            }   
            pred1 <- TrainData[, "GBM"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GBM"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GBM"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["GBM"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["GBM"]))
            weights["GBM"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["GBM"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"GBM"] <- gbm::predict.gbm(object=results, newdata=TestData.vars, n.trees=GBM.n.trees, type="response")
                if (PROBIT == T) {
                    TestData[,"GBM.step1"] <- TestData[,"GBM"]
                    TestData[,"GBM"] <- predict.glm(object=results2, newdata=TestData, type="response")
               }else{
                    TestData[which(TestData[, "GBM"] < 0), "GBM"] <- 0
                    TestData[which(TestData[, "GBM"] > 1), "GBM"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1,"GBM"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GBM"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["GBM"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["GBM"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: GBM evaluation failed","\n\n",sep = ""))
                    ws["GBM"] <- 0
                    weights["GBM"] <- 0
                    TrainData[,"GBM"] <- 0
                    TestData[,"GBM"] <- 0
                }
            }
            if(evaluations.keep ==T) {
                evaluations$GBM.trees <- results$n.trees 
                evaluations$GBM.C <- eval1
                evaluations$GBM.T <- eval2
            }
            if (models.keep==T) {
                models$GBM <- results
                models$GBM.PROBIT <- results2
                models$formulae$GBM.formula <- GBM.formula
            }
        }else{ 
            cat(paste("\n", "WARNING: GBM calibration failed", "\n", "\n"))
            ws["GBM"] <- weights["GBM"] <- 0
            TrainData[,"GBM"] <- 0
            if (no.tests == F) {TestData[,"GBM"] <- 0}
        }
    }
    if(ws["GBMSTEP"] > 0 && length(names(TrainData.vars)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "GBMSTEP model can therefore not be calibrated", "\n", sep = ""))
        ws["GBMSTEP"] <- weights["GBMSTEP"]<- 0
    }
    if (ws["GBMSTEP"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". gbm step algorithm (package: dismo)\n", sep=""))
#        require(gbm, quietly=T)        
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(GBMSTEP.OLD) == T) {
            if (CATCH.OFF == F) {   
                tryCatch(results <- dismo::gbm.step(data=TrainData, gbm.y=1, gbm.x=GBMSTEP.gbm.x, family="bernoulli",
                    site.weights=Yweights1, tree.complexity=GBMSTEP.tree.complexity, learning.rate = GBMSTEP.learning.rate, 
                    bag.fraction=GBMSTEP.bag.fraction, step.size=GBMSTEP.step.size, verbose=F, silent=F, plot.main=F),
                error= function(err) {print(paste("stepwise GBM calibration failed"))},
                silent=F)
            }else{
                results <- dismo::gbm.step(data=TrainData, gbm.y=1, gbm.x=GBMSTEP.gbm.x, family="bernoulli",
                    site.weights=Yweights1, tree.complexity=GBMSTEP.tree.complexity, learning.rate = GBMSTEP.learning.rate, 
                    bag.fraction=GBMSTEP.bag.fraction, step.size=GBMSTEP.step.size, verbose=F, silent=F, plot.main=F)
            }
        }else{ 
            results <- GBMSTEP.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "stepwise GBM trees (target > 1000)", "\n", sep = ""))
            print(results$n.trees)
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"GBMSTEP"] <- gbm::predict.gbm(object=results, newdata=TrainData.vars, n.trees=GBM.n.trees, type="response")
            if (PROBIT == T) {
                if(is.null(GBMSTEP.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ GBMSTEP"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- GBMSTEP.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"GBMSTEP.step1"] <- TrainData[,"GBMSTEP"]
                TrainData[,"GBMSTEP"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "GBMSTEP"] < 0), "GBMSTEP"] <- 0
                TrainData[which(TrainData[, "GBMSTEP"] > 1), "GBMSTEP"] <- 1               
            }   
            pred1 <- TrainData[, "GBMSTEP"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GBMSTEP"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GBMSTEP"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["GBMSTEP"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["GBMSTEP"]))
            weights["GBMSTEP"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["GBMSTEP"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"GBMSTEP"] <- gbm::predict.gbm(object=results, newdata=TestData.vars, n.trees=GBM.n.trees, type="response")
                if (PROBIT == T) {
                    TestData[,"GBMSTEP.step1"] <- TestData[,"GBMSTEP"]
                    TestData[,"GBMSTEP"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }else{
                    TestData[which(TestData[, "GBMSTEP"] < 0), "GBMSTEP"] <- 0
                    TestData[which(TestData[, "GBMSTEP"] > 1), "GBMSTEP"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1,"GBMSTEP"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GBMSTEP"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["GBMSTEP"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["GBMSTEP"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: stepwise GBM evaluation failed","\n\n",sep = ""))
                    ws["GBMSTEP"] <- 0
                    weights["GBMSTEP"] <- 0
                    TrainData[,"GBMSTEP"] <- 0
                    TestData[,"GBMSTEP"] <- 0
                }
            }
            if(evaluations.keep ==T) {
                evaluations$GBMSTEP.trees <- results$n.trees 
                evaluations$GBMSTEP.C <- eval1
                evaluations$GBMSTEP.T <- eval2
            }
            if (models.keep==T) {
                models$GBMSTEP <- results
                models$GBMSTEP.PROBIT <- results2
            }
        }else{ 
            cat(paste("\n", "WARNING: stepwise GBM calibration failed", "\n", "\n"))
            ws["GBMSTEP"] <- weights["GBMSTEP"] <- 0
            TrainData[,"GBMSTEP"] <- 0
            if (no.tests == F) {TestData[,"GBMSTEP"] <- 0}
        }
    }
    if(ws["RF"] > 0 && length(names(TrainData.vars)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "RF model can therefore not be calibrated", "\n", sep = ""))
        ws["RF"] <- weights["RF"] <- 0
    }
    if (ws["RF"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Random forest algorithm (package: randomForest)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(RF.OLD) == T) {
            if (CATCH.OFF == F) {     
                tryCatch(results <- randomForest::randomForest(formula=RF.formula, ntree=RF.ntree, mtry=RF.mtry, data=TrainData, na.action=na.omit),
                    error= function(err) {print(paste("random forest calibration failed"))},
                    silent=F)
            }else{
                results <- randomForest::randomForest(formula=RF.formula, ntree=RF.ntree, mtry=RF.mtry, data=TrainData, na.action=na.omit)
            }
        }else{ 
            results <- RF.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"RF"] <- predict.RF(object=results, newdata=TrainData.vars)
            if (PROBIT == T) {
                if(is.null(RF.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ RF"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- RF.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"RF.step1"] <- TrainData[,"RF"]
                TrainData[,"RF"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "RF"] < 0), "RF"] <- 0
                TrainData[which(TrainData[, "RF"] > 1), "RF"] <- 1               
            }   
            pred1 <- TrainData[, "RF"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"RF"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"RF"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["RF"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["RF"]))
            weights["RF"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["RF"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"RF"] <- predict.RF(object=results, newdata=TestData.vars)
                if (PROBIT == T) {
                    TestData[,"RF.step1"] <- TestData[,"RF"]
                    TestData[,"RF"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }else{
                    TestData[which(TestData[, "RF"] < 0), "RF"] <- 0
                    TestData[which(TestData[, "RF"] > 1), "RF"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1,"RF"]
                TestAbs <- TestData[TestData[,"pb"]==0,"RF"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["RF"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["RF"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: random forest evaluation failed","\n\n",sep = ""))
                    ws["RF"] <- 0
                    weights["RF"] <- 0
                    TrainData[,"RF"] <- 0
                    TestData[,"RF"] <- 0
                }
            }
            if(evaluations.keep ==T) {
                evaluations$RF.C <- eval1
                evaluations$RF.T <- eval2
            }
            if (models.keep==T) {
                models$RF <- results
                models$RF.PROBIT <- results2
                models$formulae$RF.formula <- RF.formula
            }
        }else{ 
            cat(paste("\n", "WARNING: random forest calibration failed", "\n", "\n"))
            ws["RF"] <- weights["RF"] <- 0
            TrainData[,"RF"] <- 0
            if (no.tests == F) {TestData[,"RF"] <- 0}
        }
    }
    if(ws["GLM"] > 0 && length(names(TrainData.vars)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "GLM model can therefore not be calibrated", "\n", sep = ""))
        ws["GLM"] <- weights["GLM"] <- 0
    }
    if (ws["GLM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Generalized Linear Model \n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(GLM.OLD) == T) {
            if (CATCH.OFF == F) {
                tryCatch(results <- glm(formula=GLM.formula, family=GLM.family, data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit)),
                    error= function(err) {print(paste("GLM calibration failed"))},
                    silent=F)
            }else{
                results <- glm(formula=GLM.formula, family=GLM.family, data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
            }
        }else{ 
            results <- GLM.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"GLM"] <- predict.glm(object=results, newdata=TrainData.vars, type="response")
            if (PROBIT == T) {
                if(is.null(GLM.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ GLM"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- GLM.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"GLM.step1"] <- TrainData[,"GLM"]
                TrainData[,"GLM"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "GLM"] < 0), "GLM"] <- 0
                TrainData[which(TrainData[, "GLM"] > 1), "GLM"] <- 1               
            }   
            pred1 <- TrainData[, "GLM"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GLM"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GLM"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["GLM"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["GLM"]))
            weights["GLM"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["GLM"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"GLM"] <- predict.glm(object=results, newdata=TestData.vars, type="response")
                if (PROBIT == T) {
                    TestData[,"GLM.step1"] <- TestData[,"GLM"]
                    TestData[,"GLM"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }else{
                    TestData[which(TestData[, "GLM"] < 0), "GLM"] <- 0
                    TestData[which(TestData[, "GLM"] > 1), "GLM"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1,"GLM"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GLM"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["GLM"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["GLM"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: GLM evaluation failed","\n\n",sep = ""))
                    ws["GLM"] <- 0
                    weights["GLM"] <- 0
                    TrainData[,"GLM"] <- 0
                    TestData[,"GLM"] <- 0
                }
            }
            if(evaluations.keep ==T) {
                evaluations$GLM.C <- eval1
                evaluations$GLM.T <- eval2
            }
            if (models.keep==T) {
                models$GLM <- results
                models$GLM.PROBIT <- results2
                models$formulae$GLM.formula <- GLM.formula
            }
        }else{ 
            cat(paste("\n", "WARNING: GLM calibration failed", "\n", "\n"))
            ws["GLM"] <- weights["GLM"] <- 0
            TrainData[,"GLM"] <- 0
            if (no.tests == F) {TestData[,"GLM"] <- 0}
        }
    }
    if(ws["GLMSTEP"] > 0 && length(names(TrainData.vars)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "GLMSTEP model can therefore not be calibrated", "\n", sep = ""))
        ws["GLMSTEP"] <- weights["GLMSTEP"] <-0
    }
    if (ws["GLMSTEP"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Stepwise Generalized Linear Model \n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(GLMSTEP.OLD) == T) {
            if (CATCH.OFF == F) {
                tryCatch(results <- glm(formula=STEP.formula, family=GLM.family, data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit)),
                    error= function(err) {print(paste("first step of stepwise GLM calibration failed"))},
                    silent=F)
            }else{
                results <- glm(formula=STEP.formula, family=GLM.family, data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
            }
            if (CATCH.OFF == F) {
                tryCatch(results2 <- MASS::stepAIC(results, scope=GLMSTEP.scope, direction="both", trace=F, steps=GLMSTEP.steps, k=GLMSTEP.k),
                    error= function(err) {print(paste("stepwise GLM calibration failed"))},
                    silent=F)
            }else{
                results2 <- MASS::stepAIC(results, scope=GLMSTEP.scope, direction="both", trace=F, steps=GLMSTEP.steps, k=GLMSTEP.k)
            }
        }else{ 
            results2 <- GLMSTEP.OLD
        }
        if (is.null(results2) == F) {
            results <- results2
            results2 <- NULL
            cat(paste("\n", "stepwise GLM formula","\n\n",sep = ""))
            print(formula(results))
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"GLMSTEP"] <- predict.glm(object=results, newdata=TrainData.vars, type="response")
            if (PROBIT == T) {
                if(is.null(GLMSTEP.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ GLMSTEP"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- GLMSTEP.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"GLMSTEP.step1"] <- TrainData[,"GLMSTEP"]
                TrainData[,"GLMSTEP"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "GLMSTEP"] < 0), "GLMSTEP"] <- 0
                TrainData[which(TrainData[, "GLMSTEP"] > 1), "GLMSTEP"] <- 1               
            }
            pred1 <- TrainData[, "GLMSTEP"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GLMSTEP"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GLMSTEP"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["GLMSTEP"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["GLMSTEP"]))
            weights["GLMSTEP"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["GLMSTEP"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"GLMSTEP"] <- predict.glm(object=results, newdata=TestData.vars, type="response")
                if (PROBIT == T) {
                    TestData[,"GLMSTEP.step1"] <- TestData[,"GLMSTEP"]
                    TestData[,"GLMSTEP"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }else{
                    TestData[which(TestData[, "GLMSTEP"] < 0), "GLMSTEP"] <- 0
                    TestData[which(TestData[, "GLMSTEP"] > 1), "GLMSTEP"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1,"GLMSTEP"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GLMSTEP"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["GLMSTEP"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["GLMSTEP"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: stepwise GLM evaluation failed","\n\n",sep = ""))
                    ws["GLMSTEP"] <- 0
                    weights["GLMSTEP"] <- 0
                    TrainData[,"GLMSTEP"] <- 0
                    TestData[,"MAXENT"] <- 0
                }
            }
            if(evaluations.keep ==T) {
                evaluations$GLMS.C <- eval1
                evaluations$GLMS.T <- eval2
            }
            if (models.keep==T) {
                models$GLMSTEP <- results
                models$GLMSTEP.PROBIT <- results2
                models$formulae$STEP.formula <- STEP.formula
                models$formulae$GLMSTEP.scope <- GLMSTEP.scope
            }
        }else{ 
            cat(paste("\n", "WARNING: stepwise GLM calibration failed", "\n", "\n"))
            ws["GLMSTEP"] <- weights["GLMSTEP"] <- 0
            TrainData[,"GLMSTEP"] <- 0
            if (no.tests == F) {TestData[,"GLMSTEP"] <- 0}
        }
    }
    if(ws["GAM"] > 0 && length(names(TrainData.vars)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "GAM model can therefore not be calibrated", "\n", sep = ""))
        ws["GAM"] <- weights["GAM"] <- 0
    }
    if (ws["GAM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Generalized Additive Model (package: gam)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(GAM.OLD) == T) {
            if (CATCH.OFF == F) {
                tryCatch(results <- gam::gam(formula=GAM.formula, family=GAM.family, data=TrainData, weights=Yweights1, control=gam::gam.control(maxit=maxit, bf.maxit=50)),
                    error= function(err) {print(paste("GAM (package: gam) calibration failed"))},
                    silent=F)
            }else{
                results <- gam::gam(formula=GAM.formula, family=GAM.family, data=TrainData, weights=Yweights1, control=gam::gam.control(maxit=maxit, bf.maxit=50))
            }
        }else{ 
            results <- GAM.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"GAM"] <- gam::predict.Gam(object=results, newdata=TrainData.vars, type="response")
            if (PROBIT == T) {
                if(is.null(GAM.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ GAM"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- GAM.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"GAM.step1"] <- TrainData[,"GAM"]
                TrainData[,"GAM"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "GAM"] < 0), "GAM"] <- 0
                TrainData[which(TrainData[, "GAM"] > 1), "GAM"] <- 1               
            }   
            pred1 <- TrainData[, "GAM"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GAM"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GAM"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["GAM"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["GAM"]))
            weights["GAM"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["GAM"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"GAM"] <- gam::predict.Gam(object=results, newdata=TestData.vars, type="response")
                if (PROBIT == T) {
                    TestData[,"GAM.step1"] <- TestData[,"GAM"]
                    TestData[,"GAM"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }else{
                    TestData[which(TestData[, "GAM"] < 0), "GAM"] <- 0
                    TestData[which(TestData[, "GAM"] > 1), "GAM"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1,"GAM"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GAM"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["GAM"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["GAM"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: GAM (package: gam) evaluation failed","\n\n",sep = ""))
                    ws["GAM"] <- 0
                    weights["GAM"] <- 0
                    TrainData[,"GAM"] <- 0
                    TestData[,"GAM"] <- 0
                }
            }
            if(evaluations.keep ==T) {
                evaluations$GAM.C <- eval1
                evaluations$GAM.T <- eval2
            }
            if (models.keep==T) {
                models$GAM <- results
                models$GAM.PROBIT <- results2
                models$formulae$GAM.formula <- GAM.formula
            }
        }else{ 
            cat(paste("\n", "WARNING: GAM (package: gam) calibration failed", "\n", "\n"))
            ws["GAM"] <- weights["GAM"] <- 0
            TrainData[,"GAM"] <- 0
            if (no.tests == F) {TestData[,"GAM"] <- 0}
        }
    }
    if(ws["GAMSTEP"] > 0 && length(names(TrainData.vars)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "GAMSTEP model can therefore not be calibrated", "\n", sep = ""))
        ws["GAMSTEP"] <- weights["GAMSTEP"] <- 0
    }
    if (ws["GAMSTEP"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Stepwise Generalized Additive Model (package: gam)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(GAMSTEP.OLD) == T) {
            if (CATCH.OFF == F) {
                tryCatch(results <- gam::gam(formula=STEP.formula, family=GAM.family, data=TrainData, weights=Yweights1, control=gam::gam.control(maxit=maxit, bf.maxit=50)), 
                    error= function(err) {print(paste("first step of stepwise GAM (package: gam) calibration failed"))},
                    silent=F)
            }else{
                results <- gam::gam(formula=STEP.formula, family=GAM.family, data=TrainData, weights=Yweights1, control=gam::gam.control(maxit=maxit, bf.maxit=50))
            }
            assign("TrainData", TrainData, pos=GAMSTEP.pos)
            assign("GAM.family", GAM.family, pos=GAMSTEP.pos)
            assign("maxit", maxit, pos=GAMSTEP.pos)
            if (CATCH.OFF == F) {
                tryCatch(results2 <- gam::step.Gam(results, scope=GAMSTEP.scope, direction="both", trace=F, steps=GAMSTEP.steps), 
                    error= function(err) {print(paste("stepwise GAM (package: gam) calibration failed"))},
                    silent=F)
            }else{
                results2 <- gam::step.Gam(results, scope=GAMSTEP.scope, direction="both", trace=F, steps=GAMSTEP.steps)
            }
            remove(TrainData, pos=GAMSTEP.pos)
            remove(GAM.family, pos=GAMSTEP.pos)
            remove(maxit, pos=GAMSTEP.pos)
        }else{ 
            results2 <- GAMSTEP.OLD
        }
        if (is.null(results2) == F) {
            results <- results2
            results2 <- NULL
            cat(paste("\n", "stepwise GAM formula (gam package)","\n\n",sep = ""))
            print(formula(results))
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"GAMSTEP"] <- gam::predict.Gam(object=results, newdata=TrainData.vars, type="response")
            if (PROBIT == T) {
                if(is.null(GAMSTEP.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ GAMSTEP"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- GAMSTEP.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"GAMSTEP.step1"] <- TrainData[,"GAMSTEP"]
                TrainData[,"GAMSTEP"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "GAMSTEP"] < 0), "GAMSTEP"] <- 0
                TrainData[which(TrainData[, "GAMSTEP"] > 1), "GAMSTEP"] <- 1               
            }   
            pred1 <- TrainData[, "GAMSTEP"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GAMSTEP"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GAMSTEP"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["GAMSTEP"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["GAMSTEP"]))
            weights["GAMSTEP"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["GAMSTEP"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n","\n", sep = ""))
                TestData[,"GAMSTEP"] <- gam::predict.Gam(object=results, newdata=TestData.vars, type="response")
                if (PROBIT == T) {
                    TestData[,"GAMSTEP.step1"] <- TestData[,"GAMSTEP"]
                    TestData[,"GAMSTEP"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }else{
                    TestData[which(TestData[, "GAMSTEP"] < 0), "GAMSTEP"] <- 0
                    TestData[which(TestData[, "GAMSTEP"] > 1), "GAMSTEP"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1,"GAMSTEP"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GAMSTEP"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["GAMSTEP"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["GAMSTEP"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: stepwise GAM (package: gam) evaluation failed","\n\n",sep = ""))
                    ws["GAMSTEP"] <- 0
                    weights["GAMSTEP"] <- 0
                    TrainData[,"GAMSTEP"] <- 0
                    TestData[,"GAMSTEP"] <- 0
                }
            }
            if(evaluations.keep ==T) {
                evaluations$GAMS.C <- eval1
                evaluations$GAMS.T <- eval2
            }
            if (models.keep==T) {
                models$GAMSTEP <- results
                models$GAMSTEP.PROBIT <- results2
                models$formulae$STEP.formula <- STEP.formula
                models$formulae$GAMSTEP.scope <- GAMSTEP.scope
            }
        }else{ 
            cat(paste("\n", "WARNING: stepwise GAM (package: gam) calibration failed", "\n", "\n"))
            ws["GAMSTEP"] <- weights["GAMSTEP"] <- 0
            TrainData[,"GAMSTEP"] <- 0
            if (no.tests == F) {TestData[,"GAMSTEP"] <- 0}
        }
    }
    if(ws["MGCV"] > 0 && length(names(TrainData.vars)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "MGCV model can therefore not be calibrated", "\n", sep = ""))
        ws["MGCV"] <- weights["MGCV"] <- 0
    }
    if (ws["MGCV"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Generalized Additive Model (package: mgcv)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(MGCV.OLD) == T) {
            if (CATCH.OFF == F) {
                tryCatch(results <- mgcv::gam(formula=MGCV.formula, family=GAM.family, data=TrainData, weights=Yweights1, 
                    select=MGCV.select, control=mgcv::gam.control(maxit=maxit)),
                error= function(err) {print(paste("GAM (package: mgcv) calibration failed"))},
                silent=F)
            }else{
                results <- mgcv::gam(formula=MGCV.formula, family=GAM.family, data=TrainData, weights=Yweights1, 
                    select=MGCV.select, control=mgcv::gam.control(maxit=maxit))
            }
        }else{ 
            results <- MGCV.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"MGCV"] <- predict.MGCV(object=results, newdata=TrainData.vars)
            if (PROBIT == T) {
                if(is.null(MGCV.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ MGCV"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- MGCV.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"MGCV.step1"] <- TrainData[,"MGCV"]
                TrainData[,"MGCV"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "MGCV"] < 0), "MGCV"] <- 0
                TrainData[which(TrainData[, "MGCV"] > 1), "MGCV"] <- 1               
            }   
            pred1 <- TrainData[, "MGCV"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"MGCV"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"MGCV"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["MGCV"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["MGCV"]))
            weights["MGCV"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["MGCV"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"MGCV"] <- predict.MGCV(object=results, newdata=TestData.vars)
                if (PROBIT == T) {
                    TestData[,"MGCV.step1"] <- TestData[,"MGCV"]
                    TestData[,"MGCV"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }else{
                    TestData[which(TestData[, "MGCV"] < 0), "MGCV"] <- 0
                    TestData[which(TestData[, "MGCV"] > 1), "MGCV"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1,"MGCV"]
                TestAbs <- TestData[TestData[,"pb"]==0,"MGCV"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["MGCV"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["MGCV"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: GAM (package: mgcv) evaluation failed","\n\n",sep = ""))
                    ws["MGCV"] <- 0
                    weights["MGCV"] <- 0
                    TrainData[,"MGCV"] <- 0
                    TestData[,"MGCV"] <- 0
                }
            }
            if(evaluations.keep ==T) {
                evaluations$MGCV.C <- eval1
                evaluations$MGCV.T <- eval2
            }
            if (models.keep==T) {
                models$MGCV <- results
                models$MGCV.PROBIT <- results2
                models$formulae$MGCV.formula <- MGCV.formula
            }
        }else{ 
            cat(paste("\n", "WARNING: GAM (package: mgcv) calibration failed", "\n", "\n"))
            ws["MGCV"] <- weights["MGCV"] <- 0
            TrainData[,"MGCV"] <- 0
            if (no.tests == F) {TestData[,"MGCV"] <- 0}
        }
    }
    if(ws["MGCVFIX"] > 0 && length(names(TrainData.vars)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "MGCVFIX model can therefore not be calibrated", "\n", sep = ""))
        ws["MGCVFIX"] <- weights["MGCVFIX"] <- 0
    }
    if (ws["MGCVFIX"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". GAM with fixed d.f. regression splines (package: mgcv)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(MGCVFIX.OLD) == T) {
            if (CATCH.OFF == F) {
                tryCatch(results <- mgcv::gam(formula=MGCVFIX.formula, family=GAM.family, data=TrainData, weights=Yweights1, select=FALSE, control=mgcv::gam.control(maxit=maxit)),
                    error= function(err) {print(paste("GAM with fixed d.f. regression splines (package: mgcv) calibration failed"))},
                    silent=F)
            }else{
                results <- mgcv::gam(formula=MGCVFIX.formula, family=GAM.family, data=TrainData, weights=Yweights1, select=FALSE, control=mgcv::gam.control(maxit=maxit))
            }
        }else{ 
            results <- MGCVFIX.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"MGCVFIX"] <- predict.MGCV(object=results, newdata=TrainData.vars)
            if (PROBIT == T) {
                if(is.null(MGCVFIX.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ MGCVFIX"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- MGCVFIX.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"MGCVFIX.step1"] <- TrainData[,"MGCVFIX"]
                TrainData[,"MGCVFIX"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "MGCVFIX"] < 0), "MGCVFIX"] <- 0
                TrainData[which(TrainData[, "MGCVFIX"] > 1), "MGCVFIX"] <- 1               
            }  
            pred1 <- TrainData[, "MGCVFIX"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"MGCVFIX"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"MGCVFIX"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["MGCVFIX"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["MGCVFIX"]))
            weights["MGCVFIX"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["MGCVFIX"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"MGCVFIX"] <- predict.MGCV(object=results, newdata=TestData)
                if (PROBIT == T) {
                    TestData[,"MGCVFIX.step1"] <- TestData[,"MGCVFIX"]
                    TestData[,"MGCVFIX"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }else{
                    TestData[which(TestData[, "MGCVFIX"] < 0), "MGCVFIX"] <- 0
                    TestData[which(TestData[, "MGCVFIX"] > 1), "MGCVFIX"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1,"MGCVFIX"]
                TestAbs <- TestData[TestData[,"pb"]==0,"MGCVFIX"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["MGCVFIX"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["MGCVFIX"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: GAM with fixed d.f. regression splines (package: mgcv) evaluation failed","\n\n",sep = ""))
                    ws["MGCVFIX"] <- 0
                    weights["MGCVFIX"] <- 0
                    TrainData[,"MGCVFIX"] <- 0
                    TestData[,"MGCVFOX"] <- 0
                }
            }
            if(evaluations.keep ==T) {
                evaluations$MGCVF.C <- eval1
                evaluations$MGCVF.T <- eval2
            }
            if (models.keep==T) {
                models$MGCVFIX <- results
                models$MGCVFIX.PROBIT <- results2
                models$formulae$MGCVFIX.formula <- MGCVFIX.formula
            }
        }else{ 
            cat(paste("\n", "WARNING: GAM with fixed d.f. regression splines (package: mgcv) calibration failed", "\n", "\n"))
            ws["MGCVFIX"] <- weights["MGCVFIX"] <- 0
            TrainData[,"MGCVFIX"] <- 0
            if (no.tests == F) {TestData[,"MGCVFIX"] <- 0}
        }
    }
    if(ws["EARTH"] > 0 && length(names(TrainData.numvars)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "MARS (EARTH) model can therefore not be calibrated", "\n", sep = ""))
        ws["EARTH"] <- weights["EARTH"] <- 0
    }
    if (ws["EARTH"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Multivariate Adaptive Regression Splines (package: earth)\n", sep=""))
        if (is.null(factors) == F) {
            cat(paste("\n", "NOTE: factors could not be used as explanatory variables for MARS (maybe consider dummy variables)", "\n", sep=""))
        } 
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(EARTH.OLD) == T) {
            if (CATCH.OFF == F) {
                tryCatch(results <- earth::earth(formula=EARTH.formula, glm=EARTH.glm, data=TrainData, degree=2),
                    error= function(err) {print(paste("MARS (package: earth) calibration failed"))},
                    silent=F)
            }else{
                results <- earth::earth(formula=EARTH.formula, glm=EARTH.glm, data=TrainData, degree=2)
            }
        }else{ 
            results <- EARTH.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"EARTH"] <- predict.EARTH(object=results, newdata=TrainData.vars)
            if (PROBIT == T) {
                if(is.null(EARTH.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ EARTH"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- EARTH.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"EARTH.step1"] <- TrainData[,"EARTH"]
                TrainData[,"EARTH"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "EARTH"] < 0), "EARTH"] <- 0
                TrainData[which(TrainData[, "EARTH"] > 1), "EARTH"] <- 1               
            }   
            pred1 <- TrainData[, "EARTH"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"EARTH"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"EARTH"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["EARTH"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["EARTH"]))
            weights["EARTH"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["EARTH"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"EARTH"] <- predict.EARTH(object=results, newdata=TestData.vars)
                if (PROBIT == T) {
                    TestData[,"EARTH.step1"] <- TestData[,"EARTH"]
                    TestData[,"EARTH"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }else{
                    TestData[which(TestData[, "EARTH"] < 0), "EARTH"] <- 0
                    TestData[which(TestData[, "EARTH"] > 1), "EARTH"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1,"EARTH"]
                TestAbs <- TestData[TestData[,"pb"]==0,"EARTH"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["EARTH"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["EARTH"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: MARS (package: earth) evaluation failed","\n\n",sep = ""))
                    ws["EARTH"] <- 0
                    weights["EARTH"] <- 0
                    TrainData[,"EARTH"] <- 0
                    TestData[,"EARTH"] <- 0
                }
            }
            if(evaluations.keep ==T) {
                evaluations$EARTH.C <- eval1
                evaluations$EARTH.T <- eval2
            }
            if (models.keep==T) {
                models$EARTH <- results
                models$EARTH.PROBIT <- results2
                models$formulae$EARTH.formula <- EARTH.formula
            }
        }else{ 
            cat(paste("\n", "WARNING: MARS (package: earth) calibration failed", "\n", "\n"))
            ws["EARTH"] <- weights["EARTH"] <- 0
            TrainData[,"EARTH"] <- 0
            if (no.tests == F) {TestData[,"EARTH"] <- 0}
        }
    }
    if(ws["RPART"] > 0 && length(names(TrainData.vars)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "RPART model can therefore not be calibrated", "\n", sep = ""))
        ws["RPART"] <- weights["rpart"] <- 0
    }
    if (ws["RPART"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Recursive Partitioning And Regression Trees (package: rpart)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(RPART.OLD) == T) {
            if (CATCH.OFF == F) {
                tryCatch(results <- rpart::rpart(formula=RPART.formula, data=TrainData, weights=Yweights1,
                    control=rpart::rpart.control(xval=RPART.xval, minbucket=5, minsplit=5, cp=0.001, maxdepth=25)),
                    error= function(err) {print(paste("RPART calibration failed"))},
                silent=F)
            }else{
                results <- rpart::rpart(formula=RPART.formula, data=TrainData, weights=Yweights1,
                    control=rpart::rpart.control(xval=RPART.xval, minbucket=5, minsplit=5, cp=0.001, maxdepth=25))
            }
        }else{ 
            results <- RPART.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"RPART"] <- predict(object=results, newdata=TrainData.vars, type="prob")[,2]
            if (PROBIT == T) {
                if(is.null(RPART.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ RPART"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- RPART.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"RPART.step1"] <- TrainData[,"RPART"]
                TrainData[,"RPART"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "RPART"] < 0), "RPART"] <- 0
                TrainData[which(TrainData[, "RPART"] > 1), "RPART"] <- 1               
            }  
            pred1 <- TrainData[, "RPART"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"RPART"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"RPART"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["RPART"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["RPART"]))
            weights["RPART"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["RPART"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"RPART"] <- predict(object=results, newdata=TestData.vars, type="prob")[,2]
                if (PROBIT == T) {
                    TestData[,"RPART.step1"] <- TestData[,"RPART"]
                    TestData[,"RPART"] <- predict.glm(object=results2, newdata=TestData, type="response")
               }else{
                    TestData[which(TestData[, "RPART"] < 0), "RPART"] <- 0
                    TestData[which(TestData[, "RPART"] > 1), "RPART"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1,"RPART"]
                TestAbs <- TestData[TestData[,"pb"]==0,"RPART"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(TestPres) == F && is.null(TestAbs) == F) {
                    eval2 <-  dismo::evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
                    weights["RPART"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["RPART"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: RPART evaluation failed","\n\n",sep = ""))
                    ws["RPART"] <- 0
                    weights["RPART"] <- 0
                    TrainData[,"RPART"] <- 0
                    TestData[,"RPART"] <- 0
                }
            }
            if(evaluations.keep ==T) {
                evaluations$RPART.C <- eval1
                evaluations$RPART.T <- eval2
            }
            if (models.keep==T) {
                models$RPART <- results
                models$RPART.PROBIT <- results2
                models$formulae$RPART.formula <- RPART.formula
            }
        }else{ 
            cat(paste("\n", "WARNING: RPART calibration failed", "\n", "\n"))
            ws["RPART"] <- weights["RPART"] <- 0
            TrainData[,"RPART"] <- 0
            if (no.tests == F) {TestData[,"RPART"] <- 0}
        }
    }
    if(ws["NNET"] > 0 && length(names(TrainData.vars)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "NNET model can therefore not be calibrated", "\n", sep = ""))
        ws["NNET"] <- weights["NNET"] <- 0
    }
    if (ws["NNET"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Artificial Neural Network (package: nnet)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(NNET.OLD) == T) {
            if (CATCH.OFF == F) {
                tryCatch(results <- nnet::nnet(formula=NNET.formula, size=NNET.size, decay=NNET.decay, data=TrainData, weights=Yweights1, 
                    rang=0.1, maxit=maxit, trace=F),
                error= function(err) {print(paste("Artificial Neural Network (package: nnet) calibration failed"))},
                silent=F)
            }else{
                results <- nnet::nnet(formula=NNET.formula, size=NNET.size, decay=NNET.decay, data=TrainData, weights=Yweights1, 
                    rang=0.1, maxit=maxit, trace=F)
            }
        }else{ 
            results <- NNET.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"NNET"] <- predict.NNET(object=results, newdata=TrainData.vars)
            if (PROBIT == T) {
                if(is.null(NNET.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ NNET"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- NNET.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"NNET.step1"] <- TrainData[,"NNET"]
                TrainData[,"NNET"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "NNET"] < 0), "NNET"] <- 0
                TrainData[which(TrainData[, "NNET"] > 1), "NNET"] <- 1               
            }   
            pred1 <- TrainData[, "NNET"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"NNET"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"NNET"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["NNET"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["NNET"]))
            weights["NNET"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["NNET"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"NNET"] <- predict.NNET(object=results, newdata=TestData.vars)
                if (PROBIT == T) {
                    TestData[,"NNET.step1"] <- TestData[,"NNET"]
                    TestData[,"NNET"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }else{
                    TestData[which(TestData[, "NNET"] < 0), "NNET"] <- 0
                    TestData[which(TestData[, "NNET"] > 1), "NNET"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1,"NNET"]
                TestAbs <- TestData[TestData[,"pb"]==0,"NNET"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["NNET"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["NNET"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: Artificial Neural Network (package: nnet) evaluation failed","\n\n",sep = ""))
                    ws["NNET"] <- 0
                    weights["NNET"] <- 0
                    TrainData[,"NNET"] <- 0
                    TestData[,"NNET"] <- 0
                }
            }
            if(evaluations.keep ==T) {
                evaluations$NNET.C <- eval1
                evaluations$NNET.T <- eval2
            }
            if (models.keep==T) {
                models$NNET <- results
                models$NNET.PROBIT <- results2
                models$formulae$NNET.formula <- NNET.formula
            }
        }else{ 
            cat(paste("\n", "WARNING: Artificial Neural Network (package: nnet) calibration failed", "\n", "\n"))
            ws["NNET"] <- weights["NNET"] <- 0
            TrainData[,"NNET"] <- 0
            if (no.tests == F) {TestData[,"NNET"] <- 0}
        }
    }
    if(ws["FDA"] > 0 && length(names(TrainData.vars)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "FDA model can therefore not be calibrated", "\n", sep = ""))
        ws["FDA"] <- weights["FDA"] <- 0
    }
    if (ws["FDA"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Flexible Discriminant Analysis (package: mda)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(FDA.OLD) == T) {
            if (CATCH.OFF == F) {
                tryCatch(results <- mda::fda(formula=FDA.formula, method=mda::mars, data=TrainData, weights=Yweights1),
                    error= function(err) {print(paste("Flexible Discriminant Analysis calibration failed"))},
                    silent=F)
            }else{
                results <- mda::fda(formula=FDA.formula, method=mda::mars, data=TrainData, weights=Yweights1)
            }
        }else{ 
            results <- FDA.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"FDA"] <- predict(object=results, newdata=TrainData.vars, type="posterior")[,2]
            if (PROBIT == T) {
                if(is.null(FDA.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ FDA"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- FDA.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"FDA.step1"] <- TrainData[,"FDA"]
                TrainData[,"FDA"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "FDA"] < 0), "FDA"] <- 0
                TrainData[which(TrainData[, "FDA"] > 1), "FDA"] <- 1               
            } 
            pred1 <- TrainData[, "FDA"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"FDA"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"FDA"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["FDA"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["FDA"]))
            weights["FDA"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["FDA"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"FDA"] <- predict(object=results, newdata=TestData.vars, type="posterior")[,2]
                if (PROBIT == T) {
                    TestData[,"FDA.step1"] <- TestData[,"FDA"]
                    TestData[,"FDA"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }else{
                    TestData[which(TestData[, "FDA"] < 0), "FDA"] <- 0
                    TestData[which(TestData[, "FDA"] > 1), "FDA"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1,"FDA"]
                TestAbs <- TestData[TestData[,"pb"]==0,"FDA"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["FDA"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["FDA"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: Flexible Discriminant Analysis evaluation failed","\n\n",sep = ""))
                    ws["FDA"] <- 0
                    weights["FDA"] <- 0
                    TrainData[,"FDA"] <- 0
                    TestData[,"FDA"] <- 0
                }
            }
            if(evaluations.keep ==T) {
                evaluations$FDA.C <- eval1
                evaluations$FDA.T <- eval2
            }
            if (models.keep==T) {
                models$FDA <- results
                models$FDA.PROBIT <- results2
                models$formulae$FDA.formula <- FDA.formula
            }
        }else{ 
            cat(paste("\n", "WARNING: Flexible Discriminant Analysis calibration failed", "\n", "\n"))
            ws["FDA"] <- weights["FDA"] <- 0
            TrainData[,"FDA"] <- 0
            if (no.tests == F) {TestData[,"FDA"] <- 0}
        }
    }
    if(ws["SVM"] > 0 && length(names(TrainData.vars)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "SVM model can therefore not be calibrated", "\n", sep = ""))
        ws["SVM"] <- weights["SVM"] <- 0
    }
    if (ws["SVM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Support Vector Machines (package: kernlab)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(SVM.OLD) == T) {
            if (CATCH.OFF == F) {
                tryCatch(results <- kernlab::ksvm(SVM.formula, data=TrainData, type="C-svc", prob.model=T),
                    error= function(err) {print(paste("Support Vector Machines (package: kernlab) calibration failed"))},
                    silent=F)
            }else{
                results <- kernlab::ksvm(SVM.formula, data=TrainData, type="C-svc", prob.model=T)
            }
        }else{ 
            results <- SVM.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"SVM"] <- kernlab::predict(object=results, newdata=TrainData.vars, type="probabilities")[,2]
            if (PROBIT == T) {
                if(is.null(SVM.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ SVM"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- SVM.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"SVM.step1"] <- TrainData[,"SVM"]
                TrainData[,"SVM"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "SVM"] < 0), "SVM"] <- 0
                TrainData[which(TrainData[, "SVM"] > 1), "SVM"] <- 1               
            }   
            pred1 <- TrainData[, "SVM"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"SVM"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"SVM"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["SVM"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["SVM"]))
            weights["SVM"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["SVM"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"SVM"] <- kernlab::predict(object=results, newdata=TestData.vars, type="probabilities")[,2]
                if (PROBIT == T) {
                    TestData[,"SVM.step1"] <- TestData[,"SVM"]
                    TestData[,"SVM"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }else{
                    TestData[which(TestData[, "SVM"] < 0), "SVM"] <- 0
                    TestData[which(TestData[, "SVM"] > 1), "SVM"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1,"SVM"]
                TestAbs <- TestData[TestData[,"pb"]==0,"SVM"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["SVM"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["SVM"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: Support Vector Machines (package: kernlab)  evaluation failed","\n\n",sep = ""))
                    ws["SVM"] <- 0
                    weights["SVM"] <- 0
                    TrainData[,"SVM"] <- 0
                    TestData[,"SVM"] <- 0
                }
            }
            if(evaluations.keep ==T) {
                evaluations$SVM.C <- eval1
                evaluations$SVM.T <- eval2
            }
            if (models.keep==T) {
                models$SVM <- results
                models$SVM.PROBIT <- results2
                models$formulae$SVM.formula <- SVM.formula
            }
        }else{ 
            cat(paste("\n", "WARNING: Support Vector Machines (package: kernlab) calibration failed", "\n", "\n"))
            ws["SVM"] <- weights["SVM"] <- 0
            TrainData[,"SVM"] <- 0
            if (no.tests == F) {TestData[,"SVM"] <- 0}
        }
    }
    if(ws["SVME"] > 0 && length(names(TrainData.vars)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "SVME model can therefore not be calibrated", "\n", sep = ""))
        ws["SVME"] <- weights["SVME"] <- 0
    }
    if (ws["SVME"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Support Vector Machines (package: e1071)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(SVME.OLD) == T) {
            if (CATCH.OFF == F) {
                tryCatch(results <- e1071::svm(SVME.formula, data=TrainData, type="C-classification", kernel="polynomial", degree=3, probability=TRUE),
                    error= function(err) {print(paste("Support Vector Machines (package: e1071) calibration failed"))},
                    silent=F)
            }else{
                results <- e1071::svm(SVME.formula, data=TrainData, type="C-classification", kernel="polynomial", degree=3, probability=TRUE)
            }
        }else{ 
            results <- SVME.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"SVME"] <- predict.SVME(model=results, newdata=TrainData.vars)
            if (PROBIT == T) {
                if(is.null(SVME.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ SVME"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- SVME.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"SVME.step1"] <- TrainData[,"SVME"]
                TrainData[,"SVME"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "SVME"] < 0), "SVME"] <- 0
                TrainData[which(TrainData[, "SVME"] > 1), "SVME"] <- 1               
            }  
            pred1 <- TrainData[, "SVME"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"SVME"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"SVME"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["SVME"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["SVME"]))
            weights["SVME"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["SVME"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"SVME"] <- predict.SVME(model=results, newdata=TestData.vars)
                if (PROBIT == T) {
                    TestData[,"SVME.step1"] <- TestData[,"SVME"]
                    TestData[,"SVME"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }else{
                    TestData[which(TestData[, "SVME"] < 0), "SVME"] <- 0
                    TestData[which(TestData[, "SVME"] > 1), "SVME"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1,"SVME"]
                TestAbs <- TestData[TestData[,"pb"]==0,"SVME"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["SVME"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["SVME"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: Support Vector Machines (package: e1071) evaluation failed","\n\n",sep = ""))
                    ws["SVME"] <- 0
                    weights["SVME"] <- 0
                    TrainData[,"SVME"] <- 0
                    TestData[,"SVME"] <- 0
                }
            }
            if(evaluations.keep ==T) {
                evaluations$SVME.C <- eval1
                evaluations$SVME.T <- eval2
            }
            if (models.keep==T) {
                models$SVME <- results
                models$SVME.PROBIT <- results2
                models$formulae$SVME.formula <- SVME.formula
            }
        }else{ 
            cat(paste("\n", "WARNING: Support Vector Machines (package: e1071) calibration failed", "\n", "\n"))
            ws["SVME"] <- weights["SVME"] <- 0
            TrainData[,"SVME"] <- 0
            if (no.tests == F) {TestData[,"SVME"] <- 0}
        }
    }
    if(ws["GLMNET"] > 0 && length(names(TrainData.numvars)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "GLMNET model can therefore not be calibrated", "\n", sep = ""))
        ws["GLMNET"] <- weights["GLMNET"] <- 0
    }
    if (ws["GLMNET"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". GLM with lasso or elasticnet regularization (package: glmnet)\n", sep=""))
        if (is.null(factors) == F) {
            cat(paste("\n", "NOTE: factors could not be used as explanatory variables for GLMNET (maybe consider dummy variables)", "\n", sep=""))
        } 
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(GLMNET.OLD) == T) {
            if (CATCH.OFF == F) {
                tryCatch(results <- glmnet::glmnet(x=as.matrix(TrainData.numvars), y=TrainData[, "pb"], family="binomial", weights=Yweights1, nlambda=GLMNET.nlambda),
                    error= function(err) {print(paste("GLMNET calibration failed"))},
                    silent=F)
            }else{
                results <- glmnet::glmnet(x=as.matrix(TrainData.numvars), y=TrainData[, "pb"], family="binomial", weights=Yweights1, nlambda=GLMNET.nlambda)
            }
        }else{ 
            results <- GLMNET.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"GLMNET"] <- predict.GLMNET(model=results, newdata=TrainData.numvars, GLMNET.class=GLMNET.class)
            if (PROBIT == T) {
                if(is.null(GLMNET.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ GLMNET"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- GLMNET.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"GLMNET.step1"] <- TrainData[,"GLMNET"]
                TrainData[,"GLMNET"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "GLMNET"] < 0), "GLMNET"] <- 0
                TrainData[which(TrainData[, "GLMNET"] > 1), "GLMNET"] <- 1               
            } 
            pred1 <- TrainData[, "GLMNET"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GLMNET"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GLMNET"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["GLMNET"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["GLMNET"]))
            weights["GLMNET"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["GLMNET"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"GLMNET"] <- predict.GLMNET(model=results, newdata=TestData.numvars, GLMNET.class=GLMNET.class)
                if (PROBIT == T) {
                    TestData[,"GLMNET.step1"] <- TestData[,"GLMNET"]
                    TestData[,"GLMNET"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }else{
                    TestData[which(TestData[, "GLMNET"] < 0), "GLMNET"] <- 0
                    TestData[which(TestData[, "GLMNET"] > 1), "GLMNET"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1,"GLMNET"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GLMNET"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["GLMNET"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["GLMNET"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: GLMNET evaluation failed","\n\n",sep = ""))
                    ws["GLMNET"] <- 0
                    weights["GLMNET"] <- 0
                    TrainData[,"GLMNET"] <- 0
                    TestData[,"GLMNET"] <- 0
                }
            }
            if(evaluations.keep == T) {
                evaluations$GLMNET.C <- eval1
                evaluations$GLMNET.T <- eval2
            }
            if (models.keep==T) {
                models$GLMNET <- results
                models$GLMNET.PROBIT <- results2
                if (GLMNET.class == F) {models$formulae$GLMNET.class <- FALSE}
                if (GLMNET.class == T) {models$formulae$GLMNET.class <- TRUE}
            }
        }else{ 
            cat(paste("\n", "WARNING: GLMNET calibration failed", "\n", "\n"))
            ws["GLMNET"] <- weights["GLMNET"] <- 0
            TrainData[,"GLMNET"] <- 0
            if (no.tests == F) {TestData[,"GLMNET"] <- 0}
        }
    }
    if(ws["BIOCLIM.O"] > 0 && length(names(TrainData.pres)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "BIOCLIM.O model can therefore not be calibrated", "\n", sep = ""))
        ws["BIOCLIM.O"] <- weights["BIOCLIM.O"] <- 0
    }
    if (ws["BIOCLIM.O"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". original BIOCLIM algorithm (package: BiodiversityR)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(BIOCLIM.O.OLD) == T) {
            if (CATCH.OFF == F) {
                tryCatch(results <- BiodiversityR::ensemble.bioclim.object(x=TrainData.pres, fraction=BIOCLIM.O.fraction, species.name=species.name, factors=factors),
                    error= function(err) {print(paste("original BIOCLIM calibration failed"))},
                    silent=F)
            }else{
                results <- BiodiversityR::ensemble.bioclim.object(x=TrainData.pres, fraction=BIOCLIM.O.fraction, species.name=species.name, factors=factors)
            }
        }else{ 
            results <- BIOCLIM.O.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"BIOCLIM.O"] <- predict.BIOCLIM.O(object=results, newdata=TrainData.vars)
            if (PROBIT == T) {
                if(is.null(BIOCLIM.O.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ BIOCLIM.O"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- BIOCLIM.O.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"BIOCLIM.O.step1"] <- TrainData[,"BIOCLIM.O"]
                TrainData[,"BIOCLIM.O"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "BIOCLIM.O"] < 0), "BIOCLIM.O"] <- 0
                TrainData[which(TrainData[, "BIOCLIM.O"] > 1), "BIOCLIM.O"] <- 1               
            }   
            pred1 <- TrainData[, "BIOCLIM.O"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"BIOCLIM.O"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"BIOCLIM.O"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["BIOCLIM.O"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["BIOCLIM.O"]))
            weights["BIOCLIM.O"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["BIOCLIM.O"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"BIOCLIM.O"] <- predict.BIOCLIM.O(object=results, newdata=TestData.vars)
                if (PROBIT == T) {
                    TestData[,"BIOCLIM.O.step1"] <- TestData[,"BIOCLIM.O"]
                    TestData[,"BIOCLIM.O"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }else{
                    TestData[which(TestData[, "BIOCLIM.O"] < 0), "BIOCLIM.O"] <- 0
                    TestData[which(TestData[, "BIOCLIM.O"] > 1), "BIOCLIM.O"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1,"BIOCLIM.O"]
                TestAbs <- TestData[TestData[,"pb"]==0,"BIOCLIM.O"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["BIOCLIM.O"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["BIOCLIM.O"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: original BIOCLIM evaluation failed","\n\n",sep = ""))
                    ws["BIOCLIM.O"] <- 0
                    weights["BIOCLIM.O"] <- 0
                    TrainData[,"BIOCLIM.O"] <- 0
                    TestData[,"BIOCLIM.O"] <- 0
                }
            }
            if(evaluations.keep ==T) {
                evaluations$BIOCLIM.O.C <- eval1
                evaluations$BIOCLIM.O.T <- eval2
            }
            if (models.keep==T) {
                models$BIOCLIM.O <- results
                models$BIOCLIM.O.PROBIT <- results2
            }
        }else{ 
            cat(paste("\n", "WARNING: original BIOCLIM calibration failed", "\n", "\n"))
            ws["BIOCLIM.O"] <- weights["BIOCLIM.O"] <- 0
            TrainData[,"BIOCLIM.O"] <- 0
            if (no.tests == F) {TestData[,"BIOCLIM.O"] <- 0}
        }
    }
    if (ws["BIOCLIM"] > 0 || ws["DOMAIN"] > 0 || ws["MAHAL"] > 0 || ws["MAHAL01"] > 0) {
        if(is.null(factors) == F) {
            for (i in 1:length(factors)) {
                TrainData.vars <- TrainData.vars[, which(names(TrainData.vars) != factors[i]), drop=F]
                TrainData.pres <- TrainData.pres[, which(names(TrainData.pres) != factors[i]), drop=F]
                if (no.tests == F) {TestData.vars <- TestData.vars[, which(names(TestData.vars) != factors[i]), drop=F]}
            }
        }
    }
    if(ws["BIOCLIM"] > 0 && length(names(TrainData.pres)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "BIOCLIM model can therefore not be calibrated", "\n", sep = ""))
        ws["BIOCLIM"] <- weights["BIOCLIM"] <- 0
    }
    if (ws["BIOCLIM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". BIOCLIM algorithm (package: dismo)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(BIOCLIM.OLD) == T) {
            if (CATCH.OFF == F) {
                tryCatch(results <- dismo::bioclim(x=TrainData.pres),
                    error= function(err) {print(paste("BIOCLIM calibration failed"))},
                    silent=F)
            }else{
                results <- dismo::bioclim(x=TrainData.pres)
            }
        }else{ 
            results <- BIOCLIM.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"BIOCLIM"] <- dismo::predict(object=results, x=TrainData.vars)
            if (PROBIT == T) {
                if(is.null(BIOCLIM.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ BIOCLIM"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- BIOCLIM.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"BIOCLIM.step1"] <- TrainData[,"BIOCLIM"]
                TrainData[,"BIOCLIM"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "BIOCLIM"] < 0), "BIOCLIM"] <- 0
                TrainData[which(TrainData[, "BIOCLIM"] > 1), "BIOCLIM"] <- 1               
            }   
            pred1 <- TrainData[, "BIOCLIM"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"BIOCLIM"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"BIOCLIM"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["BIOCLIM"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["BIOCLIM"]))
            weights["BIOCLIM"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["BIOCLIM"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"BIOCLIM"] <- dismo::predict(object=results, x=TestData.vars)
                if (PROBIT == T) {
                    TestData[,"BIOCLIM.step1"] <- TestData[,"BIOCLIM"]
                    TestData[,"BIOCLIM"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }else{
                    TestData[which(TestData[, "BIOCLIM"] < 0), "BIOCLIM"] <- 0
                    TestData[which(TestData[, "BIOCLIM"] > 1), "BIOCLIM"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1,"BIOCLIM"]
                TestAbs <- TestData[TestData[,"pb"]==0,"BIOCLIM"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["BIOCLIM"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["BIOCLIM"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: BIOCLIM evaluation failed","\n\n",sep = ""))
                    ws["BIOCLIM"] <- 0
                    weights["BIOCLIM"] <- 0
                    TrainData[,"BIOCLIM"] <- 0
                    TestData[,"BIOCLIM"] <- 0
                }
            }
            if(evaluations.keep ==T) {
                evaluations$BIOCLIM.C <- eval1
                evaluations$BIOCLIM.T <- eval2
            }
            if (models.keep==T) {
                models$BIOCLIM <- results
                models$BIOCLIM.PROBIT <- results2
            }
        }else{ 
            cat(paste("\n", "WARNING: BIOCLIM calibration failed", "\n", "\n"))
            ws["BIOCLIM"] <- weights["BIOCLIM"] <- 0
            TrainData[,"BIOCLIM"] <- 0
            if (no.tests == F) {TestData[,"BIOCLIM"] <- 0}
        }
    }
    if (ws["DOMAIN"] > 0) {
        dummy.vars.noDOMAIN <- as.character(NULL)
        if(is.null(dummy.vars) == F) {
            for (i in 1:length(dummy.vars)) {
                max.var <- max(TrainData.pres[, which(names(TrainData.pres) == dummy.vars[i])], na.rm=T)
                min.var <- min(TrainData.pres[, which(names(TrainData.pres) == dummy.vars[i])], na.rm=T)
                if (max.var == min.var) {
                    dummy.vars.noDOMAIN <- c(dummy.vars.noDOMAIN, dummy.vars[i])
                    cat(paste("\n", "WARNING: dummy variable ", dummy.vars[i], " was removed for DOMAIN because it has no variation for the training points", sep=""))
                    TrainData.vars <- TrainData.vars[, which(names(TrainData.vars) != dummy.vars[i]), drop=F]
                    TrainData.pres <- TrainData.pres[, which(names(TrainData.pres) != dummy.vars[i]), drop=F]
                    if (no.tests == F) {TestData.vars <- TestData.vars[, which(names(TestData.vars) != dummy.vars[i]), drop=F]}
                }
            }
            if (length(dummy.vars.noDOMAIN) == 0) {dummy.vars.noDOMAIN <- NULL}
            cat(paste("\n"))
        }
    }
    if(ws["DOMAIN"] > 0 && length(names(TrainData.vars)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "DOMAIN model can therefore not be calibrated", "\n", sep = ""))
        ws["DOMAIN"] <- weights["DOMAIN"] <- 0
    }
    if (ws["DOMAIN"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". DOMAIN algorithm (package: dismo)\n", sep=""))

        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(DOMAIN.OLD) == T) {
            if (CATCH.OFF == F) {
                tryCatch(results <- dismo::domain(x=TrainData.pres),
                    error= function(err) {print(paste("DOMAIN calibration failed"))},
                    silent=F)
            }else{
                results <- dismo::domain(x=TrainData.pres)
            }
        }else{ 
            results <- DOMAIN.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"DOMAIN"] <- dismo::predict(object=results, x=TrainData.vars)
            if (PROBIT == T) {
                if(is.null(DOMAIN.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ DOMAIN"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- DOMAIN.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"DOMAIN.step1"] <- TrainData[,"DOMAIN"]
                TrainData[,"DOMAIN"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "DOMAIN"] < 0), "DOMAIN"] <- 0
                TrainData[which(TrainData[, "DOMAIN"] > 1), "DOMAIN"] <- 1               
            }  
            pred1 <- TrainData[, "DOMAIN"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"DOMAIN"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"DOMAIN"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["DOMAIN"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["DOMAIN"]))
            weights["DOMAIN"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["DOMAIN"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"DOMAIN"] <- dismo::predict(object=results, x=TestData.vars)
                if (PROBIT == T) {
                    TestData[,"DOMAIN.step1"] <- TestData[,"DOMAIN"]
                    TestData[,"DOMAIN"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }else{
                    TestData[which(TestData[, "DOMAIN"] < 0), "DOMAIN"] <- 0
                    TestData[which(TestData[, "DOMAIN"] > 1), "DOMAIN"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1,"DOMAIN"]
                TestAbs <- TestData[TestData[,"pb"]==0,"DOMAIN"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["DOMAIN"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["DOMAIN"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: DOMAIN evaluation failed","\n\n",sep = ""))
                    ws["DOMAIN"] <- 0
                    weights["DOMAIN"] <- 0
                    TrainData[,"DOMAIN"] <- 0
                    TestData[,"DOMAIN"] <- 0
                }
            }
            if(evaluations.keep ==T) {
                evaluations$DOMAIN.C <- eval1
                evaluations$DOMAIN.T <- eval2
            }
            if (models.keep==T) {
                models$DOMAIN <- results
                if (length(dummy.vars.noDOMAIN) > 0) {
                    models$dummy.vars.noDOMAIN <- dummy.vars.noDOMAIN
                    evaluations$dummy.vars.noDOMAIN <- dummy.vars.noDOMAIN
                }
                models$DOMAIN.PROBIT <- results2
            }
        }else{ 
            cat(paste("\n", "WARNING: DOMAIN calibration failed", "\n", "\n"))
            ws["DOMAIN"] <- weights["DOMAIN"] <- 0
            TrainData[,"DOMAIN"] <- 0
            if (no.tests == F) {TestData[,"DOMAIN"] <- 0}
        }
    }
    if (ws["MAHAL"] > 0) {
        if(is.null(dummy.vars) == F) {
            cat(paste("\n", "NOTE: all dummy variables were removed for Mahalanobis algorithm", "\n", sep=""))
            for (i in 1:length(dummy.vars)) {
                TrainData.vars <- TrainData.vars[, which(names(TrainData.vars) != dummy.vars[i]), drop=F]
                TrainData.pres <- TrainData.pres[, which(names(TrainData.pres) != dummy.vars[i]), drop=F]
                if (no.tests == F) {TestData.vars <- TestData.vars[, which(names(TestData.vars) != dummy.vars[i]), drop=F]}
            }
        }
    }
    if(ws["MAHAL"] > 0 && length(names(TrainData.pres)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "MAHAL model can therefore not be calibrated", "\n", sep = ""))
        ws["MAHAL"] <-  weights["MAHAL"] <- 0
    }
    if (ws["MAHAL"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Mahalanobis algorithm (package: dismo)\n", sep=""))

        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(MAHAL.OLD) == T) {
            if (CATCH.OFF == F) {
                tryCatch(results <- dismo::mahal(x=TrainData.pres),
                    error= function(err) {print(paste("Mahalanobis calibration failed"))},
                    silent=F)
            }else{
                results <- dismo::mahal(x=TrainData.pres)
            }
        }else{ 
            results <- MAHAL.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"MAHAL"] <- predict.MAHAL(model=results, newdata=TrainData.vars, PROBIT=PROBIT)
            if (PROBIT == T) {
                if(is.null(MAHAL.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ MAHAL"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- MAHAL.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"MAHAL.step1"] <- TrainData[,"MAHAL"]
                TrainData[,"MAHAL"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "MAHAL"] < 0), "MAHAL"] <- 0
                TrainData[which(TrainData[, "MAHAL"] > 1), "MAHAL"] <- 1               
            }  
            pred1 <- TrainData[, "MAHAL"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"MAHAL"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"MAHAL"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["MAHAL"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["MAHAL"]))
            weights["MAHAL"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["MAHAL"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"MAHAL"] <- predict.MAHAL(model=results, newdata=TestData.vars, PROBIT=PROBIT)
                if (PROBIT == T) {
                    TestData[,"MAHAL.step1"] <- TestData[,"MAHAL"]
                    TestData[,"MAHAL"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }else{
                    TestData[which(TestData[, "MAHAL"] < 0), "MAHAL"] <- 0
                    TestData[which(TestData[, "MAHAL"] > 1), "MAHAL"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1,"MAHAL"]
                TestAbs <- TestData[TestData[,"pb"]==0,"MAHAL"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["MAHAL"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["MAHAL"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: Mahalanobis evaluation failed","\n\n",sep = ""))
                    ws["MAHAL"] <- 0
                    weights["MAHAL"] <- 0
                    TrainData[,"MAHAL"] <- 0
                    TestData[,"MAHAL"] <- 0
                }
            }
            if(evaluations.keep ==T) {
                evaluations$MAHAL.C <- eval1
                evaluations$MAHAL.T <- eval2
            }
            if (models.keep==T) {
                models$MAHAL <- results
                models$MAHAL.PROBIT <- results2
            }
        }else{ 
            cat(paste("\n", "WARNING: Mahalanobis calibration failed", "\n", "\n"))
            ws["MAHAL"] <- weights["MAHAL"] <- 0
            TrainData[,"MAHAL"] <- 0
            if (no.tests == F) {TestData[,"MAHAL"] <- 0}
        }
    }
    if (ws["MAHAL01"] > 0) {
        if(is.null(dummy.vars) == F) {
            cat(paste("\n", "NOTE: all dummy variables were removed for Mahalanobis algorithm", "\n", sep=""))
            for (i in 1:length(dummy.vars)) {
                TrainData.vars <- TrainData.vars[, which(names(TrainData.vars) != dummy.vars[i]), drop=F]
                TrainData.pres <- TrainData.pres[, which(names(TrainData.pres) != dummy.vars[i]), drop=F]
                if (no.tests == F) {TestData.vars <- TestData.vars[, which(names(TestData.vars) != dummy.vars[i]), drop=F]}
            }
        }
    }
    if(ws["MAHAL01"] > 0 && length(names(TrainData.pres)) == 0) {
        cat(paste("\n", "WARNING: no explanatory variables available", sep = ""))
        cat(paste("\n", "MAHAL01 model can therefore not be calibrated", "\n", sep = ""))
        ws["MAHAL01"] <- weights["MAHAL01"] <- 0
    }
    if (ws["MAHAL01"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Mahalanobis algorithm (transformed within 0 to 1 interval)", "\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(MAHAL01.OLD) == T) {
            if (CATCH.OFF == F) {
                tryCatch(results <- dismo::mahal(x=TrainData.pres),
                    error= function(err) {print(paste("transformed Mahalanobis calibration failed"))},
                    silent=F)
            }else{
                results <- dismo::mahal(x=TrainData.pres)
            }
        }else{ 
            results <- MAHAL01.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"MAHAL01"] <- predict.MAHAL01(model=results, newdata=TrainData.vars, MAHAL.shape=MAHAL.shape)
            if (PROBIT == T) {
                if(is.null(MAHAL01.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ MAHAL01"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights1, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- MAHAL01.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"MAHAL01.step1"] <- TrainData[,"MAHAL01"]
                TrainData[,"MAHAL01"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }else{
                TrainData[which(TrainData[, "MAHAL01"] < 0), "MAHAL01"] <- 0
                TrainData[which(TrainData[, "MAHAL01"] > 1), "MAHAL01"] <- 1               
            }   
            pred1 <- TrainData[, "MAHAL01"]
            pred1[pred1 < 0.0000000001] <- 0.0000000001
            pred1[pred1 > 0.9999999999] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1, "MAHAL01"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0, "MAHAL01"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["MAHAL01"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["MAHAL01"]))
            weights["MAHAL01"] <- max(c(eval1@auc, 0), na.rm=T)
            AUC.calibration["MAHAL01"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"MAHAL01"] <- predict.MAHAL01(model=results, newdata=TestData.vars, MAHAL.shape=MAHAL.shape)
                if (PROBIT == T) {
                    TestData[,"MAHAL01.step1"] <- TestData[,"MAHAL01"]
                    TestData[,"MAHAL01"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }else{
                    TestData[which(TestData[, "MAHAL01"] < 0), "MAHAL01"] <- 0
                    TestData[which(TestData[, "MAHAL01"] > 1), "MAHAL01"] <- 1   
                }
                TestPres <- TestData[TestData[,"pb"]==1, "MAHAL01"]
                TestAbs <- TestData[TestData[,"pb"]==0, "MAHAL01"]
                eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["MAHAL01"] <- max(c(eval2@auc, 0), na.rm=T)
                    AUC.testing["MAHAL01"] <- max(c(eval2@auc, 0), na.rm=T)
                }else{
                    cat(paste("\n", "WARNING: transformed Mahalanobis evaluation failed","\n\n",sep = ""))
                    ws["MAHAL01"] <- 0
                    weights["MAHAL01"] <- 0
                    TrainData[,"MAHAL01"] <- 0
                    TestData[,"MAHAL01"] <- 0
                }
            }
            if(evaluations.keep ==T) {
                evaluations$MAHAL01.C <- eval1
                evaluations$MAHAL01.T <- eval2
            }
            if (models.keep==T) {
                models$MAHAL01 <- results
                models$formulae$MAHAL.shape <- MAHAL.shape
                models$MAHAL01.PROBIT <- results2
            }
        }else{ 
            cat(paste("\n", "WARNING: transformed Mahalanobis calibration failed", "\n", "\n"))
            ws["MAHAL01"] <- weights["MAHAL01"]  <- 0
            TrainData[,"MAHAL01"] <- 0
            if (no.tests == F) {TestData[,"MAHAL01"] <- 0}
        }
    }
#
    models$thresholds <- thresholds
#
    if(ENSEMBLE.tune == F) {
        if (sum(ws, na.rm=T) > 0) {
            cat(paste("\n", "Ensemble weights based directly on input weights scaled to sum up to 1", "\n", sep = ""))
            print(ws)
        }else{
            cat(paste("\n", "NOTE: no positive input weights", "\n", sep = ""))
        }
        if(evaluations.keep == T) {evaluations$ensemble.weights <- ws}
        if(models.keep==T) {models$output.weights <- ws}
    }else{
# use different strategies for calculating the ensemble model
# similar to using test data for calculating input AUC, use internal test data for calculating best ensemble
# recalculating AUC does not require much computing time - initial calculations kept to spot problems for specific algorithms
        cat(paste("\n", "Weights tuned by ensemble.strategy function", sep=""))
        strategy.results <- ensemble.strategy(TrainData=TrainData, TestData=TestData,
            ENSEMBLE.exponent=ENSEMBLE.exponent, ENSEMBLE.best=ENSEMBLE.best, ENSEMBLE.min=ENSEMBLE.min)
        ws <- strategy.results$weights
        if (sum(ws, na.rm=T) > 0) {
            cat(paste("\n", "Minimum input weight is ", ENSEMBLE.weight.min, "\n", sep=""))
            ws2 <- ws
            while(min(ws2) < ENSEMBLE.weight.min) {
                ws2 <- ws2[-which.min(ws2)]
                ws2 <- ensemble.weights(weights=ws2, exponent=1, best=0, min.weight=0)
            }
            ws[] <- 0
            for (i in 1:length(ws2)) {ws[which(names(ws) == names(ws2)[i])] <- ws2[i]}
            cat(paste("\n", "Weights for ensemble forecasting", "\n", sep = ""))
            print(ws)
        }else{
            ENSEMBLE.tune <- FALSE
        }
        if(evaluations.keep == T) {evaluations$STRATEGY.weights <- ws}
        if(models.keep==T) {models$output.weights <- ws}
    }
# do not return ensemble tests for no models when tuning is not implemented 
    if((sum(ws > 0, na.rm=T) > 0) || (ENSEMBLE.tune == T)) {
        TrainData[,"ENSEMBLE"] <- ws["MAXENT"]*TrainData[,"MAXENT"] + ws["MAXLIKE"]*TrainData[,"MAXLIKE"] + ws["GBM"]*TrainData[,"GBM"] +
            ws["GBMSTEP"]*TrainData[,"GBMSTEP"] + ws["RF"]*TrainData[,"RF"] + ws["GLM"]*TrainData[,"GLM"] +
            ws["GLMSTEP"]*TrainData[,"GLMSTEP"] + ws["GAM"]*TrainData[,"GAM"] + ws["GAMSTEP"]*TrainData[,"GAMSTEP"] +
            ws["MGCV"]*TrainData[,"MGCV"] + ws["MGCVFIX"]*TrainData[,"MGCVFIX"] + ws["EARTH"]*TrainData[,"EARTH"] +
            ws["RPART"]*TrainData[,"RPART"] + ws["NNET"]*TrainData[,"NNET"] + ws["FDA"]*TrainData[,"FDA"] +
            ws["SVM"]*TrainData[,"SVM"] + ws["SVME"]*TrainData[,"SVME"] + ws["GLMNET"]*TrainData[,"GLMNET"] +
            ws["BIOCLIM.O"]*TrainData[,"BIOCLIM.O"] + ws["BIOCLIM"]*TrainData[,"BIOCLIM"] +
            ws["DOMAIN"]*TrainData[,"DOMAIN"] + ws["MAHAL"]*TrainData[,"MAHAL"] + ws["MAHAL01"]*TrainData[,"MAHAL01"]
        pred1 <- TrainData[, "ENSEMBLE"]
        pred1[pred1 < 0.0000000001] <- 0.0000000001
        pred1[pred1 > 0.9999999999] <- 0.9999999999
        pred2 <- rep(mean(obs1), times=length(pred1))
        if (no.tests == F) {
            TestData[,"ENSEMBLE"] <- ws["MAXENT"]*TestData[,"MAXENT"] + ws["MAXLIKE"]*TestData[,"MAXLIKE"] + ws["GBM"]*TestData[,"GBM"] +
                ws["GBMSTEP"]*TestData[,"GBMSTEP"] + ws["RF"]*TestData[,"RF"] + ws["GLM"]*TestData[,"GLM"] +
                ws["GLMSTEP"]*TestData[,"GLMSTEP"] + ws["GAM"]*TestData[,"GAM"] + ws["GAMSTEP"]*TestData[,"GAMSTEP"] +
                ws["MGCV"]*TestData[,"MGCV"] + ws["MGCVFIX"]*TestData[,"MGCVFIX"] + ws["EARTH"]*TestData[,"EARTH"] +
                ws["RPART"]*TestData[,"RPART"] + ws["NNET"]*TestData[,"NNET"] + ws["FDA"]*TestData[,"FDA"] +
                ws["SVM"]*TestData[,"SVM"] + ws["SVME"]*TestData[,"SVME"] + ws["GLMNET"]*TestData[,"GLMNET"] + 
                ws["BIOCLIM.O"]*TestData[,"BIOCLIM.O"] + ws["BIOCLIM"]*TestData[,"BIOCLIM"] +
                ws["DOMAIN"]*TestData[,"DOMAIN"] + ws["MAHAL"]*TestData[,"MAHAL"] + ws["MAHAL01"]*TestData[,"MAHAL01"]
        }
        mc <- mc+1
        cat(paste("\n\n", mc, ". Ensemble algorithm\n", sep=""))
        eval1 <- eval2 <- NULL
        cat(paste("\n", "Ensemble evaluation with calibration data", "\n\n", sep = ""))
        cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n", sep = ""))
        cat(paste("Residual deviance if all predictions were ", mean(obs1), " (prevalence): ", dismo::calc.deviance(obs=obs1, pred=pred2, calc.mean=F), "\n", sep = ""))
# worst possible predictions
        numpres <- sum(TrainData[, "pb"])
        numabs <- nrow(TrainData) - numpres
        pred1 <- rep(0.000000001, numpres)
        pred2 <- rep(0.999999999, numabs)
        pred3 <- c(pred1, pred2)
        null.dev.cal <- dismo::calc.deviance(obs=obs1, pred=pred3, calc.mean=F)
        cat(paste("Residual deviance for worst possible predictions: ", null.dev.cal, "\n\n", sep = ""))

        TrainPres <- as.numeric(TrainData[TrainData[,"pb"]==1, "ENSEMBLE"])
        TrainAbs <- as.numeric(TrainData[TrainData[,"pb"]==0, "ENSEMBLE"])
        if (sum(TrainPres, na.rm=T) <= 0 || sum(TrainAbs, na.rm=T) <= 0) {
            cat(paste("\n", "NOTE: not possible to evaluate the ensemble model since calibration probabilities not available", "\n", sep = ""))
        }else{
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["ENSEMBLE"] <- ensemble.threshold(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            AUC.calibration["ENSEMBLE"] <- max(c(eval1@auc, 0), na.rm=T)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["ENSEMBLE"]))
            if (models.keep == T) {models$thresholds <- thresholds}
            if (no.tests == F) {
                cat(paste("\n", "Ensemble evaluation with testing data", "\n\n", sep = ""))
                TestPres <- as.numeric(TestData[TestData[,"pb"]==1,"ENSEMBLE"])
                TestAbs <- as.numeric(TestData[TestData[,"pb"]==0,"ENSEMBLE"])
                if (sum(TestPres, na.rm=T) <= 0 || sum(TestAbs, na.rm=T) <= 0) {
                    cat(paste("\n", "NOTE: not possible to evaluate the ensemble model since evaluation probabilities not available", "\n", sep = ""))
                }else{
                    eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
                    AUC.testing["ENSEMBLE"] <- max(c(eval2@auc, 0), na.rm=T)
                }
            }
            if(evaluations.keep==T) {
                evaluations$ENSEMBLE.C <- eval1
                evaluations$ENSEMBLE.T <- eval2
            }
        }
    }
    if(models.keep==T) {
        models$TrainData <- TrainData
        if (no.tests == F) {models$TestData <- TestData}
        models$var.names <- var.names
        models$p <- p
        models$pt <- pt
        models$a <- a
        models$at <- at
        models$MAXENT.a <- MAXENT.a
        models$var.names <- var.names
        models$factors <- factors
        models$dummy.vars <- dummy.vars
        models$dummy.vars.noDOMAIN <- dummy.vars.noDOMAIN
    }
    if(evaluations.keep == T) {
        evaluations$TrainData <- TrainData
        if (no.tests == F) {evaluations$TestData <- TestData}
        evaluations$p <- p
        evaluations$pt <- pt
        evaluations$a <- a
        evaluations$at <- at
        evaluations$MAXENT.a <- MAXENT.a
        evaluations$var.names <- var.names
        evaluations$factors <- factors
        evaluations$dummy.vars <- dummy.vars
        evaluations$dummy.vars.noDOMAIN <- dummy.vars.noDOMAIN
    }
    remove(Yweights1, envir=.BiodiversityR)
    remove(TrainData, envir=.BiodiversityR)
    remove(TrainData.vars, envir=.BiodiversityR)
    remove(TrainData.numvars, envir=.BiodiversityR)
    remove(TrainData.pres, envir=.BiodiversityR)
    if (no.tests == F) {
        remove(TestData, envir=.BiodiversityR)
        remove(TestData.vars, envir=.BiodiversityR)
        remove(TestData.numvars, envir=.BiodiversityR)
    }
    if (models.save==T && models.keep==T) {
        ensemble.models <- models
        save(ensemble.models, file=paste(getwd(), "/models/", models$species.name, "_models", sep=""), compress="xz")
    }
    if (models.keep == F) {models <- NULL}     
    result <- list(evaluations=evaluations, AUC.calibration=AUC.calibration, AUC.testing=AUC.testing, models=models, VIF=newVIF, call=match.call() )
    cat(paste("\n\n"))
    if (SINK==T && OLD.SINK==F) {sink(file=NULL, append=T)}
    cat(paste("\n", "finished model calibrations (function ensemble.calibrate.models)", "\n\n", sep = ""))
    return(result)
}

