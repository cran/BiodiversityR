`ensemble.test.splits` <- function(
    x=NULL, p=NULL, a=NULL, an=1000, excludep=FALSE, ext=NULL, k=5, 
    TrainData=NULL,
    layer.drops=NULL, VIF=FALSE,
    PLOTS=FALSE, data.keep=FALSE,
    ENSEMBLE.decay=1, ENSEMBLE.best=0, ENSEMBLE.min=0.7,
    input.weights=NULL,
    MAXENT=1, GBM=1, GBMSTEP=1, RF=1, GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=1, MGCV=1, MGCVFIX=0, 
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, BIOCLIM=1, DOMAIN=1, MAHAL=1, 
    GEODIST=0, 
    Yweights="BIOMOD", factors=NULL, dummy.vars=NULL,
    formulae.defaults=TRUE, maxit=100,
    MAXENT.a=NULL, MAXENT.an=10000, MAXENT.BackData=NULL, MAXENT.path=paste(getwd(), "/models/maxent", sep=""),
    GBM.formula=NULL, GBM.n.trees=2001, 
    GBMSTEP.gbm.x=2:(ncol(TrainData.orig)), GBMSTEP.tree.complexity=5, GBMSTEP.learning.rate=0.005, 
    GBMSTEP.bag.fraction=0.5, GBMSTEP.step.size=100, 
    RF.formula=NULL, RF.ntree=751, RF.mtry=floor(sqrt(ncol(TrainData.orig)-1)), 
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
    MAHAL.shape=1, 
    GEODIST.file.name="Species001", RASTER.format="raster"
)
{
    .BiodiversityR <- new.env()
    k <- as.integer(k)
    if (k < 2) {
        cat(paste("\n", "NOTE: parameter k was set to be smaller than 2", sep = ""))
        cat(paste("\n", "default value of 5 therefore set for parameter k", "\n", sep = ""))
        k <- 5
    }
    if (! require(dismo)) {stop("Please install the dismo package")}
# check data
    if (is.null(TrainData) == T) {
        if(is.null(x) == T) {stop("value for parameter x is missing (RasterStack object)")}
        if(inherits(x,"RasterStack") == F) {stop("x is not a RasterStack object")}
        if(is.null(p) == T) {stop("presence locations are missing (parameter p)")}
    }
# geoDist requires presence locations
    if ((GEODIST > 0) && (is.null(p) == T)) {stop("presence locations are missing for geoDist")}
# check TrainData
    if (is.null(TrainData) == F) {
        TrainData <- data.frame(TrainData)
        if (colnames(TrainData)[1] !="pb") {stop("first column for TrainData should be 'pb' containing presence (1) and absence (0) data")}
        if ((is.null(x) == F) && (nlayers(x) != (ncol(TrainData)-1))) {
            cat(paste("\n", "WARNING: different number of explanatory variables in rasterStack and TrainData", sep = ""))
        }
    }
# modify rasterStack x only if TrainData is not used
    if (is.null(TrainData) == T) {
        if (is.null(layer.drops) == F) {
            vars <- names(x)
            layer.drops <- as.character(layer.drops)
            factors <- as.character(factors)
            dummy.vars <- as.character(dummy.vars)
            nd <- length(layer.drops)
            for (i in 1:nd) {     
                if (any(vars==layer.drops[i])==FALSE) {
                    cat(paste("\n", "WARNING: variable to exclude '", layer.drops[i], "' not among grid layers", "\n", sep = ""))
                }else{
                    cat(paste("\n", "NOTE: variable ", layer.drops[i], " will not be included as explanatory variable", "\n", sep = ""))
                    x <- dropLayer(x, which(names(x) == layer.drops[i]))
                    vars <- names(x)
                    if (is.null(factors) == F) {
                        factors <- factors[factors != layer.drops[i]]
                        if(length(factors) == 0) {factors <- NULL}
                    }
                    if (is.null(dummy.vars) == F) {
                        dummy.vars <- dummy.vars[dummy.vars != layer.drops[i]]
                        if(length(dummy.vars) == 0) {dummy.vars <- NULL}
                    }
                }
            }
        }
        if (is.null(factors) == F) {
            vars <- names(x)
            factors <- as.character(factors)
            nf <- length(factors)
            for (i in 1:nf) {
                if (any(vars==factors[i])==FALSE) {
                    cat(paste("\n", "WARNING: categorical variable '", factors[i], "' not among grid layers", "\n", sep = ""))
                }
            }
        }
        if (is.null(dummy.vars) == F) {
            vars <- names(x)
            dummy.vars <- as.character(dummy.vars)
            nf <- length(dummy.vars)
            for (i in 1:nf) {
                if (any(vars==dummy.vars[i])==FALSE) {
                    cat(paste("\n", "WARNING: dummy variable '", dummy.vars[i], "' not among grid layers", "\n", sep = ""))
                }
            }
        }
        # set minimum and maximum values
            for (i in 1:nlayers(x)) {
                x[[i]] <- setMinMax(x[[i]])
            }
        # declare factor layers
        if(is.null(factors)==F) {
            for (i in 1:length(factors)) {
                j <- which(names(x) == factors[i])
                x[[j]] <- raster::as.factor(x[[j]])
            }
        }
# modify TrainData
    }else{
        if (is.null(layer.drops) == F) {
            vars <- colnames(TrainData)
            layer.drops <- as.character(layer.drops)
            factors <- as.character(factors)
            dummy.vars <- as.character(dummy.vars)
            nd <- length(layer.drops)
            for (i in 1:nd) {     
                if (any(vars==layer.drops[i])==FALSE) {
                    cat(paste("\n", "WARNING: variable to exclude '", layer.drops[i], "' not among columns of TrainData", "\n", sep = ""))
                }else{
                    cat(paste("\n", "NOTE: variable ", layer.drops[i], " will not be included as explanatory variable", "\n", sep = ""))
                    TrainData <- TrainData[, which(colnames(TrainData) != layer.drops[i])]
                    vars <- colnames(TrainData)
                    if (is.null(factors) == F) {
                        factors <- factors[factors != layer.drops[i]]
                        if(length(factors) == 0) {factors <- NULL}
                    }
                    if (is.null(dummy.vars) == F) {
                        dummy.vars <- dummy.vars[dummy.vars != layer.drops[i]]
                        if(length(dummy.vars) == 0) {dummy.vars <- NULL}
                    }
                }
            }
        }
        if (is.null(factors) == F) {
            vars <- colnames(TrainData)
            factors <- as.character(factors)
            nf <- length(factors)
            for (i in 1:nf) {
                if (any(vars==factors[i])==FALSE) {
                     cat(paste("\n", "WARNING: categorical variable '", factors[i], "' not among columns of TrainData", "\n", sep = ""))
                }
            }
        }
        if (is.null(dummy.vars) == F) {
            vars <- colnames(TrainData)
            dummy.vars <- as.character(dummy.vars)
            nf <- length(dummy.vars)
            for (i in 1:nf) {
                if (any(vars==dummy.vars[i])==FALSE) {
                    cat(paste("\n", "WARNING: dummy variable '", dummy.vars[i], "' not among columns of TrainData", "\n", sep = ""))
                }
            }
        }
    }
    if (is.null(input.weights)==F) {
        MAXENT <- max(c(input.weights["MAXENT"], -1), na.rm=T)
        GBM <- max(c(input.weights["GBM"], -1), na.rm=T)
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
        GEODIST <- max(c(input.weights["GEODIST"], -1), na.rm=T)
    }
    if (MAXENT > 0) {
	jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
	if (!file.exists(jar)) {stop('maxent program is missing: ', jar, '\nPlease download it here: http://www.cs.princeton.edu/~schapire/maxent/')}
    }
    if (formulae.defaults == T) {
        if (is.null(TrainData) == T) {
            formulae <- ensemble.formulae(x, factors=factors, dummy.vars=dummy.vars)
        }else{
            formulae <- ensemble.formulae(TrainData, factors=factors, dummy.vars=dummy.vars)
        }
    }
    if (GBM > 0) {
        if (! require(gbm)) {stop("Please install the gbm package")}
        if (is.null(GBM.formula) == T && formulae.defaults == T) {GBM.formula <- formulae$GBM.formula}
        if (is.null(GBM.formula) == T) {stop("Please provide the GBM.formula (hint: use ensemble.formulae function)")}
        environment(GBM.formula) <- .BiodiversityR
    }
    if (GBMSTEP > 0) {
        if (! require(gbm)) {stop("Please install the gbm package")}
    }
    if (RF > 0) {
        if (! require(randomForest)) {stop("Please install the randomForest package")}
        if (is.null(RF.formula) == T && formulae.defaults == T) {RF.formula <- formulae$RF.formula}
        if (is.null(RF.formula) == T) {stop("Please provide the RF.formula (hint: use ensemble.formulae function)")}
        environment(RF.formula) <- .BiodiversityR
        if (identical(RF.ntree, trunc(RF.ntree/2)) == F) {RF.ntree <- RF.ntree + 1}
    }
    if (GLM > 0) {
        if (is.null(GLM.formula) == T && formulae.defaults == T) {GLM.formula <- formulae$GLM.formula}
        if (is.null(GLM.formula) == T) {stop("Please provide the GLM.formula (hint: use ensemble.formulae function)")}
        environment(GLM.formula) <- .BiodiversityR
        assign("GLM.family", GLM.family, envir=.BiodiversityR)
    }
    if (GLMSTEP > 0) {
        if (! require(MASS)) {stop("Please install the MASS package")}
        if (is.null(STEP.formula) == T && formulae.defaults == T) {STEP.formula <- formulae$STEP.formula}
        if (is.null(GLMSTEP.scope) == T && formulae.defaults == T) {GLMSTEP.scope <- formulae$GLMSTEP.scope}
        if (is.null(STEP.formula) == T) {stop("Please provide the STEP.formula (hint: use ensemble.formulae function)")}
        if (is.null(GLMSTEP.scope) == T) {stop("Please provide the GLMSTEP.scope (hint: use ensemble.formulae function)")}
        environment(STEP.formula) <- .BiodiversityR
        assign("GLM.family", GLM.family, envir=.BiodiversityR)
    }
    if (GAM > 0) {
        if (! require(gam)) {stop("Please install the gam package")}
        if (is.null(GAM.formula) == T && formulae.defaults == T) {GAM.formula <- formulae$GAM.formula}
        if (is.null(GAM.formula) == T) {stop("Please provide the GAM.formula (hint: use ensemble.formulae function)")}
        environment(GAM.formula) <- .BiodiversityR
        assign("GAM.family", GAM.family, envir=.BiodiversityR)
        detach(package:gam)   
    }
    if (GAMSTEP > 0) {
        if (! require(gam)) {stop("Please install the gam package")}
        if (is.null(STEP.formula) == T && formulae.defaults == T) {STEP.formula <- formulae$STEP.formula}
        if (is.null(GAMSTEP.scope) == T && formulae.defaults == T) {GAMSTEP.scope <- formulae$GAMSTEP.scope}
        if (is.null(STEP.formula) == T) {stop("Please provide the STEP.formula (hint: use ensemble.formulae function)")}
        if (is.null(GAMSTEP.scope) == T) {stop("Please provide the GAMSTEP.scope (hint: use ensemble.formulae function)")}
        environment(STEP.formula) <- .BiodiversityR
        assign("GAM.family", GAM.family, envir=.BiodiversityR)
        detach(package:gam)   
    }
    if (MGCV > 0) {
        if (! require(mgcv)) {stop("Please install the mgcv package")}
        if (is.null(MGCV.formula) == T && formulae.defaults == T) {MGCV.formula <- formulae$MGCV.formula}
        if (is.null(MGCV.formula) == T) {stop("Please provide the MGCV.formula (hint: use ensemble.formulae function)")}
        environment(MGCV.formula) <- .BiodiversityR
        assign("GAM.family", GAM.family, envir=.BiodiversityR)
        detach(package:mgcv)
#         get the probabilities from MGCV
            predict.mgcv <- function(object, newdata, type="response") {
                p <- predict(object=object, newdata=newdata, type=type)
                return(as.numeric(p))
            } 
    }
    if (MGCVFIX > 0) {
        if (! require(mgcv)) {stop("Please install the mgcv package")}
        if (is.null(MGCVFIX.formula) == T && formulae.defaults == T) {MGCVFIX.formula <- formulae$MGCVFIX.formula}
        if (is.null(MGCVFIX.formula) == T) {stop("Please provide the MGCVFIX.formula (hint: use ensemble.formulae function)")}
        environment(MGCVFIX.formula) <- .BiodiversityR
        assign("GAM.family", GAM.family, envir=.BiodiversityR)
        detach(package:mgcv)
#         get the probabilities from MGCV
            predict.mgcv <- function(object, newdata, type="response") {
                p <- predict(object=object, newdata=newdata, type=type)
                return(as.numeric(p))
            }
    }
    if (EARTH > 0) {
        if (! require(earth)) {stop("Please install the earth package")}
        if (is.null(EARTH.formula) == T && formulae.defaults == T) {EARTH.formula <- formulae$EARTH.formula}
        if (is.null(EARTH.formula) == T) {stop("Please provide the EARTH.formula (hint: use ensemble.formulae function)")}
        environment(EARTH.formula) <- .BiodiversityR
#         get the probabilities from earth
            predict.earth2 <- function(object, newdata, type="response") {
                p <- predict(object=object, newdata=newdata, type=type)
                return(as.numeric(p))
            }
    }
    if (RPART > 0) {
        if (! require(rpart)) {stop("Please install the rpart package")}
        if (is.null(RPART.formula) == T && formulae.defaults == T) {RPART.formula <- formulae$RPART.formula}
        if (is.null(RPART.formula) == T) {stop("Please provide the RPART.formula (hint: use ensemble.formulae function)")}
        environment(RPART.formula) <- .BiodiversityR
    }
    if (NNET > 0) {
        if (! require(nnet)) {stop("Please install the nnet package")}
        if (is.null(NNET.formula) == T && formulae.defaults == T) {NNET.formula <- formulae$NNET.formula}
        if (is.null(NNET.formula) == T) {stop("Please provide the NNET.formula (hint: use ensemble.formulae function)")}
        environment(NNET.formula) <- .BiodiversityR
#         get the probabilities from nnet
            predict.nnet2 <- function(object, newdata, type="raw") {
                p <- predict(object=object, newdata=newdata, type=type)
                return(as.numeric(p))
            }
    }
    if (FDA > 0) {
        if (! require(mda)) {stop("Please install the mda package")}
        if (is.null(FDA.formula) == T && formulae.defaults == T) {FDA.formula <- formulae$FDA.formula}
        if (is.null(FDA.formula) == T) {stop("Please provide the FDA.formula (hint: use ensemble.formulae function)")}
        environment(FDA.formula) <- .BiodiversityR
    }
    if (SVM > 0) {
        if (! require(kernlab)) {stop("Please install the kernlab package")}
        if (is.null(SVM.formula) == T && formulae.defaults == T) {SVM.formula <- formulae$SVM.formula}
        if (is.null(SVM.formula) == T) {stop("Please provide the SVM.formula (hint: use ensemble.formulae function)")}
        environment(SVM.formula) <- .BiodiversityR
    }
    if (SVME > 0) {
        if (! require(e1071)) {stop("Please install the e1071 package")}
        if (is.null(SVME.formula) == T && formulae.defaults == T) {SVME.formula <- formulae$SVME.formula}
        if (is.null(SVME.formula) == T) {stop("Please provide the SVME.formula (hint: use ensemble.formulae function)")}
        environment(SVME.formula) <- .BiodiversityR
#         get the probabilities from svm
            predict.svme <- function(model, newdata) {
                p <- predict(model, newdata, probability=T)
                return(attr(p, "probabilities")[,1])
             }
    }
    if (MAHAL > 0) {
#         get the probabilities from mahal
            predict.mahal <- function(model, newdata, MAHAL.shape) {
                p <- dismo::predict(object=model, x=newdata)
                p <- p - 1 - MAHAL.shape
                p <- abs(p)
                p <- MAHAL.shape / p
                return(p)
             }
    }
# create TrainData and TestData
    if (is.null(TrainData) == F) {
        if(any(is.na(TrainData))) {
            cat(paste("\n", "WARNING: sample units with missing data removed from calibration data","\n\n",sep = ""))
        }
        TrainValid <- complete.cases(TrainData)
        TrainData <- TrainData[TrainValid,]
    }else{
        if (is.null(a)==T) {            
            if (excludep == T) {
                a <- randomPoints(x[[1]], n=an, p=p, ext=ext, excludep=T)
            }else{
                a <- randomPoints(x[[1]], n=an, p=NULL, ext=ext, excludep=F)
            }  
        }
        TrainData <- prepareData(x, p, b=a, factors=factors, xy=FALSE)
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
        TrainData <- prepareData(x, p, b=a, factors=factors, xy=FALSE)
    }
    TrainData.orig <- TrainData
    assign("TrainData.orig", TrainData.orig, envir=.BiodiversityR)
# create background data for MAXENT
# these same background data is used in each cross-validation run
    if (MAXENT > 0) {
        if (is.null(MAXENT.BackData) == T) {
            if (is.null(x) == T) {
                cat(paste("\n", "WARNING: not possible to create MAXENT.BackData as RasterStack x is missing", sep = "")) 
                cat(paste("\n", "MAXENT model will not be calibrated", "\n", sep = "")) 
                MAXENT <- 0
            }else{
# default option of MAXENT is to exclude presence locations
                if (is.null(MAXENT.a) == T) {
                    MAXENT.a <- randomPoints(x[[1]], n=MAXENT.an, p=p, ext=ext, excludep=T)
                }
                MAXENT.BackData <- extract(x, y=MAXENT.a)
            }
        }
        MAXENT.BackData <- data.frame(MAXENT.BackData)
        TestValid <- complete.cases(MAXENT.BackData)
        MAXENT.BackData <- MAXENT.BackData[TestValid,]
    }
#
    output.rownames <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL", "GEODIST", "ENSEMBLE")
    ENSEMBLE.tune <- TRUE
    if((length(ENSEMBLE.decay) == 1) && (length(ENSEMBLE.best) == 1) && (length(ENSEMBLE.min) == 1)) {ENSEMBLE.tune <- FALSE}
    if(ENSEMBLE.tune == F) {
        output <- array(0, dim=c(length(output.rownames), k+1))
        rownames(output) <- output.rownames
        colnames(output) <- c(paste("T_", c(1:k), sep=""),"MEAN")
    }else{
        output <- array(0, dim=c(length(output.rownames), 2*k+2))
        rownames(output) <- output.rownames
        colnames(output) <- c(paste("T_", c(1:k), sep=""), "MEAN.T", paste("S_", c(1:k), sep=""), "MEAN")
    }
# keep data for final checks with suggested weights
    TestData.all <- vector("list", k)
# start cross-validations
    groupp <- kfold(TrainData, k=k, by=TrainData[,"pb"])
    for (i in 1:k){
        cat(paste("\n", "K-FOLD CROSS-VALIDATION RUN: ", i, "\n\n", sep = ""))
        TrainData.c <- TrainData[groupp != i,]
        TestData.c <- TrainData[groupp == i,]
# avoid problems with different factor levels
        if(is.null(factors)==F) {
            for (i in 1:length(factors)) {
                if (identical(levels(TrainData.c[,factors[i]]),levels(TestData.c[,factors[i]]))==F) {
                    cat(paste("\n", "WARNING: differences in factor levels between calibration and evaluation data (variable ", factors[i], ")", "\n", sep = ""))
                    cat(paste("same levels set for both data sets to avoid problems with evaluations for RF and SVM","\n", sep = ""))
                    uniquelevels <- unique(c(levels(TrainData.c[,factors[i]]), levels(TestData.c[,factors[i]])))
                    levels(TrainData.c[,factors[i]]) <- uniquelevels
                    levels(TestData.c[,factors[i]]) <- uniquelevels
                } 
            }
        }
        assign("TrainData.c", TrainData.c, envir=.BiodiversityR)
        assign("TestData.c", TestData.c, envir=.BiodiversityR)
        tests <- ensemble.test(x=x, 
            TrainData=TrainData.c, TestData=TestData.c,
            VIF=VIF,
            PLOTS=PLOTS, evaluations.keep=T, models.keep=F,
            ENSEMBLE.decay=ENSEMBLE.decay, ENSEMBLE.best=ENSEMBLE.best, ENSEMBLE.min=ENSEMBLE.min, 
            MAXENT=MAXENT, GBM=GBM, GBMSTEP=GBMSTEP, RF=RF, GLM=GLM, GLMSTEP=GLMSTEP, 
            GAM=GAM, GAMSTEP=GAMSTEP, MGCV=MGCV, MGCVFIX=MGCVFIX, EARTH=EARTH, RPART=RPART, 
            NNET=NNET, FDA=FDA, SVM=SVM, SVME=SVME, BIOCLIM=BIOCLIM, DOMAIN=DOMAIN, MAHAL=MAHAL,
            GEODIST=GEODIST,   
            Yweights=Yweights, factors=factors, dummy.vars=dummy.vars,
            maxit=maxit,
            MAXENT.a=MAXENT.a, MAXENT.an=MAXENT.an, MAXENT.BackData=MAXENT.BackData, MAXENT.path=MAXENT.path,
            GBM.formula=GBM.formula, GBM.n.trees=GBM.n.trees, 
            GBMSTEP.gbm.x=GBMSTEP.gbm.x, GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
            GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction,
            GBMSTEP.step.size=GBMSTEP.step.size, 
            RF.formula=RF.formula, RF.ntree=RF.ntree, RF.mtry=RF.mtry, 
            GLM.formula=GLM.formula, GLM.family=GLM.family, 
            GLMSTEP.k=GLMSTEP.k, GLMSTEP.steps=GLMSTEP.steps, STEP.formula=STEP.formula, 
            GLMSTEP.scope=GLMSTEP.scope, 
            GAM.formula=GAM.formula, GAM.family=GAM.family, 
            GAMSTEP.steps=GAMSTEP.steps, GAMSTEP.scope=GAMSTEP.scope, GAMSTEP.pos=GAMSTEP.pos, 
            MGCV.formula=MGCV.formula, MGCV.select=MGCV.select, 
            MGCVFIX.formula=MGCVFIX.formula, 
            EARTH.formula=EARTH.formula, EARTH.glm=EARTH.glm, 
            RPART.formula=RPART.formula, RPART.xval=RPART.xval, 
            NNET.formula=NNET.formula, NNET.size=NNET.size, NNET.decay=NNET.decay,
            FDA.formula=FDA.formula, 
            SVM.formula=SVM.formula, 
            SVME.formula=SVME.formula, 
            MAHAL.shape=MAHAL.shape, 
            GEODIST.file.name=GEODIST.file.name, RASTER.format=RASTER.format)
        if(is.null(tests$evaluations$MAXENT.T)==F) {output["MAXENT",i] <- tests$evaluations$MAXENT.T@auc}
        if(is.null(tests$evaluations$GBM.T)==F) {output["GBM",i] <- tests$evaluations$GBM.T@auc} 
        if(is.null(tests$evaluations$GBMSTEP.T)==F) {output["GBMSTEP",i] <- tests$evaluations$GBMSTEP.T@auc} 
        if(is.null(tests$evaluations$RF.T)==F) {output["RF",i] <- tests$evaluations$RF.T@auc}
        if(is.null(tests$evaluations$GLM.T)==F) {output["GLM",i] <- tests$evaluations$GLM.T@auc} 
        if(is.null(tests$evaluations$GLMS.T)==F) {output["GLMSTEP",i] <- tests$evaluations$GLMS.T@auc}
        if(is.null(tests$evaluations$GAM.T)==F) {output["GAM",i] <- tests$evaluations$GAM.T@auc} 
        if(is.null(tests$evaluations$GAMS.T)==F) {output["GAMSTEP",i] <- tests$evaluations$GAMS.T@auc}
        if(is.null(tests$evaluations$MGCV.T)==F) {output["MGCV",i] <- tests$evaluations$MGCV.T@auc} 
        if(is.null(tests$evaluations$MGCVF.T)==F) {output["MGCVFIX",i] <- tests$evaluations$MGCVF.T@auc} 
        if(is.null(tests$evaluations$EARTH.T)==F) {output["EARTH",i] <- tests$evaluations$EARTH.T@auc} 
        if(is.null(tests$evaluations$RPART.T)==F) {output["RPART",i] <- tests$evaluations$RPART.T@auc}
        if(is.null(tests$evaluations$NNET.T)==F) {output["NNET",i] <- tests$evaluations$NNET.T@auc} 
        if(is.null(tests$evaluations$FDA.T)==F) {output["FDA",i] <- tests$evaluations$FDA.T@auc}
        if(is.null(tests$evaluations$SVM.T)==F) {output["SVM",i] <- tests$evaluations$SVM.T@auc}
        if(is.null(tests$evaluations$SVME.T)==F) {output["SVME",i] <- tests$evaluations$SVME.T@auc}
        if(is.null(tests$evaluations$BIOCLIM.T)==F) {output["BIOCLIM",i] <- tests$evaluations$BIOCLIM.T@auc}
        if(is.null(tests$evaluations$DOMAIN.T)==F) {output["DOMAIN",i] <- tests$evaluations$DOMAIN.T@auc}
        if(is.null(tests$evaluations$MAHAL.T)==F) {output["MAHAL",i] <- tests$evaluations$MAHAL.T@auc}
        if(is.null(tests$evaluations$GEODIST.T)==F) {output["GEODIST",i] <- tests$evaluations$GEODIST.T@auc}
        if(is.null(tests$evaluations$ENSEMBLE.T)==F) {output["ENSEMBLE",i] <- tests$evaluations$ENSEMBLE.T@auc}
        if(ENSEMBLE.tune == T) {
            output["MAXENT",k+1+i] <- tests$evaluations$STRATEGY.weights["MAXENT"]
            output["GBM",k+1+i] <- tests$evaluations$STRATEGY.weights["GBM"]
            output["GBMSTEP",k+1+i] <- tests$evaluations$STRATEGY.weights["GBMSTEP"]
            output["RF",k+1+i] <- tests$evaluations$STRATEGY.weights["RF"]
            output["GLM",k+1+i] <- tests$evaluations$STRATEGY.weights["GLM"]
            output["GLMSTEP",k+1+i] <- tests$evaluations$STRATEGY.weights["GLMSTEP"]
            output["GAM",k+1+i] <- tests$evaluations$STRATEGY.weights["GAM"]
            output["GAMSTEP",k+1+i] <- tests$evaluations$STRATEGY.weights["GAMSTEP"]
            output["MGCV",k+1+i] <- tests$evaluations$STRATEGY.weights["MGCV"]
            output["MGCVFIX",k+1+i] <- tests$evaluations$STRATEGY.weights["MGCVFIX"]
            output["EARTH",k+1+i] <- tests$evaluations$STRATEGY.weights["EARTH"]
            output["RPART",k+1+i] <- tests$evaluations$STRATEGY.weights["RPART"]
            output["NNET",k+1+i] <- tests$evaluations$STRATEGY.weights["NNET"]
            output["FDA",k+1+i] <- tests$evaluations$STRATEGY.weights["FDA"]
            output["SVM",k+1+i] <- tests$evaluations$STRATEGY.weights["SVM"]
            output["SVME",k+1+i] <- tests$evaluations$STRATEGY.weights["SVME"]
            output["BIOCLIM",k+1+i] <- tests$evaluations$STRATEGY.weights["BIOCLIM"]
            output["DOMAIN",k+1+i] <- tests$evaluations$STRATEGY.weights["DOMAIN"]
            output["MAHAL",k+1+i] <- tests$evaluations$STRATEGY.weights["MAHAL"]
        }
        TestData.all[[i]] <- tests$evaluations$TestData
    }
    output[,k+1] <- rowMeans(output[,c(1:k)], na.rm=T)
    output[is.na(output[,k+1]),(k+1)] <- 0
    if(ENSEMBLE.tune == T) {
        output[,2*k+2] <- rowMeans(output[,c((k+2):(2*k+1))], na.rm=T)
        output[is.na(output[,2*k+2]),(2*k+2)] <- 0
    }
    output.weights <- output[,"MEAN"]
    output.weights <- output.weights[names(output.weights) != "ENSEMBLE"]
# weights multiplied by 100 to avoid problems with default ENSEMBLE.min = 0.7
    if(ENSEMBLE.tune == T) {output.weights <- output.weights * 100}
    output <- output[order(output[,k+1], decreasing=T),]
    cat(paste("\n", "Results of ensemble.test.splits sorted by average AUC for tests T_1 to T_", k, "\n", sep = ""))
    if(ENSEMBLE.tune == T) {
        cat(paste("S_1 to S_", k, " show the weights for the ensemble model with best AUC", "\n", sep = ""))
        cat(paste("column MEAN shows the mean of these weights", "\n\n", sep = ""))
    }
    print(output)
    cat(paste("\n", "suggested input weights for ensemble modelling (percentages based on MEAN column)",  "\n\n", sep = ""))
    print(output.weights)
# test with suggested weights
    output2 <- numeric(length=k)
    names(output2) <- paste("T_", c(1:k), sep="")
    ws <- ensemble.weights(output.weights, decay=1, best=0, min.weight=0, scale=T, multiply=T) 
    for (i in 1:k) {
        TestData <- TestData.all[[i]]
        TestData[,"ENSEMBLE"] <- ws["MAXENT"]*TestData[,"MAXENT"] + ws["GBM"]*TestData[,"GBM"] +
            ws["GBMSTEP"]*TestData[,"GBMSTEP"] + ws["RF"]*TestData[,"RF"] + ws["GLM"]*TestData[,"GLM"] +
            ws["GLMSTEP"]*TestData[,"GLMSTEP"] + ws["GAM"]*TestData[,"GAM"] + ws["GAMSTEP"]*TestData[,"GAMSTEP"] +
            ws["MGCV"]*TestData[,"MGCV"] + ws["MGCVFIX"]*TestData[,"MGCVFIX"] + ws["EARTH"]*TestData[,"EARTH"] +
            ws["RPART"]*TestData[,"RPART"] + ws["NNET"]*TestData[,"NNET"] + ws["FDA"]*TestData[,"FDA"] +
            ws["SVM"]*TestData[,"SVM"] + ws["SVME"]*TestData[,"SVME"] + ws["BIOCLIM"]*TestData[,"BIOCLIM"] +
            ws["DOMAIN"]*TestData[,"DOMAIN"] + ws["MAHAL"]*TestData[,"MAHAL"]
#        TestData[,"ENSEMBLE"] <- trunc(TestData[,"ENSEMBLE"])
        eval1 <- eval2 <- NULL
        TestPres <- as.numeric(TestData[TestData[,"pb"]==1,"ENSEMBLE"])
        TestAbs <- as.numeric(TestData[TestData[,"pb"]==0,"ENSEMBLE"])
        eval1 <- evaluate(p=TestPres, a=TestAbs)
        output2[i] <- eval1@auc
    }
    cat(paste("\n", "AUC for ensemble models based on suggested input weights",  "\n\n", sep = ""))
    print(output2)
    remove(TrainData.c, envir=.BiodiversityR)
    remove(TestData.c, envir=.BiodiversityR)
    if (data.keep == F) {
        cat(paste("\n\n"))
        return(list(table=output, output.weights=output.weights, AUC.with.suggested.weights=output2, call=match.call()))
    }else{
        cat(paste("\n", "(note that output data are integer values representing probabilities multiplied by 1000)",  "\n\n", sep = ""))
        return(list(table=output, output.weights=output.weights, AUC.with.suggested.weights=output2, data=TestData.all, call=match.call()))
    }

}



