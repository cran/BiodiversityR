`ensemble.raster` <- function(
    x=NULL, p=NULL, a=NULL, an=1000, excludep=FALSE, ext=NULL, k=0, pt=NULL, at=NULL, xn=x,
    TrainData=NULL, TestData=NULL,
    ENSEMBLE.rasters=TRUE, 
    layer.drops=NULL,
    models.keep=FALSE,
    RASTER.species.name="Species001", RASTER.stack.name=xn@title, RASTER.format="raster", RASTER.datatype="INT2S", RASTER.NAflag=-32767,
    RASTER.models.overwrite=TRUE,
    threshold.method="spec_sens",
    ENSEMBLE.decay=1, ENSEMBLE.best=0, ENSEMBLE.min=0.7,
    input.weights=NULL, models.list=NULL,
    MAXENT=1, GBM=1, GBMSTEP=1, RF=1, GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=1, MGCV=1, MGCVFIX=0,
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, BIOCLIM=1, DOMAIN=1, MAHAL=1,  
    Yweights="BIOMOD", factors=NULL, dummy.vars=NULL,
    evaluation.strip=TRUE, 
    formulae.defaults=TRUE, maxit=100,
    MAXENT.a=NULL, MAXENT.an=10000, MAXENT.BackData=NULL, MAXENT.path=paste(getwd(), "/models/maxent", sep=""), MAXENT.OLD=NULL,
    GBM.formula=NULL, GBM.n.trees=2001, GBM.OLD=NULL,    
    GBMSTEP.gbm.x=2:(ncol(TrainData.vars)+1), GBMSTEP.tree.complexity=5, GBMSTEP.learning.rate=0.005, 
    GBMSTEP.bag.fraction=0.5, GBMSTEP.step.size=100, GBMSTEP.OLD=NULL,
    RF.formula=NULL, RF.ntree=751, RF.mtry=floor(sqrt(ncol(TrainData.vars))), RF.OLD=NULL,
    GLM.formula=NULL, GLM.family=binomial(link="logit"), GLM.OLD=NULL,
    GLMSTEP.steps=1000, STEP.formula=NULL, GLMSTEP.scope=NULL, GLMSTEP.k=2, GLMSTEP.OLD=NULL,
    GAM.formula=NULL, GAM.family=binomial(link="logit"), GAM.OLD=NULL, 
    GAMSTEP.steps=1000, GAMSTEP.scope=NULL, GAMSTEP.pos=1, GAMSTEP.OLD=NULL, 
    MGCV.formula=NULL, MGCV.select=FALSE, MGCV.OLD=NULL,
    MGCVFIX.formula=NULL, MGCVFIX.OLD=NULL,
    EARTH.formula=NULL, EARTH.glm=list(family=binomial(link="logit"), maxit=maxit), EARTH.OLD=NULL,
    RPART.formula=NULL, RPART.xval=50, RPART.OLD=NULL,
    NNET.formula=NULL, NNET.size=8, NNET.decay=0.01, NNET.OLD=NULL,
    FDA.formula=NULL, FDA.OLD=NULL, 
    SVM.formula=NULL, SVM.OLD=NULL, 
    SVME.formula=NULL, SVME.OLD=NULL,
    BIOCLIM.OLD=NULL, DOMAIN.OLD=NULL, 
    MAHAL.shape=1, MAHAL.OLD=NULL    
)
{
    .BiodiversityR <- new.env()
    if (! require(dismo)) {stop("Please install the dismo package")}
    k <- as.integer(k)
    if (is.null(xn) == T) {
        if (is.null(x) == F) {
            cat(paste("\n", "NOTE: new rasterStack assumed to be equal to the base rasterStack", sep = ""))
            xn <- x
        }else{
            stop("value for parameter xn is missing (RasterStack object)")
        }
    }
    if (is.null(ext) == F) {
        if(length(xn@title) == 0) {xn@title <- "stack1"}
        title.old <- xn@title
        xn <- crop(xn, y=ext, snap="in")
        xn@title <- title.old
    }
# set minimum and maximum values for xn
    for (i in 1:nlayers(xn)) {
        xn[[i]] <- setMinMax(xn[[i]])
    }
#
# check data
    if (is.null(TrainData) == T) {
        if(is.null(x) == T) {stop("value for parameter x is missing (RasterStack object)")}
        if(inherits(x,"RasterStack") == F) {stop("x is not a RasterStack object")}
        if(is.null(p) == T) {stop("presence locations are missing (parameter p)")}
    }
# check TrainData
    if (is.null(TrainData) == F) {
        TrainData <- data.frame(TrainData)
        if (colnames(TrainData)[1] !="pb") {stop("first column for TrainData should be 'pb' containing presence (1) and absence (0) data")}
        if (is.null(x) == F) {
            if (nlayers(x) != (ncol(TrainData)-1)) {
                cat(paste("\n", "WARNING: different number of explanatory variables in rasterStack x and TrainData", sep = ""))
            }
        }
        if ((is.null(xn) == F) && (nlayers(xn) != (ncol(TrainData)-1))) {
            cat(paste("\n", "WARNING: different number of explanatory variables in rasterStack xn and TrainData", sep = ""))
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
                    xn <- dropLayer(xn, which(names(xn) == layer.drops[i]))
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
                j <- which(names(xn) == factors[i])
                xn[[j]] <- raster::as.factor(xn[[j]])
            }
        }
# modify TrainData and xn
    }else{
        if (is.null(layer.drops) == F) {
            vars <- colnames(TrainData)
            layer.drops <- as.character(layer.drops)
            factors <- as.character(factors)
            dummy.vars <- as.character(dummy.vars)
            nd <- length(layer.drops)
            for (i in 1:nd) {
                xn <- dropLayer(xn, which(names(xn) == layer.drops[i]))     
                if (any(vars==layer.drops[i])==FALSE) {
                    cat(paste("\n", "WARNING: variable to exclude '", layer.drops[i], "' not among columns of TrainData", "\n", sep = ""))
                }else{
                    cat(paste("\n", "NOTE: variable ", layer.drops[i], " will not be included as explanatory variable", "\n", sep = ""))
                    TrainData <- TrainData[, which(colnames(TrainData) != layer.drops[i])]
                    if (is.null(TestData) == F) {TestData <- TestData[, which(colnames(TestData) != layer.drops[i])]}
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
        if (identical(RF.ntree, trunc(RF.ntree/2)) == F) {RF.ntree <- RF.ntree + 1}
        environment(RF.formula) <- .BiodiversityR
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
# do not use locations if TrainData is provided
        p <- NULL
        a <- NULL
        pt <- NULL
        at <- NULL
        if(is.null(TestData) == T) {
            TestData <- TrainData
            if (k > 1) {
                groupp <- kfold(TrainData, k=k, by=TrainData[,"pb"])
                TrainData.c <- TrainData[groupp != 1,]
                TestData <- TrainData[groupp == 1,]
                TrainData <- TrainData.c
            }
        }else{
            TestValid <- complete.cases(TestData)
            TestData <- TestData[TestValid,]
        } 
    }else{
        if (is.null(a)==T) {
            if (excludep == T) {
                a <- randomPoints(x[[1]], n=an, p=p, ext=ext, excludep=T)
            }else{
                a <- randomPoints(x[[1]], n=an, p=NULL, ext=ext, excludep=F)
            }        
        }
        if (is.null(pt)==T && is.null(TestData)) {pt <- p}
        if (k > 1 && identical(pt, p) == T) {
            groupp <- kfold(p, k=k)
            pc <- p[groupp != 1,]
            pt <- p[groupp == 1,]
            p <- pc
        }
        if (is.null(at)==T && is.null(TestData)) {at <- a}
        if (k > 1 && identical(at, a) == T) {
            groupa <- kfold(a, k=k)
            ac <- a[groupa != 1,]
            at <- a[groupa == 1,]
            a <- ac
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
    if (is.null(TestData) == F) {
        TestData <- data.frame(TestData)
        if(any(is.na(TestData))) {
            cat(paste("\n", "WARNING: sample units with missing data removed from testing data","\n\n",sep = ""))
        }
        TestValid <- complete.cases(TestData)
        TestData <- TestData[TestValid,]
        if (all(colnames(TestData)!="pb") == T) {stop("one column needed of 'pb' with presence and absence for TestData")} 
    }else{
        TestData <- prepareData(x, p=pt, b=at, factors=factors, xy=FALSE)
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
        TestData <- prepareData(x, p=pt, b=at, factors=factors, xy=FALSE)
    }
# set same factor levels
    if(is.null(factors)==F) {
        for (i in 1:length(factors)) {
            j <- which(names(xn) == factors[i])
            tabulation <- data.frame(freq(xn[[j]]))
            NA.index <- !is.na(tabulation[,"value"])
            tabulation <- tabulation[NA.index,] 
            rasterlevels <- levels(as.factor(tabulation[,1]))
            trainlevels <- levels(TrainData[,factors[i]])
            testlevels <- levels(TestData[,factors[i]])
            uniquelevels <- unique(c(rasterlevels, trainlevels, testlevels))
            levels(TrainData[,factors[i]]) <- uniquelevels
            levels(TestData[,factors[i]]) <- uniquelevels
        } 
    }
# check if TestData is different from TrainData
    no.tests <- FALSE
    if (identical(TrainData, TestData) == T) {no.tests <- TRUE}
#
    if(models.keep==T) {
        models <- list(MAXENT=NULL, GBM=NULL, GBMSTEP=NULL, RF=NULL, GLM=NULL, 
            GLMSTEP=NULL, GAM=NULL, GAMSTEP=NULL, MGCV=NULL, MGCVFIX=NULL, EARTH=NULL, RPART=NULL, 
            NNET=NULL, FDA=NULL, SVM=NULL, SVME=NULL, BIOCLIM=NULL, DOMAIN=NULL, MAHAL=NULL,
            formulae=NULL, TrainData=NULL, TestData=NULL, p=NULL, a=NULL, pt=NULL, at=NULL, MAXENT.BackData=NULL)
        models$TrainData <- TrainData
        models$p <- p
        models$a <- a
        models$pt <- pt  
        models$at <- at
        if (no.tests == F) {models$TestData <- TestData}
    }else{
        models <- NULL
    }
# Data frames for distance-based methods and SVME
    TrainData.vars <- TrainData[,colnames(TrainData) != "pb"]
    assign("TrainData.vars", TrainData.vars, envir=.BiodiversityR)
    TrainData.pres <- TrainData[TrainData[,"pb"]==1,]
    TrainData.pres <- TrainData.pres[,colnames(TrainData.pres) != "pb"]
    assign("TrainData.pres", TrainData.pres, envir=.BiodiversityR)
    TestData.vars <- TestData[,colnames(TestData) != "pb"]
    assign("TestData.vars", TestData.vars, envir=.BiodiversityR)
#
# separate data set to calibrate MAXENT
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
        if (is.null(layer.drops) == F) {
            nd <- length(layer.drops)
            for (i in 1:nd) {     
                MAXENT.BackData <- MAXENT.BackData[, which(colnames(MAXENT.BackData) != layer.drops[i])]
            }
        }
        TestValid <- complete.cases(MAXENT.BackData)
        MAXENT.BackData <- MAXENT.BackData[TestValid,]
        if(models.keep == T) {models$MAXENT.BackData <- MAXENT.BackData}
        if (identical(colnames(TrainData.pres), colnames(MAXENT.BackData)) == F) {
                cat(paste("\n", "WARNING: MAXENT.BackData has different (sequence of) colnames than TrainData", sep = ""))
                cat(paste("\n", "MAXENT model will not be calibrated", "\n", sep = "")) 
                MAXENT <- 0
        }else{
            MAXENT.TrainData <- rbind(TrainData.pres, MAXENT.BackData)
            MAXENT.pa <- c(rep(1, nrow(TrainData.pres)), rep(0, nrow(MAXENT.BackData)))
            assign("MAXENT.TrainData", MAXENT.TrainData, envir=.BiodiversityR)
            assign("MAXENT.pa", MAXENT.pa, envir=.BiodiversityR)            
        }
    }
#
    modelresults <- data.frame(array(dim=c(nrow(TrainData), 20), 0))
    colnames(modelresults) <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL", "ENSEMBLE")
    TrainData <- cbind(TrainData, modelresults)
    modelresults <- data.frame(array(dim=c(nrow(TestData), 20), 0))
    colnames(modelresults) <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL", "ENSEMBLE")
    TestData <- cbind(TestData, modelresults)
    assign("TestData", TestData, envir=.BiodiversityR)
    assign("TrainData", TrainData, envir=.BiodiversityR)
#
    ws <- as.numeric(c(MAXENT, GBM, GBMSTEP, RF, GLM, GLMSTEP, GAM, GAMSTEP, MGCV, MGCVFIX, 
        EARTH, RPART, NNET, FDA, SVM, SVME, BIOCLIM, DOMAIN, MAHAL))
    names(ws) <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX", 
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL")
    if((length(ENSEMBLE.decay) == 1) && (length(ENSEMBLE.best) == 1) && (length(ENSEMBLE.min) == 1)) {
        ws <- ensemble.weights(ws, decay=ENSEMBLE.decay, 
            best=ENSEMBLE.best, min.weight=ENSEMBLE.min, scale=TRUE, multiply=TRUE)
            cat(paste("\n", "Weights for ensemble forecasting", sep = ""))
            cat(paste("\n", "(minimum weights for ensemble forecasting was: ", ENSEMBLE.min, ")\n", sep = ""))
            print(ws)
    }
    thresholds <- ws
    prediction.failures <- FALSE
#
    evaluations <- array(NA, dim=c(20, 2))
    rownames(evaluations) <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL", "ENSEMBLE")
    colnames(evaluations) <- c("calibration","test")
#
    if(evaluation.strip==T) {
        ResponseData <- evaluation.strip.data(xn, ext=ext, vars=names(xn), factors=factors)
        evaluations <- list(evaluations=evaluations, evaluation.strip=ResponseData)
        vars <- names(xn)
        ResponseData.vars <- ResponseData[,vars]
        assign("ResponseData.vars", ResponseData.vars, envir=.BiodiversityR)
        assign("ResponseData", ResponseData, envir=.BiodiversityR)
    }else{    
        evaluations <- list(evaluations=evaluations)
    }
    dir.create("models", showWarnings = F)
    dir.create("ensembles", showWarnings = F)
    rasterfull <- paste("ensembles/", RASTER.species.name, "_ENSEMBLE_", xn@title , sep="")
    rastercount <- paste("ensembles/", RASTER.species.name, "_COUNT_", xn@title , sep="")
    rasterpresence <- paste("ensembles/", RASTER.species.name, "_PRESENCE_", xn@title, sep="")
    RASTER.species.orig <- RASTER.species.name
    if (RASTER.models.overwrite==T) {
        RASTER.species.name <- "working"
    }else{
        RASTER.species.name <- paste(RASTER.species.name, "_", RASTER.stack.name, sep="")
    }
    Yweights1 <- Yweights
    if (Yweights == "BIOMOD") {
        #have equal weight of presence vs. background
        Yweights1 <- numeric(length = nrow(TrainData))
        pres <- length(TrainData[, 1] == 1)
        abs <- length(TrainData[, 1] == 0)
        Yweights1[which(TrainData[, 1] == 1)] <- 1
        Yweights1[which(TrainData[, 1] == 0)] <- pres/abs
    }
    if (Yweights == "equal") {
        Yweights1 <- numeric(length = nrow(TrainData))
        Yweights1[] <- 1
    }
    Yweights <- Yweights1
    assign("Yweights", Yweights, envir=.BiodiversityR)
    cat(paste("\n", "Start of modelling for organism: ", RASTER.species.orig, "\n", sep = ""))
    ensemble.statistics <- NULL
    if (ENSEMBLE.rasters == T) {
        cat(paste("ensemble raster layers will be saved in folder ", getwd(), "/ensembles", "\n\n", sep = ""))
        statistics.names <- c("n.models", "ensemble.threshold", "ensemble.min", "ensemble.max", "count.min", "count.max") 
        ensemble.statistics <- numeric(6)
        names(ensemble.statistics) <- statistics.names
    }
# start raster layer creations
    if (ws["MAXENT"] > 0) {
        # Put the file 'maxent.jar' in the 'java' folder of dismo
        # the file 'maxent.jar' can be obtained from from http://www.cs.princeton.edu/~schapire/maxent/.
        jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
        eval1 <- eval2 <- results <- pmaxent <- TestPres <- TestAbs <- TrainPres <- TrainAbs <- NULL
        if(is.null(MAXENT.OLD) == T) {
            tryCatch(results <- maxent(x=MAXENT.TrainData, p=MAXENT.pa, factors=factors, path=MAXENT.path),
                error= function(err) {print(paste("MAXENT calibration failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["MAXENT"] <- -1
            }
        }else{ results <- MAXENT.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with calibration data for MAXENT", sep = ""))
            cat(paste("\n", "(tested with calibration data of other algorithms)", "\n\n",sep = ""))
            TrainData[,"MAXENT"] <- dismo::predict(object=results, x=TrainData.vars)
            TrainData[,"MAXENT"] <- trunc(1000*TrainData[,"MAXENT"])
            TrainPres <- TrainData[TrainData[,"pb"]==1,"MAXENT"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"MAXENT"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["MAXENT"] <- threshold(eval1, threshold.method)
            evaluations$evaluations["MAXENT",1] <-  eval1@auc
            if (no.tests == F) {
                cat(paste("\n", "MAXENT evaluation with test data","\n\n", sep = ""))
                TestData[,"MAXENT"] <- dismo::predict(object=results, x=TestData.vars)
                TestData[,"MAXENT"] <- trunc(1000*TestData[,"MAXENT"])
                TestPres <- TestData[TestData[,"pb"]==1,"MAXENT"]
                TestAbs <- TestData[TestData[,"pb"]==0,"MAXENT"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("MAXENT evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                    evaluations$evaluations["MAXENT",2] <-  eval2@auc
                }
            }
            if (models.keep==T) {models$MAXENT <- results}
            if (evaluation.strip == T) {
                ResponseData$MAXENT <- dismo::predict(object=results, x=ResponseData.vars)           
            }
            fullname <- paste("models/", RASTER.species.name, "_MAXENT", sep="")
            tryCatch(pmaxent <- raster::predict(object=results, x=xn, ext=ext, na.rm=TRUE, 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("MAXENT prediction failed"))},
                silent=T)
            if (is.null(pmaxent) == F) {
                pmaxent <- trunc(1000*pmaxent)
                writeRaster(x=pmaxent, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            }else{
                cat(paste("\n", "WARNING: MAXENT prediction failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["MAXENT"] <- -1 
            }
        }
    }
    if (ws["GBM"] > 0) {
#       require(gbm, quietly=T) 
        eval1 <- eval2 <- results <- pgbm <- TestPres <- TestAbs <- TrainPres <- TrainAbs <- NULL
        if(is.null(GBM.OLD) == T) {       
            tryCatch(results <- gbm(formula=GBM.formula, data=TrainData, weights=Yweights, distribution="bernoulli", 
                    interaction.depth=7, shrinkage=0.001, bag.fraction=0.5, train.fraction=1, 
                    n.trees=GBM.n.trees, verbose=F, cv.folds=5),
                error= function(err) {print(paste("GBM calibration failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["GBM"] <- -1
            }
        }else{ results <- GBM.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with calibration data for GBM",  "\n\n", sep = ""))
            TrainData[,"GBM"] <- predict(object=results, newdata=TrainData.vars, n.trees=results$n.trees, type="response")
            TrainData[,"GBM"] <- trunc(1000*TrainData[,"GBM"])
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GBM"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GBM"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["GBM"] <- threshold(eval1, threshold.method)
            evaluations$evaluations["GBM",1] <-  eval1@auc
            if (no.tests == F) {
                cat(paste("\n", "GBM evaluation with test data","\n\n", sep = ""))
                TestData[,"GBM"] <- predict(object=results, newdata=TestData.vars, n.trees=results$n.trees, type="response")
                TestData[,"GBM"] <- trunc(1000*TestData[,"GBM"])
                TestPres <- TestData[TestData[,"pb"]==1,"GBM"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GBM"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("GBM evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                    evaluations$evaluations["GBM",2] <-  eval2@auc
                }
            }
            if (models.keep==T) {
                models$GBM <- results
                models$formulae$GBM.formula <- GBM.formula
            }
            if (evaluation.strip == T) {
                ResponseData$GBM <- predict(results, newdata=ResponseData.vars, n.trees=results$n.trees, type="response")
            }
            fullname <- paste("models/", RASTER.species.name, "_GBM", sep="")
            tryCatch(pgbm <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, n.trees=results$n.trees, type="response", 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("GBM prediction failed"))},
                silent=T)
            if (is.null(pgbm) == F) {
                pgbm <- trunc(1000*pgbm)
                writeRaster(x=pgbm, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            }else{
                cat(paste("\n", "WARNING: GBM prediction failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["GBM"] <- -1 
            }
        }    
#       detach(package:gbm)
    }
    if (ws["GBMSTEP"] > 0) {
#       require(gbm, quietly=T) 
        eval1 <- eval2 <- results <- pgbms <- TestPres <- TestAbs <- TrainPres <- TrainAbs <- NULL
        if(is.null(GBMSTEP.OLD) == T) {       
            tryCatch(results <- gbm.step(data=TrainData, gbm.y=1, gbm.x=GBMSTEP.gbm.x, family="bernoulli",
                    site.weights=Yweights, tree.complexity=GBMSTEP.tree.complexity, learning.rate = GBMSTEP.learning.rate, 
                    bag.fraction=GBMSTEP.bag.fraction, step.size=GBMSTEP.step.size, verbose=F, silent=T, plot.main=F),
                error= function(err) {print(paste("stepwise GBM calibration failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["GBMSTEP"] <- -1
            }
        }else{ results <- GBMSTEP.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", "stepwise GBM trees (target > 1000)", "\n", sep = ""))
            print(results$n.trees)
            cat(paste("\n", RASTER.species.orig, ": Evaluation with calibration data for stepwise GBM",  "\n\n", sep = ""))
            TrainData[,"GBMSTEP"] <- predict(object=results, newdata=TrainData.vars, n.trees=results$n.trees, type="response")
            TrainData[,"GBMSTEP"] <- trunc(1000*TrainData[,"GBMSTEP"])
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GBMSTEP"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GBMSTEP"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["GBMSTEP"] <- threshold(eval1, threshold.method)
            evaluations$evaluations["GBMSTEP",1] <-  eval1@auc
            if (no.tests == F) {
                cat(paste("\n", "stepwise GBM evaluation with test data","\n\n", sep = ""))
                TestData[,"GBMSTEP"] <- predict(object=results, newdata=TestData.vars, n.trees=results$n.trees, type="response")
                TestData[,"GBMSTEP"] <- trunc(1000*TestData[,"GBMSTEP"])
                TestPres <- TestData[TestData[,"pb"]==1,"GBMSTEP"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GBMSTEP"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("stepwise GBM evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                    evaluations$evaluations["GBMSTEP",2] <-  eval2@auc
                }
            }
            if (models.keep==T) {models$GBMSTEP <- results}
            if (evaluation.strip == T) {
                ResponseData$GBMSTEP <- predict(results, newdata=ResponseData.vars, n.trees=results$n.trees, type="response")
            }
            fullname <- paste("models/", RASTER.species.name, "_GBMSTEP", sep="")
            tryCatch(pgbms <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, n.trees=results$n.trees, type="response", 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("stepwise GBM prediction failed"))},
                silent=T)
            if (is.null(pgbms) == F) {
                pgbms <- trunc(1000*pgbms)
                writeRaster(x=pgbms, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            }else{
                cat(paste("\n", "WARNING: stepwise GBM prediction failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["GBMSTEP"] <- -1 
            }
        }
#       detach(package:gbm)
    }
    if (ws["RF"] > 0) {
#       require(randomForest, quietly=T)
        eval1 <- eval2 <- results <- prf <- TestPres <- TestAbs <- TrainPres <- TrainAbs <- NULL
        if(is.null(RF.OLD) == T) {        
            tryCatch(results <- randomForest(formula=RF.formula, ntree=RF.ntree, mtry=RF.mtry, data=TrainData, na.action=na.omit),
                error= function(err) {print(paste("random forest calibration failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["RF"] <- -1
            }
        }else{ results <- RF.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with calibration data for RF",  "\n\n", sep = ""))
            TrainData[,"RF"] <- predict(object=results, newdata=TrainData.vars, type="response")
            TrainData[,"RF"] <- trunc(1000*TrainData[,"RF"])
            TrainPres <- TrainData[TrainData[,"pb"]==1,"RF"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"RF"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["RF"] <- threshold(eval1, threshold.method)
            evaluations$evaluations["RF",1] <-  eval1@auc
            if (no.tests == F) {
                cat(paste("\n", "RF evaluation with test data","\n\n", sep = ""))
                TestData[,"RF"] <- predict(object=results, newdata=TestData.vars, type="response")
                TestData[,"RF"] <- trunc(1000*TestData[,"RF"])
                TestPres <- TestData[TestData[,"pb"]==1,"RF"]
                TestAbs <- TestData[TestData[,"pb"]==0,"RF"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("RF evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                    evaluations$evaluations["RF",2] <-  eval2@auc
                }
            }
            if (models.keep==T) {
                models$RF <- results
                models$formulae$RF.formula <- RF.formula
            }
            if (evaluation.strip == T) {
                ResponseData$RF <- predict(results, newdata=ResponseData.vars, type="response")
            }
            fullname <- paste("models/", RASTER.species.name, "_RF", sep="")
            tryCatch(prf <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("random forest prediction failed"))},
                silent=T)
            if (is.null(prf) == F) {
                prf <- trunc(1000*prf)
                writeRaster(x=prf, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            }else{
                cat(paste("\n", "WARNING: random forest prediction failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["RF"] <- -1 
            }
        }
#       detach(package:randomForest)
    } 
    if (ws["GLM"] > 0) {
        eval1 <- eval2 <- results <- pglm <- TestPres <- TestAbs <- TrainPres <- TrainAbs <- NULL
        if(is.null(GLM.OLD) == T) { 
            tryCatch(results <- glm(formula=GLM.formula, family=GLM.family, data=TrainData, weights=Yweights, control=glm.control(maxit=maxit)),
                error= function(err) {print(paste("GLM calibration failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["GLM"] <- -1
            }
        }else{ results <- GLM.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with calibration data for GLM",  "\n\n", sep = ""))
            TrainData[,"GLM"] <- predict(object=results, newdata=TrainData.vars, type="response")
            TrainData[,"GLM"] <- trunc(1000*TrainData[,"GLM"])
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GLM"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GLM"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["GLM"] <- threshold(eval1, threshold.method)
            evaluations$evaluations["GLM",1] <-  eval1@auc
            if (no.tests == F) {
                cat(paste("\n", "GLM evaluation with test data","\n\n", sep = ""))
                TestData[,"GLM"] <- predict(object=results, newdata=TestData.vars, type="response")
                TestData[,"GLM"] <- trunc(1000*TestData[,"GLM"])
                TestPres <- TestData[TestData[,"pb"]==1,"GLM"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GLM"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("GLM evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                    evaluations$evaluations["GLM",2] <-  eval2@auc
                }
            }
            if (models.keep==T) {
                models$GLM <- results
                models$formulae$GLM.formula <- GLM.formula
            }
            if (evaluation.strip == T) {
                ResponseData$GLM <- predict(results, newdata=ResponseData.vars, type="response")
            }
            fullname <- paste("models/", RASTER.species.name, "_GLM", sep="")
            tryCatch(pglm <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, type="response", 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("GLM prediction failed"))},
                silent=T)
            if (is.null(pglm) == F) {
                pglm <- trunc(1000*pglm)
                writeRaster(x=pglm, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            }else{
                cat(paste("\n", "WARNING: GLM prediction failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["GLM"] <- -1 
            }
        }
    }
    if (ws["GLMSTEP"] > 0) {
#        require(MASS, quietly=T)
        eval1 <- eval2 <- results <- results2 <- pglms <- TestPres <- TestAbs <- TrainPres <- TrainAbs <- NULL
        if(is.null(GLMSTEP.OLD) == T) {
            tryCatch(results <- glm(formula=STEP.formula, family=GLM.family, data=TrainData, weights=Yweights, control=glm.control(maxit=maxit)),
                error= function(err) {print(paste("first step of stepwise GLM calibration failed"))},
                silent=T)
            tryCatch(results2 <- stepAIC(results, scope=GLMSTEP.scope, direction="both", trace=F, steps=GLMSTEP.steps, k=GLMSTEP.k),
                error= function(err) {print(paste("stepwise GLM calibration failed"))},
                silent=T)
            if(is.null(results2) == T) {
                prediction.failures <- TRUE
                ws["GLMSTEP"] <- -1
            }
        }else{ results2 <- GLMSTEP.OLD}
        if (is.null(results2) == F) {
            results <- results2
            cat(paste("\n", "stepwise GLM formula","\n\n",sep = ""))
            print(formula(results))
            cat(paste("\n", RASTER.species.orig, ": Evaluation with calibration data for stepwise GLM", "\n\n", sep = ""))
            TrainData[,"GLMSTEP"] <- predict(object=results, newdata=TrainData.vars, type="response")
            TrainData[,"GLMSTEP"] <- trunc(1000*TrainData[,"GLMSTEP"])
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GLMSTEP"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GLMSTEP"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["GLMSTEP"] <- threshold(eval1, threshold.method)
            evaluations$evaluations["GLMSTEP",1] <-  eval1@auc
            if (no.tests == F) {
                cat(paste("\n", "stepwise GLM evaluation with test data","\n\n", sep = ""))
                TestData[,"GLMSTEP"] <- predict(object=results, newdata=TestData.vars, type="response")
                TestData[,"GLMSTEP"] <- trunc(1000*TestData[,"GLMSTEP"])
                TestPres <- TestData[TestData[,"pb"]==1,"GLMSTEP"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GLMSTEP"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("stepwise GLM evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                    evaluations$evaluations["GLMSTEP",2] <-  eval2@auc
                }
            }
            if (models.keep==T) {
                models$GLMSTEP <- results
                models$formulae$STEP.formula <- STEP.formula
                models$formulae$GLMSTEP.scope <- GLMSTEP.scope
            }
            if (evaluation.strip == T) {
                ResponseData$GLMSTEP <- predict(results, newdata=ResponseData.vars, type="response")
            }
            fullname <- paste("models/", RASTER.species.name, "_GLMSTEP", sep="")
            tryCatch(pglms <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, type="response", 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("stepwise GLM prediction failed"))},
                silent=T)
            if (is.null(pglms) == F) {
                pglms <- trunc(1000*pglms)
                writeRaster(x=pglms, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            }else{
                cat(paste("\n", "WARNING: stepwise GLM prediction failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["GLMSTEP"] <- -1 
            } 
        }
#       detach(package:MASS)
    }
    if (ws["GAM"] > 0) {
        cat(paste("\n\n"))
        require(gam, quietly=T)
        eval1 <- eval2 <- results <- pgam <- TestPres <- TestAbs <- TrainPres <- TrainAbs <- NULL
        if(is.null(GAM.OLD) == T) {
            tryCatch(results <- gam(formula=GAM.formula, family=GAM.family, data=TrainData, weights=Yweights, control=gam.control(maxit=maxit, bf.maxit=50)),
                error= function(err) {print(paste("GAM calibration (gam package) failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["GAM"] <- -1
            }
        }else{ results <- GAM.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with calibration data for GAM (gam package)", "\n\n", sep = ""))
            TrainData[,"GAM"] <- predict(object=results, newdata=TrainData.vars, type="response")
            TrainData[,"GAM"] <- trunc(1000*TrainData[,"GAM"])
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GAM"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GAM"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["GAM"] <- threshold(eval1, threshold.method)
            evaluations$evaluations["GAM",1] <-  eval1@auc
            if (no.tests == F) {
                cat(paste("\n", "GAM evaluation with test data","\n\n", sep = ""))
                TestData[,"GAM"] <- predict(object=results, newdata=TestData.vars, type="response")
                TestData[,"GAM"] <- trunc(1000*TestData[,"GAM"])
                TestPres <- TestData[TestData[,"pb"]==1,"GAM"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GAM"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("GAM evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                    evaluations$evaluations["GAM",2] <-  eval2@auc
                }
            }
            if (models.keep==T) {
                models$GAM <- results
                models$formulae$GAM.formula <- GAM.formula
            }
            if (evaluation.strip == T) {
                ResponseData$GAM <- predict(results, newdata=ResponseData.vars, type="response")
            }
            fullname <- paste("models/", RASTER.species.name, "_GAM", sep="")
            tryCatch(pgam <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, type="response", 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("GAM prediction (gam package) failed"))},
                silent=T)
            if (is.null(pgam) == F) {
                pgam <- trunc(1000*pgam)
                writeRaster(x=pgam, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            }else{
                cat(paste("\n", "WARNING: GAM prediction (gam package) failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["GAM"] <- -1 
            } 
        }
        detach(package:gam)
    }
    if (ws["GAMSTEP"] > 0) {
        cat(paste("\n\n"))
        require(gam, quietly=T)
        eval1 <- eval2 <- results <- results2 <- pgams <- TestPres <- TestAbs <- TrainPres <- TrainAbs <- NULL
        if(is.null(GAMSTEP.OLD) == T) {
            tryCatch(results <- gam(formula=STEP.formula, family=GAM.family, data=TrainData, weights=Yweights, control=gam.control(maxit=maxit, bf.maxit=50)), 
                error= function(err) {print(paste("first step of stepwise GAM calibration (gam package) failed"))},
                silent=T)
            assign("TrainData", TrainData, pos=GAMSTEP.pos)
            assign("GAM.family", GAM.family, pos=GAMSTEP.pos)
            assign("maxit", maxit, pos=GAMSTEP.pos)   
            tryCatch(results2 <- step.gam(results, scope=GAMSTEP.scope, direction="both", trace=F, steps=GAMSTEP.steps), 
                error= function(err) {print(paste("stepwise GAM calibration (gam package) failed"))},
                silent=T)
            remove(TrainData, pos=GAMSTEP.pos)
            remove(GAM.family, pos=GAMSTEP.pos)
            remove(maxit, pos=GAMSTEP.pos)
            if(is.null(results2) == T) {
                prediction.failures <- TRUE
                ws["GAMSTEP"] <- -1
            }
        }else{ results2 <- GAMSTEP.OLD}
        if (is.null(results2) == F) {
            results <- results2
            cat(paste("\n", "stepwise GAM formula (gam package)","\n\n",sep = ""))
            print(formula(results))
            cat(paste("\n", RASTER.species.orig, ": Evaluation with calibration data for stepwise GAM (gam package)", "\n\n", sep = ""))
            TrainData[,"GAMSTEP"] <- predict(object=results, newdata=TrainData.vars, type="response")
            TrainData[,"GAMSTEP"] <- trunc(1000*TrainData[,"GAMSTEP"])
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GAMSTEP"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GAMSTEP"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["GAMSTEP"] <- threshold(eval1, threshold.method) 
            evaluations$evaluations["GAMSTEP",1] <-  eval1@auc
            if (no.tests == F) {
                cat(paste("\n", "stepwise GAM evaluation with test data","\n\n", "\n", sep = ""))
                TestData[,"GAMSTEP"] <- predict(object=results, newdata=TestData.vars, type="response")
                TestData[,"GAMSTEP"] <- trunc(1000*TestData[,"GAMSTEP"])
                TestPres <- TestData[TestData[,"pb"]==1,"GAMSTEP"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GAMSTEP"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("stepwise GAM evaluation (gam package) failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                    evaluations$evaluations["GAMSTEP",2] <-  eval2@auc
                }
            }
            if (models.keep==T) {
                models$GAMSTEP <- results
                models$formulae$STEP.formula <- STEP.formula
                models$formulae$GAMSTEP.scope <- GAMSTEP.scope
            }
            if (evaluation.strip == T) {
                ResponseData$GAMSTEP <- predict(results, newdata=ResponseData.vars, type="response")
            }
            fullname <- paste("models/", RASTER.species.name, "_GAMSTEP", sep="")
            tryCatch(pgams <- raster::predict(object=xn, model=results, ext=ext, type="response", na.rm=TRUE, 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("stepwise GAM prediction (gam package) failed"))},
                silent=T)
            if (is.null(pgams) == F) {
                pgams <- trunc(1000*pgams)
                writeRaster(x=pgams, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            }else{
                cat(paste("\n", "WARNING: stepwise GAM prediction (gam package) failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["GAMSTEP"] <- -1 
            } 
        }
        detach(package:gam)
    }
    if (ws["MGCV"] > 0) {
        cat(paste("\n\n"))
        eval1 <- eval2 <- results <- pmgcv <- TestPres <- TestAbs <- TrainPres <- TrainAbs <- NULL
        require(mgcv, quietly=T)
        if(is.null(MGCV.OLD) == T) {
            tryCatch(results <- gam(formula=MGCV.formula, family=GAM.family, data=TrainData, weights=Yweights, 
                        select=MGCV.select, control=gam.control(maxit=maxit)),
                    error= function(err) {print(paste("GAM calibration (mgcv package) failed"))},
                    silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["MGCV"] <- -1
            }
        }else{ results <- MGCV.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with calibration data for GAM (mgcv package)", "\n\n", sep = ""))
            TrainData[,"MGCV"] <- predict.mgcv(object=results, newdata=TrainData.vars)
            TrainData[,"MGCV"] <- trunc(1000*TrainData[,"MGCV"])
            TrainPres <- TrainData[TrainData[,"pb"]==1,"MGCV"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"MGCV"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["MGCV"] <- threshold(eval1, threshold.method) 
            evaluations$evaluations["MGCV",1] <-  eval1@auc
            if (no.tests == F) {
                cat(paste("\n", "MGCV evaluation with test data","\n\n", sep = ""))
                TestData[,"MGCV"] <- predict.mgcv(object=results, newdata=TestData.vars)
                TestData[,"MGCV"] <- trunc(1000*TestData[,"MGCV"])
                TestPres <- TestData[TestData[,"pb"]==1,"MGCV"]
                TestAbs <- TestData[TestData[,"pb"]==0,"MGCV"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("GAM evaluation (mgcv package) failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                    evaluations$evaluations["MGCV",2] <-  eval2@auc
                }
            }
            if (models.keep==T) {
                models$MGCV <- results
                models$formulae$MGCV.formula <- MGCV.formula
            }
            if (evaluation.strip == T) {
                ResponseData$MGCV <- predict(results, newdata=ResponseData.vars, type="response")
            }
            fullname <- paste("models/", RASTER.species.name, "_MGCV", sep="")
            tryCatch(pmgcv <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, type="response", 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("GAM prediction (mgcv package) failed"))},
                silent=T)
            if (is.null(pmgcv) == F) {
                pmgcv <- trunc(1000*pmgcv)
                writeRaster(x=pmgcv, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            }else{
                cat(paste("\n", "WARNING: GAM prediction (mgcv package) failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["MGCV"] <- -1 
            } 
        }
        detach(package:mgcv)
    }
    if (ws["MGCVFIX"] > 0) {
        cat(paste("\n\n"))
        eval1 <- eval2 <- results <- pmgcvF <- TestPres <- TestAbs <- TrainPres <- TrainAbs <- NULL
        require(mgcv, quietly=T)
        if(is.null(MGCVFIX.OLD) == T) {
            tryCatch(results <- gam(formula=MGCVFIX.formula, family=GAM.family, data=TrainData, weights=Yweights, select=FALSE, control=gam.control(maxit=maxit)),
                error= function(err) {print(paste("GAM calibration with fixed d.f. regression splines (mgcv package) failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["MGCVFIX"] <- -1
            }
        }else{ results <- MGCVFIX.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with calibration data for MGCVFIX (mgcv package)", "\n\n", sep = ""))
            TrainData[,"MGCVFIX"] <- predict.mgcv(object=results, newdata=TrainData.vars)
            TrainData[,"MGCVFIX"] <- trunc(1000*TrainData[,"MGCVFIX"])
            TrainPres <- TrainData[TrainData[,"pb"]==1,"MGCVFIX"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"MGCVFIX"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["MGCVFIX"] <- threshold(eval1, threshold.method)
            evaluations$evaluations["MGCVFIX",1] <-  eval1@auc
            if (no.tests == F) {
                cat(paste("\n", "MGCVFIX evaluation with test data","\n\n", sep = ""))
                TestData[,"MGCVFIX"] <- predict.mgcv(object=results, newdata=TestData)
                TestData[,"MGCVFIX"] <- trunc(1000*TestData[,"MGCVFIX"])
                TestPres <- TestData[TestData[,"pb"]==1,"MGCVFIX"]
                TestAbs <- TestData[TestData[,"pb"]==0,"MGCVFIX"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("GAMFIX evaluation (mgcv package) failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                    evaluations$evaluations["MGCVFIX",2] <-  eval2@auc
                }
            }
            if (models.keep==T) {
                models$MGCVFIX <- results
                models$formulae$MGCVFIX.formula <- MGCVFIX.formula
            }
            if (evaluation.strip == T) {
                ResponseData$MGCVFIX <- predict(results, newdata=ResponseData.vars, type="response")
            }
            fullname <- paste("models/", RASTER.species.name, "_MGCVFIX", sep="")
            tryCatch(pmgcvf <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, type="response", 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("MGCVFIX prediction (mgcv package) failed"))},
                silent=T)
            if (is.null(pmgcvf) == F) {
                pmgcvf <- trunc(1000*pmgcvf)
                writeRaster(x=pmgcvf, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            }else{
                cat(paste("\n", "WARNING: GAM prediction (mgcv package) failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["MGCVFIX"] <- -1 
            } 
        }
        detach(package:mgcv)
    }
    if (!is.null(factors) && EARTH > 0) {
        cat(paste("\n", "MARS model (earth package) with factors may require explicit dummy variables", "\n", sep=""))
    } 
    if (ws["EARTH"] > 0) {
#        require(earth, quietly=T)
        eval1 <- eval2 <- results <- pearth <- TestPres <- TestAbs <- TrainPres <- TrainAbs <- NULL
        if(is.null(EARTH.OLD) == T) {
            tryCatch(results <- earth(formula=EARTH.formula, glm=EARTH.glm, data=TrainData, degree=2),
                error= function(err) {print(paste("MARS calibration (earth package) failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["EARTH"] <- -1
            }
        }else{ results <- EARTH.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with calibration data for MARS (earth package)", "\n\n", sep = ""))
            TrainData[,"EARTH"] <- predict.earth2(object=results, newdata=TrainData.vars)
            TrainData[,"EARTH"] <- trunc(1000*TrainData[,"EARTH"])
            TrainPres <- TrainData[TrainData[,"pb"]==1,"EARTH"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"EARTH"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["EARTH"] <- threshold(eval1, threshold.method) 
            evaluations$evaluations["EARTH",1] <-  eval1@auc
            if (no.tests == F) {
                cat(paste("\n", "MARS evaluation with test data","\n\n", sep = ""))
                TestData[,"EARTH"] <- predict.earth2(object=results, newdata=TestData.vars)
                TestData[,"EARTH"] <- trunc(1000*TestData[,"EARTH"])
                TestPres <- TestData[TestData[,"pb"]==1,"EARTH"]
                TestAbs <- TestData[TestData[,"pb"]==0,"EARTH"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("MARS evaluation (earth package) failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                    evaluations$evaluations["EARTH",2] <-  eval2@auc
                }
            }
            if (models.keep==T) {
                models$EARTH <- results
                models$formulae$EARTH.formula <- EARTH.formula
            }
            if (evaluation.strip == T) {
                ResponseData$EARTH <- predict.earth2(results, newdata=ResponseData.vars, type="response")
            }
            fullname <- paste("models/", RASTER.species.name, "_EARTH", sep="")
            tryCatch(pearth <- raster::predict(object=xn, model=results, fun=predict.earth2, ext=ext, na.rm=TRUE, type="response", 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("MARS prediction (earth package) failed"))},
                silent=T)
            if (is.null(pearth) == F) {
                pearth <- trunc(1000*pearth)
                writeRaster(x=pearth, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            }else{
                cat(paste("\n", "WARNING: MARS prediction (earth package) failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["EARTH"] <- -1 
            } 
        }
#       detach(package:earth)
    }
    if (ws["RPART"] > 0) {
#        require(rpart, quietly=T)
        eval1 <- eval2 <- results <- prpart <- TestPres <- TestAbs <- TrainPres <- TrainAbs <- NULL
        if(is.null(RPART.OLD) == T) {
            tryCatch(results <- rpart(formula=RPART.formula, data=TrainData, weights=Yweights,
                    control=rpart.control(xval=RPART.xval, minbucket=5, minsplit=5, cp=0.001, maxdepth=25)),
                error= function(err) {print(paste("RPART calibration failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["RPART"] <- -1
            }
        }else{ results <- RPART.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with calibration data for RPART", "\n\n", sep = ""))
            TrainData[,"RPART"] <- predict(object=results, newdata=TrainData.vars, type="prob")[,2]
            TrainData[,"RPART"] <- trunc(1000*TrainData[,"RPART"])
            TrainPres <- TrainData[TrainData[,"pb"]==1,"RPART"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"RPART"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["RPART"] <- threshold(eval1, threshold.method) 
            evaluations$evaluations["RPART",1] <-  eval1@auc
            if (no.tests == F) {
                cat(paste("\n", "RPART evaluation with test data","\n\n", sep = ""))
                TestData[,"RPART"] <- predict(object=results, newdata=TestData.vars, type="prob")[,2]
                TestData[,"RPART"] <- trunc(1000*TestData[,"RPART"])
                TestPres <- TestData[TestData[,"pb"]==1,"RPART"]
                TestAbs <- TestData[TestData[,"pb"]==0,"RPART"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("RPART evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    evaluations$evaluations["RPART",2] <-  eval2@auc
                }
            }
            if (models.keep==T) {
                models$RPART <- results
                models$formulae$RPART.formula <- RPART.formula
            }
            if (evaluation.strip == T) {
                ResponseData$RPART <- predict(results, newdata=ResponseData.vars, type="prob")[,2]
            }
            fullname <- paste("models/", RASTER.species.name, "_RPART", sep="")
            tryCatch(prpart <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, type="prob", index=2,
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("RPART prediction failed"))},
                silent=T)
            if (is.null(prpart) == F) {
                prpart <- trunc(1000*prpart)
                writeRaster(x=prpart, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            }else{
                cat(paste("\n", "WARNING: RPART prediction failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["RPART"] <- -1 
            } 
        }
#       detach(package:rpart)
    }
    if (!is.null(factors) && NNET > 0) {
        cat(paste("\n", "ANN model with factors may require explicit dummy variables", "\n", sep=""))
    } 
    if (ws["NNET"] > 0) {
#        require(nnet, quietly=T)
        eval1 <- eval2 <- results <- pnnet <- TestPres <- TestAbs <- TrainPres <- TrainAbs <- NULL
        if(is.null(NNET.OLD) == T) {
            tryCatch(results <- nnet(formula=NNET.formula, size=NNET.size, decay=NNET.decay, data=TrainData, weights=Yweights, 
                    rang=0.1, maxit=maxit, trace=F),
                error= function(err) {print(paste("ANN calibration (nnet package) failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["NNET"] <- -1
            }
        }else{ results <- NNET.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with calibration data for ANN (nnet package)", "\n\n", sep = ""))
            TrainData[,"NNET"] <- predict.nnet2(object=results, newdata=TrainData.vars)
            TrainData[,"NNET"] <- trunc(1000*TrainData[,"NNET"])
            TrainPres <- TrainData[TrainData[,"pb"]==1,"NNET"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"NNET"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["NNET"] <- threshold(eval1, threshold.method) 
            evaluations$evaluations["NNET",1] <-  eval1@auc
            if (no.tests == F) {
                cat(paste("\n", "ANN evaluation with test data","\n\n", sep = ""))
                TestData[,"NNET"] <- predict.nnet2(object=results, newdata=TestData.vars)
                TestData[,"NNET"] <- trunc(1000*TestData[,"NNET"])
                TestPres <- TestData[TestData[,"pb"]==1,"NNET"]
                TestAbs <- TestData[TestData[,"pb"]==0,"NNET"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("ANN evaluation (nnet package) failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                    evaluations$evaluations["NNET",2] <-  eval2@auc
                }
            }
            if (models.keep==T) {
                models$NNET <- results
                models$formulae$NNET.formula <- NNET.formula
            }
            if (evaluation.strip == T) {
                ResponseData$NNET <- predict.nnet2(results, newdata=ResponseData.vars)
            }
            fullname <- paste("models/", RASTER.species.name, "_NNET", sep="")
            tryCatch(pnnet <- raster::predict(object=xn, model=results, fun=predict.nnet2, ext=ext, na.rm=TRUE, 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("ANN prediction (nnet package) failed"))},
                silent=T)
            if (is.null(pnnet) == F) {
                pnnet <- trunc(1000*pnnet)
                writeRaster(x=pnnet, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            }else{
                cat(paste("\n", "WARNING: ANN prediction (nnet package) failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["NNET"] <- -1 
            } 
        }
#       detach(package:nnet)
    }
    if (ws["FDA"] > 0) {
#        require(mda, quietly=T)
        eval1 <- eval2 <- results <- pfda <- TestPres <- TestAbs <- TrainPres <- TrainAbs <- NULL
        if(is.null(FDA.OLD) == T) {
            tryCatch(results <- fda(formula=FDA.formula, method=mars, data=TrainData, weights=Yweights),
                error= function(err) {print(paste("FDA calibration failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["FDA"] <- -1
            }
        }else{ results <- FDA.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with calibration data for FDA", "\n\n", sep = "")) 
            TrainData[,"FDA"] <- predict(object=results, newdata=TrainData.vars, type="posterior")[,2]
            TrainData[,"FDA"] <- trunc(1000*TrainData[,"FDA"])
            TrainPres <- TrainData[TrainData[,"pb"]==1,"FDA"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"FDA"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["FDA"] <- threshold(eval1, threshold.method)
            evaluations$evaluations["FDA",1] <-  eval1@auc
            if (no.tests == F) {
                cat(paste("\n", "FDA evaluation with test data","\n\n", sep = ""))
                TestData[,"FDA"] <- predict(object=results, newdata=TestData.vars, type="posterior")[,2]
                TestData[,"FDA"] <- trunc(1000*TestData[,"FDA"])
                TestPres <- TestData[TestData[,"pb"]==1,"FDA"]
                TestAbs <- TestData[TestData[,"pb"]==0,"FDA"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("FDA evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                    evaluations$evaluations["FDA",2] <-  eval2@auc
                }
            }
            if (models.keep==T) {
                models$FDA <- results
                models$formulae$FDA.formula <- FDA.formula
            }
            if (evaluation.strip == T) {
                ResponseData$FDA <- predict(results, newdata=ResponseData.vars, type="posterior")[,2]
            }
            fullname <- paste("models/", RASTER.species.name, "_FDA", sep="")
            tryCatch(pfda <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, type="posterior", index=2,
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("FDA prediction failed"))},
                silent=T)
            if (is.null(pfda) == F) {
                pfda <- trunc(1000*pfda)
                writeRaster(x=pfda, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            }else{
                cat(paste("\n", "WARNING: FDA prediction failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["FDA"] <- -1 
            } 
        }
#       detach(package:mda)
    }
    if (!is.null(factors) && SVM > 0) {
        cat(paste("\n", "SVM model with factors may require explicit dummy variables", "\n", sep=""))
    } 
    if (ws["SVM"] > 0) {
        cat(paste("\n\n"))
#        require(kernlab, quietly=T)
        eval1 <- eval2 <- results <- psvm <- TestPres <- TestAbs <- TrainPres <- TrainAbs <- NULL
        if(is.null(SVM.OLD) == T) {
            tryCatch(results <- ksvm(SVM.formula, data=TrainData, type="C-svc", prob.model=T),
                error= function(err) {print(paste("SVM calibration failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["SVM"] <- -1
            }
        }else{ results <- SVM.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with calibration data for SVM (kernlab package)", "\n\n", sep = "")) 
            TrainData[,"SVM"] <- kernlab::predict(object=results, newdata=TrainData.vars, type="probabilities")[,2]
            TrainData[,"SVM"] <- trunc(1000*TrainData[,"SVM"])
            TrainPres <- TrainData[TrainData[,"pb"]==1,"SVM"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"SVM"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["SVM"] <- threshold(eval1, threshold.method) 
            evaluations$evaluations["SVM",1] <-  eval1@auc
            if (no.tests == F) {
                cat(paste("\n", "SVM evaluation with test data","\n\n", sep = ""))
                TestData[,"SVM"] <- kernlab::predict(object=results, newdata=TestData.vars, type="probabilities")[,2]
                TestData[,"SVM"] <- trunc(1000*TestData[,"SVM"])
                TestPres <- TestData[TestData[,"pb"]==1,"SVM"]
                TestAbs <- TestData[TestData[,"pb"]==0,"SVM"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("SVM evaluation (kernlab package) failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                    evaluations$evaluations["SVM",2] <-  eval2@auc
                }
            }
            if (models.keep==T) {
                models$SVM <- results
                models$formulae$SVM.formula <- SVM.formula
            }
            if (evaluation.strip == T) {
                ResponseData$SVM <- kernlab::predict(results, newdata=ResponseData.vars, type="probabilities")[,2]
            }
            fullname <- paste("models/", RASTER.species.name, "_SVM", sep="")
            predfun <- as.function(kernlab::predict)
            tryCatch(psvm <- raster::predict(object=xn, model=results, fun=predfun, ext=ext, na.rm=TRUE, type="probabilities", index=2,
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("SVM prediction (kernlab package) failed"))},
                silent=T)
            if (is.null(psvm) == F) {
                psvm <- trunc(1000*psvm)
                writeRaster(x=psvm, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            }else{
                cat(paste("\n", "WARNING: SVM prediction (kernlab package) failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["SVM"] <- -1 
            } 
        }
#       detach(package:kernlab)
    }
    if (!is.null(factors) && SVME > 0) {
        cat(paste("\n", "SVME model with factors may require explicit dummy variables", "\n", sep=""))
    }
    if (ws["SVME"] > 0) {
#        require(e1071, quietly=T)
        eval1 <- eval2 <- results <- psvme <- TestPres <- TestAbs <- TrainPres <- TrainAbs <- NULL
        if(is.null(SVME.OLD) == T) {
            tryCatch(results <- svm(SVME.formula, data=TrainData, type="C-classification", kernel="polynomial", degree=3, probability=TRUE),
                error= function(err) {print(paste("SVM calibration (e1071 package) failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["SVME"] <- -1
            }
        }else{ results <- SVME.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with calibration data for SVM (e1071 package)", "\n\n", sep = ""))
            TrainData[,"SVME"] <- predict.svme(model=results, newdata=TrainData.vars)
            TrainData[,"SVME"] <- trunc(1000*TrainData[,"SVME"])
            TrainPres <- TrainData[TrainData[,"pb"]==1,"SVME"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"SVME"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["SVME"] <- threshold(eval1, threshold.method) 
            evaluations$evaluations["SVME",1] <-  eval1@auc
            if (no.tests == F) {
                cat(paste("\n", "SVME evaluation with test data","\n\n", sep = ""))
                TestData[,"SVME"] <- predict.svme(model=results, newdata=TestData.vars)
                TestData[,"SVME"] <- trunc(1000*TestData[,"SVME"])
                TestPres <- TestData[TestData[,"pb"]==1,"SVME"]
                TestAbs <- TestData[TestData[,"pb"]==0,"SVME"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("SVM evaluation (e1071 package) failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                    evaluations$evaluations["SVME",2] <-  eval2@auc
                }
            }
            if (models.keep==T) {
                models$SVME <- results
                models$formulae$SVME.formula <- SVME.formula
            }
            if (evaluation.strip == T) {
                ResponseData$SVME <- predict.svme(results, newdata=ResponseData.vars)
            }
            fullname <- paste("models/", RASTER.species.name, "_SVME", sep="")
            tryCatch(psvme <- raster::predict(object=xn, model=results, fun=predict.svme, ext=ext, na.rm=TRUE,
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("SVM prediction (e1071 package) failed"))},
                warning= function(war) {print(paste("SVM prediction (e1071 package) failed"))},
                silent=T)
            if (is.null(psvme) == F) {
                psvme <- trunc(1000*psvme)
                writeRaster(x=psvme, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            }else{
                cat(paste("\n", "WARNING: SVM prediction (e1071 package) failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["SVME"] <- -1 
            }
        }
#       detach(package:e1071)
    }
    if (BIOCLIM > 0 || DOMAIN > 0 || MAHAL > 0) {
        if(is.null(factors)==F) {
            for (i in 1:length(factors)) {
                TrainData.vars <- TrainData.vars[, which(colnames(TrainData.vars) != factors[i])]
                TrainData.pres <- TrainData.pres[, which(colnames(TrainData.pres) != factors[i])]
                TestData.vars <- TestData.vars[, which(colnames(TestData.vars) != factors[i])]
                xn <- dropLayer(xn, which(names(xn) == factors[i]))               
            }
        }
    }
    if (ws["BIOCLIM"] > 0) {  
        eval1 <- eval2 <- results <- pbio <- TestPres <- TestAbs <- TrainPres <- TrainAbs <- NULL
        if(is.null(BIOCLIM.OLD) == T) {
            tryCatch(results <- bioclim(x=TrainData.pres),
                error= function(err) {print(paste("BIOCLIM calibration failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["BIOCLIM"] <- -1
            }
            }else{ results <- BIOCLIM.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with calibration data for BIOCLIM", "\n\n", sep = ""))
            TrainData[,"BIOCLIM"] <- dismo::predict(object=results, x=TrainData.vars)
            TrainData[,"BIOCLIM"] <- trunc(1000*TrainData[,"BIOCLIM"])
            TrainPres <- TrainData[TrainData[,"pb"]==1,"BIOCLIM"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"BIOCLIM"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["BIOCLIM"] <- threshold(eval1, threshold.method) 
            evaluations$evaluations["BIOCLIM",1] <-  eval1@auc
            if (no.tests == F) {
                cat(paste("\n", "BIOCLIM evaluation with test data","\n\n", sep = ""))
                TestData[,"BIOCLIM"] <- dismo::predict(object=results, x=TestData.vars)
                TestData[,"BIOCLIM"] <- trunc(1000*TestData[,"BIOCLIM"])
                TestPres <- TestData[TestData[,"pb"]==1,"BIOCLIM"]
                TestAbs <- TestData[TestData[,"pb"]==0,"BIOCLIM"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("BIOCLIM evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                    evaluations$evaluations["BIOCLIM",2] <-  eval2@auc
                }
            }
            if (models.keep==T) {models$BIOCLIM <- results}
            if (evaluation.strip == T) {
                ResponseData$BIOCLIM <- dismo::predict(results, ResponseData.vars)
            }
            fullname <- paste("models/", RASTER.species.name, "_BIOCLIM", sep="")
            tryCatch(pbio <- dismo::predict(object=results, x=xn, ext=ext, na.rm=TRUE, 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("BIOCLIM prediction failed"))},
                silent=T)
            if (is.null(pbio) == F) {
                pbio <- trunc(1000*pbio)
                writeRaster(x=pbio, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            }else{
                cat(paste("\n", "WARNING: BIOCLIM prediction failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["BIOCLIM"] <- -1 
            }
        }
    }
    if (ws["DOMAIN"] > 0) {
        eval1 <- eval2 <- results <- pdom <- TestPres <- TestAbs <- TrainPres <- TrainAbs <- NULL
        if(is.null(DOMAIN.OLD) == T) {
            tryCatch(results <- domain(x=TrainData.pres),
                error= function(err) {print(paste("DOMAIN calibration failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["DOMAIN"] <- -1
            }
        }else{ results <- DOMAIN.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with calibration data for DOMAIN", "\n\n", sep = ""))
            TrainData[,"DOMAIN"] <- dismo::predict(object=results, x=TrainData.vars)
            TrainData[,"DOMAIN"] <- trunc(1000*TrainData[,"DOMAIN"])
            TrainPres <- TrainData[TrainData[,"pb"]==1,"DOMAIN"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"DOMAIN"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["DOMAIN"] <- threshold(eval1, threshold.method) 
            evaluations$evaluations["DOMAIN",1] <-  eval1@auc
            if (no.tests == F) {
                cat(paste("\n", "DOMAIN evaluation with test data","\n\n", sep = ""))
                TestData[,"DOMAIN"] <- dismo::predict(object=results, x=TestData.vars)
                TestData[,"DOMAIN"] <- trunc(1000*TestData[,"DOMAIN"])
                TestPres <- TestData[TestData[,"pb"]==1,"DOMAIN"]
                TestAbs <- TestData[TestData[,"pb"]==0,"DOMAIN"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("DOMAIN evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                    evaluations$evaluations["DOMAIN",2] <-  eval2@auc
                }
            }
            if (models.keep==T) {models$DOMAIN <- results}
            if (evaluation.strip == T) {
                ResponseData$DOMAIN <- dismo::predict(results, ResponseData.vars)
            }
            fullname <- paste("models/", RASTER.species.name, "_DOMAIN", sep="")
            tryCatch(pdom <- dismo::predict(object=results, x=xn, ext=ext, na.rm=TRUE, 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("DOMAIN prediction failed"))},
                silent=T)
            if (is.null(pdom) == F) {
                pdom <- trunc(1000*pdom)
                writeRaster(x=pdom, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            }else{
                cat(paste("\n", "WARNING: DOMAIN prediction failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["DOMAIN"] <- -1 
            }
        }
    }
    if (ws["MAHAL"] > 0) {  
        eval1 <- eval2 <- results <- pmahal <- TestPres <- TestAbs <- TrainPres <- TrainAbs <- NULL
        if(is.null(MAHAL.OLD) == T) {
            tryCatch(results <- mahal(x=TrainData.pres),
                error= function(err) {print(paste("Mahalanobis calibration failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["MAHAL"] <- -1
            }
        }else{ results <- MAHAL.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with calibration data for Mahalanobis", "\n\n", sep = ""))
            TrainData[,"MAHAL"] <- predict.mahal(model=results, newdata=TrainData.vars, MAHAL.shape=MAHAL.shape)
            TrainData[,"MAHAL"] <- trunc(1000*TrainData[,"MAHAL"])
            TrainPres <- TrainData[TrainData[,"pb"]==1,"MAHAL"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"MAHAL"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["MAHAL"] <- threshold(eval1, threshold.method)
            evaluations$evaluations["MAHAL",1] <-  eval1@auc
            if (no.tests == F) {
                cat(paste("\n", "Mahalanobis evaluation with test data","\n\n", sep = ""))
                TestData[,"MAHAL"] <- predict.mahal(model=results, newdata=TestData.vars, MAHAL.shape=MAHAL.shape)
                TestData[,"MAHAL"] <- trunc(1000*TestData[,"MAHAL"])
                TestPres <- TestData[TestData[,"pb"]==1,"MAHAL"]
                TestAbs <- TestData[TestData[,"pb"]==0,"MAHAL"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("MAHAL evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                    evaluations$evaluations["MAHAL",2] <-  eval2@auc
                }
            }
            if (models.keep==T) {models$MAHAL <- results}
            if (evaluation.strip == T) {
                ResponseData$MAHAL <- predict.mahal(model=results, newdata=ResponseData.vars, MAHAL.shape=MAHAL.shape)
            }
            fullname <- paste("models/", RASTER.species.name, "_MAHAL", sep="")
# not possible to use the predict.mahal function as raster::predict automatically reverts to dismo::predict for 'DistModel' objects
            tryCatch(pmahal <- dismo::predict(object=results, x=xn, ext=ext, na.rm=TRUE,
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("Mahalanobis prediction failed"))},
                silent=T)
            if (is.null(pmahal) == F) {
                pmahal <- pmahal - 1 - MAHAL.shape
                pmahal <- abs(pmahal)
                pmahal <- MAHAL.shape / pmahal
                pmahal <- trunc(1000*pmahal)
                writeRaster(x=pmahal, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            }else{
                cat(paste("\n", "WARNING: Mahalanobis prediction failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["MAHAL"] <- -1 
            }
        }
    }
    if (prediction.failures == T) {
        cat(paste("\n", "WARNING: some predictions failed","\n", sep = ""))
        cat(paste("\n", "actual weights that were used were (-1 indicates failed predictions):","\n", sep = ""))
        print(ws)
        ws[which(ws==-1)] <- 0
    }
# tune ensemble
    if((length(ENSEMBLE.decay) > 1) || (length(ENSEMBLE.best) > 1) || (length(ENSEMBLE.min) > 1)) {
        strategy.results <- ensemble.strategy(TrainData=TrainData, TestData=TestData,
            ENSEMBLE.decay=ENSEMBLE.decay, ENSEMBLE.best=ENSEMBLE.best, ENSEMBLE.min=ENSEMBLE.min)
        ws <- strategy.results$weights
#        names(ws) <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX", 
#            "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL")
        cat(paste("\n", "Weights for ensemble forecasting", "\n\n", sep = ""))
        print(ws)
    }
    if(models.keep==T) {models$output.weights <- ws}
    TrainData[,"ENSEMBLE"] <- ws["MAXENT"]*TrainData[,"MAXENT"] + ws["GBM"]*TrainData[,"GBM"] +
        ws["GBMSTEP"]*TrainData[,"GBMSTEP"] + ws["RF"]*TrainData[,"RF"] + ws["GLM"]*TrainData[,"GLM"] +
        ws["GLMSTEP"]*TrainData[,"GLMSTEP"] + ws["GAM"]*TrainData[,"GAM"] + ws["GAMSTEP"]*TrainData[,"GAMSTEP"] +
        ws["MGCV"]*TrainData[,"MGCV"] + ws["MGCVFIX"]*TrainData[,"MGCVFIX"] + ws["EARTH"]*TrainData[,"EARTH"] +
        ws["RPART"]*TrainData[,"RPART"] + ws["NNET"]*TrainData[,"NNET"] + ws["FDA"]*TrainData[,"FDA"] +
        ws["SVM"]*TrainData[,"SVM"] + ws["SVME"]*TrainData[,"SVME"] + ws["BIOCLIM"]*TrainData[,"BIOCLIM"] +
        ws["DOMAIN"]*TrainData[,"DOMAIN"] + ws["MAHAL"]*TrainData[,"MAHAL"]
    TrainData[,"ENSEMBLE"] <- trunc(TrainData[,"ENSEMBLE"])  
    TestData[,"ENSEMBLE"] <- ws["MAXENT"]*TestData[,"MAXENT"] + ws["GBM"]*TestData[,"GBM"] +
        ws["GBMSTEP"]*TestData[,"GBMSTEP"] + ws["RF"]*TestData[,"RF"] + ws["GLM"]*TestData[,"GLM"] +
        ws["GLMSTEP"]*TestData[,"GLMSTEP"] + ws["GAM"]*TestData[,"GAM"] + ws["GAMSTEP"]*TestData[,"GAMSTEP"] +
        ws["MGCV"]*TestData[,"MGCV"] + ws["MGCVFIX"]*TestData[,"MGCVFIX"] + ws["EARTH"]*TestData[,"EARTH"] +
        ws["RPART"]*TestData[,"RPART"] + ws["NNET"]*TestData[,"NNET"] + ws["FDA"]*TestData[,"FDA"] +
        ws["SVM"]*TestData[,"SVM"] + ws["SVME"]*TestData[,"SVME"] + ws["BIOCLIM"]*TestData[,"BIOCLIM"] +
        ws["DOMAIN"]*TestData[,"DOMAIN"] + ws["MAHAL"]*TestData[,"MAHAL"]
    TestData[,"ENSEMBLE"] <- trunc(TestData[,"ENSEMBLE"])
#
    if (evaluation.strip == T) {
        ResponseData[,"ENSEMBLE"] <- ws["MAXENT"]*ResponseData[,"MAXENT"] + ws["GBM"]*ResponseData[,"GBM"] +
            ws["GBMSTEP"]*ResponseData[,"GBMSTEP"] + ws["RF"]*ResponseData[,"RF"] + ws["GLM"]*ResponseData[,"GLM"] +
            ws["GLMSTEP"]*ResponseData[,"GLMSTEP"] + ws["GAM"]*ResponseData[,"GAM"] + ws["GAMSTEP"]*ResponseData[,"GAMSTEP"] +
            ws["MGCV"]*ResponseData[,"MGCV"] + ws["MGCVFIX"]*ResponseData[,"MGCVFIX"] + ws["EARTH"]*ResponseData[,"EARTH"] +
            ws["RPART"]*ResponseData[,"RPART"] + ws["NNET"]*ResponseData[,"NNET"] + ws["FDA"]*ResponseData[,"FDA"] +
            ws["SVM"]*ResponseData[,"SVM"] + ws["SVME"]*ResponseData[,"SVME"] + ws["BIOCLIM"]*ResponseData[,"BIOCLIM"] +
            ws["DOMAIN"]*ResponseData[,"DOMAIN"] + ws["MAHAL"]*ResponseData[,"MAHAL"]
        evaluations$evaluation.strip <- ResponseData
    }
#
# create ensembles
    if (ENSEMBLE.rasters == T) {
        ensemble.statistics["n.models"] <- sum(as.numeric(ws > 0))
        ensemble <- xn[[1]] == NAvalue(xn[[1]])
        setMinMax(ensemble)
        if (is.null(ext) == F) {ensemble <- crop(ensemble, y=ext, snap="in")}
        writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        enscount <- ensemble
        setMinMax(enscount)
        writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
        enspresence <- ensemble
        setMinMax(enspresence)
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
        writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        setMinMax(enscount)
        ensemble.statistics["count.min"] <- minValue(enscount)
        ensemble.statistics["count.max"] <- maxValue(enscount)
    }
    eval1 <- eval2 <- TestPres <- TestAbs <- TrainPres <- TrainAbs <- NULL
    cat(paste("\n", RASTER.species.orig, ": Ensemble evaluation with calibration data", "\n\n", sep = ""))
    TrainPres <- TrainData[TrainData[,"pb"]==1,"ENSEMBLE"]
    TrainAbs <- TrainData[TrainData[,"pb"]==0,"ENSEMBLE"]
    eval1 <- evaluate(p=TrainPres, a=TrainAbs)
    print(eval1)
    evaluations$evaluations["ENSEMBLE",1] <-  eval1@auc
    if (ENSEMBLE.rasters == T) {
        l3 <- threshold(eval1, threshold.method)
        ensemble.statistics["ensemble.threshold"] <- l3
        cat(paste("\n", "Suggested threshold for presence (calculated as: ", threshold.method, ")", "\n", sep = ""))    
        cat(paste(l3, "\n", sep=""))
        enspresence <- ensemble > l3
        writeRaster(x=enspresence, filename=rasterpresence, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
        if(is.null(a) == F) {
            cat(paste("\n", "Evaluation of created ensemble raster layer (", rasterfull, ") at locations p and a", "\n\n", sep = ""))
            pres_consensus <- extract(ensemble, p)
            abs_consensus <- extract(ensemble, a)
            evalc <- evaluate(p=pres_consensus, a=abs_consensus)
            print(evalc)
        }
    }
    if (no.tests == F) {
        cat(paste("\n", RASTER.species.orig, ": Ensemble evaluation with testing data", "\n\n", sep = ""))
        TestPres <- TestData[TestData[,"pb"]==1,"ENSEMBLE"]
        TestAbs <- TestData[TestData[,"pb"]==0,"ENSEMBLE"]
        eval2 <- evaluate(p=TestPres, a=TestAbs)
        print(eval2)
        evaluations$evaluations["ENSEMBLE",2] <-  eval2@auc
        if(is.null(pt) == F) {
            cat(paste("\n", "Evaluation of created ensemble raster layer (", rasterfull, ") at locations pt and at", "\n\n", sep = ""))
            pres_consensus <- extract(ensemble, pt)
            abs_consensus <- extract(ensemble, at)
            evalc <- evaluate(p=pres_consensus, a=abs_consensus)
            print(evalc)
        }
    }
    cat(paste("\n", "Summary of evaluations (AUC as percentage)", "\n", sep = ""))
    if (no.tests == F) {
        evaluations$evaluations <- evaluations$evaluations[order(evaluations$evaluations[,2], decreasing=T),]
    }else{
        evaluations$evaluations <- evaluations$evaluations[order(evaluations$evaluations[,1], decreasing=T),]
    }
    print(evaluations$evaluations)
    cat(paste("\n", "End of modelling for organism: ", RASTER.species.orig, "\n\n", sep = ""))      
    remove(Yweights, envir=.BiodiversityR)
    remove(TrainData, envir=.BiodiversityR)
    remove(TestData,envir=.BiodiversityR)
    if (evaluation.strip==T) {remove(ResponseData, envir=.BiodiversityR)}
    return(list(evaluations=evaluations, models=models, ensemble.statistics=ensemble.statistics, call=match.call() ))
}




