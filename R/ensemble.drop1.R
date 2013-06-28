`ensemble.drop1` <- function(
    x=NULL, p=NULL, a=NULL, an=1000, excludep=FALSE, ext=NULL, k=0, pt=NULL, at=NULL,
    TrainData=NULL, TestData=NULL,
    layer.drops=NULL, VIF=FALSE, COR=FALSE,
    difference=FALSE,
    ENSEMBLE.decay=1, ENSEMBLE.best=0, ENSEMBLE.min=0.7,
    input.weights=NULL,
    MAXENT=1, GBM=1, GBMSTEP=1, RF=1, GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=1, MGCV=1, MGCVFIX=0, 
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, BIOCLIM=1, DOMAIN=1, MAHAL=1, 
    Yweights="BIOMOD", factors=NULL, dummy.vars=NULL,
    maxit=100,
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
    MAHAL.shape=1
)
{
    .BiodiversityR <- new.env()
    if (! require(dismo)) {stop("Please install the dismo package")}
    k <- as.integer(k)
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
# 
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
        GEODIST <- max(c(input.weights["GEODIST"], -1), na.rm=T)
    }
# no point for GEODIST
    GEODIST <- -1
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
        if (identical(RF.ntree, trunc(RF.ntree/2)) == F) {RF.ntree <- RF.ntree + 1}
    }
    if (GLMSTEP > 0) {
        if (! require(MASS)) {stop("Please install the MASS package")}
    }
    if (GAM > 0) {
        if (! require(gam)) {stop("Please install the gam package")}
        detach(package:gam)
    }
    if (GAMSTEP > 0) {
        if (! require(gam)) {stop("Please install the gam package")}
        detach(package:gam)   
    }
    if (MGCV > 0) {
        if (! require(mgcv)) {stop("Please install the mgcv package")}
        detach(package:mgcv)   
#         get the probabilities from MGCV
            predict.mgcv <- function(object, newdata, type="response") {
                p <- predict(object=object, newdata=newdata, type=type)
                return(as.numeric(p))
            }
    }
    if (MGCVFIX > 0) {
        if (! require(mgcv)) {stop("Please install the mgcv package")}
        detach(package:mgcv)   
#         get the probabilities from MGCV
            predict.mgcv <- function(object, newdata, type="response") {
                p <- predict(object=object, newdata=newdata, type=type)
                return(as.numeric(p))
            }
    }
    if (EARTH > 0) {
        if (! require(earth)) {stop("Please install the earth package")}
    }
    if (RPART > 0) {
        if (! require(rpart)) {stop("Please install the rpart package")}
    }
    if (NNET > 0) {
        if (! require(nnet)) {stop("Please install the nnet package")}
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
        if(is.null(TestData) == T) {
            TestData <- TrainData
            if (k > 1) {
                groupp <- kfold(TrainData, k=k, by=TrainData[,"pb"])
                TrainData.c <- TrainData[groupp != 1,]
                TestData <- TrainData[groupp == 1,]
                TrainData <- TrainData.c
            }
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
    TrainData.orig <- TrainData
    assign("TrainData.orig", TrainData.orig, envir=.BiodiversityR)
    TestData.orig <- TestData
    assign("TestData.orig", TestData.orig, envir=.BiodiversityR)
    vars <- colnames(TrainData.orig)
    vars <- vars[which(vars != "pb")]
    nv <- length(vars)
#
# create background data for MAXENT
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
        MAXENT.BackData.orig <- data.frame(MAXENT.BackData)
        if (is.null(layer.drops) == F) {
            nd <- length(layer.drops)
            for (i in 1:nd) {     
                MAXENT.BackData.orig <- MAXENT.BackData.orig[, which(colnames(MAXENT.BackData.orig) != layer.drops[i])]
            }
        }
        TestValid <- complete.cases(MAXENT.BackData.orig)
        MAXENT.BackData.orig <- MAXENT.BackData.orig[TestValid,]  
        assign("MAXENT.BackData.orig", MAXENT.BackData.orig, envir=.BiodiversityR)
    }
# 
    output <- array(NA, dim=c(20, nv+1))
    rownames(output) <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL", "ENSEMBLE")
    colnames(output) <- c("all_vars", paste("without_", vars, sep=""))
# first fit with all variables
    tests <- ensemble.test(x=x,
        TrainData=TrainData.orig, TestData=TestData.orig,
        PLOTS=FALSE, evaluations.keep=T, models.keep=F,
        layer.drops=NULL, VIF=VIF, COR=COR,
        formulae.defaults=T, maxit=maxit,
        ENSEMBLE.decay=ENSEMBLE.decay, ENSEMBLE.best=ENSEMBLE.best, ENSEMBLE.min=ENSEMBLE.min,
        MAXENT=MAXENT, GBM=GBM, GBMSTEP=GBMSTEP, RF=RF, GLM=GLM, GLMSTEP=GLMSTEP, 
        GAM=GAM, GAMSTEP=GAMSTEP, MGCV=MGCV, MGCVFIX=MGCVFIX, EARTH=EARTH, RPART=RPART, 
        NNET=NNET, FDA=FDA, SVM=SVM, SVME=SVME, BIOCLIM=BIOCLIM, DOMAIN=DOMAIN, MAHAL=MAHAL,
        GEODIST=0, 
        Yweights=Yweights, factors=factors, dummy.vars=dummy.vars,
        MAXENT.BackData=MAXENT.BackData.orig, MAXENT.path=MAXENT.path,
        GBM.formula=NULL, GBM.n.trees=GBM.n.trees,
        GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
        GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction,
        GBMSTEP.step.size=GBMSTEP.step.size,
        RF.formula=NULL, RF.ntree=RF.ntree, RF.mtry=RF.mtry, 
        GLM.formula=NULL, GLM.family=GLM.family, 
        GLMSTEP.k=GLMSTEP.k, GLMSTEP.steps=GLMSTEP.steps, STEP.formula=NULL, GLMSTEP.scope=NULL, 
        GAM.formula=NULL, GAM.family=GAM.family, 
        GAMSTEP.steps=GAMSTEP.steps, GAMSTEP.scope=NULL, GAMSTEP.pos=GAMSTEP.pos,
        MGCV.formula=NULL, MGCV.select=MGCV.select,
        MGCVFIX.formula=NULL, 
        EARTH.formula=NULL, EARTH.glm=EARTH.glm,
        RPART.formula=NULL, RPART.xval=RPART.xval, 
        NNET.formula=NULL, NNET.size=NNET.size, NNET.decay=NNET.decay,
        FDA.formula=NULL, SVM.formula=NULL, SVME.formula=NULL,
        MAHAL.shape=MAHAL.shape)
    if(is.null(tests$evaluations$MAXENT.T)==F) {output["MAXENT",1] <- tests$evaluations$MAXENT.T@auc}
    if(is.null(tests$evaluations$GBM.T)==F) {output["GBM",1] <- tests$evaluations$GBM.T@auc} 
    if(is.null(tests$evaluations$GBMSTEP.T)==F) {output["GBMSTEP",1] <- tests$evaluations$GBMSTEP.T@auc} 
    if(is.null(tests$evaluations$RF.T)==F) {output["RF",1] <- tests$evaluations$RF.T@auc}
    if(is.null(tests$evaluations$GLM.T)==F) {output["GLM",1] <- tests$evaluations$GLM.T@auc} 
    if(is.null(tests$evaluations$GLMS.T)==F) {output["GLMSTEP",1] <- tests$evaluations$GLMS.T@auc}
    if(is.null(tests$evaluations$GAM.T)==F) {output["GAM",1] <- tests$evaluations$GAM.T@auc} 
    if(is.null(tests$evaluations$GAMS.T)==F) {output["GAMSTEP",1] <- tests$evaluations$GAMS.T@auc}
    if(is.null(tests$evaluations$MGCV.T)==F) {output["MGCV",1] <- tests$evaluations$MGCV.T@auc} 
    if(is.null(tests$evaluations$MGCVF.T)==F) {output["MGCVFIX",1] <- tests$evaluations$MGCVF.T@auc} 
    if(is.null(tests$evaluations$EARTH.T)==F) {output["EARTH",1] <- tests$evaluations$EARTH.T@auc} 
    if(is.null(tests$evaluations$RPART.T)==F) {output["RPART",1] <- tests$evaluations$RPART.T@auc}
    if(is.null(tests$evaluations$NNET.T)==F) {output["NNET",1] <- tests$evaluations$NNET.T@auc} 
    if(is.null(tests$evaluations$FDA.T)==F) {output["FDA",1] <- tests$evaluations$FDA.T@auc}
    if(is.null(tests$evaluations$SVM.T)==F) {output["SVM",1] <- tests$evaluations$SVM.T@auc}
    if(is.null(tests$evaluations$SVME.T)==F) {output["SVME",1] <- tests$evaluations$SVME.T@auc}
    if(is.null(tests$evaluations$BIOCLIM.T)==F) {output["BIOCLIM",1] <- tests$evaluations$BIOCLIM.T@auc}
    if(is.null(tests$evaluations$DOMAIN.T)==F) {output["DOMAIN",1] <- tests$evaluations$DOMAIN.T@auc}
    if(is.null(tests$evaluations$MAHAL.T)==F) {output["MAHAL",1] <- tests$evaluations$MAHAL.T@auc}
    if(is.null(tests$evaluations$ENSEMBLE.T)==F) {output["ENSEMBLE",1] <- tests$evaluations$ENSEMBLE.T@auc}

# sequentially leave out the focal variable, then fit again

    vars <- colnames(TrainData.orig)
    vars <- vars[which(vars != "pb")]
    nv <- length(vars)
    for (i in 1:nv) {
        var.f <- vars[i]
        cat(paste("\n", "Leaving out variable ", var.f, "\n\n", sep = ""))
    tests <- ensemble.test(x=x, 
        TrainData=TrainData.orig, TestData=TestData.orig, 
        PLOTS=FALSE, evaluations.keep=T, models.keep=F,
        layer.drops=var.f, VIF=VIF, COR=COR,
        formulae.defaults=T, maxit=maxit,
        ENSEMBLE.decay=ENSEMBLE.decay, ENSEMBLE.best=ENSEMBLE.best, ENSEMBLE.min=ENSEMBLE.min,
        MAXENT=MAXENT, GBM=GBM, GBMSTEP=GBMSTEP, RF=RF, GLM=GLM, GLMSTEP=GLMSTEP, 
        GAM=GAM, GAMSTEP=GAMSTEP, MGCV=MGCV, MGCVFIX=MGCVFIX, EARTH=EARTH, RPART=RPART, 
        NNET=NNET, FDA=FDA, SVM=SVM, SVME=SVME, BIOCLIM=BIOCLIM, DOMAIN=DOMAIN, MAHAL=MAHAL,
        GEODIST=0,
        Yweights=Yweights, factors=factors, dummy.vars=dummy.vars,
        MAXENT.BackData=MAXENT.BackData.orig, MAXENT.path=MAXENT.path,
        GBM.formula=NULL, GBM.n.trees=GBM.n.trees,
        GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
        GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction,
        GBMSTEP.step.size=GBMSTEP.step.size,
        RF.formula=NULL, RF.ntree=RF.ntree, RF.mtry=RF.mtry, 
        GLM.formula=NULL, GLM.family=GLM.family, 
        GLMSTEP.k=GLMSTEP.k, GLMSTEP.steps=GLMSTEP.steps, STEP.formula=NULL, GLMSTEP.scope=NULL, 
        GAM.formula=NULL, GAM.family=GAM.family, 
        GAMSTEP.steps=GAMSTEP.steps, GAMSTEP.scope=NULL, GAMSTEP.pos=GAMSTEP.pos,
        MGCV.formula=NULL, MGCV.select=MGCV.select,
        MGCVFIX.formula=NULL, 
        EARTH.formula=NULL, EARTH.glm=EARTH.glm,
        RPART.formula=NULL, RPART.xval=RPART.xval, 
        NNET.formula=NULL, NNET.size=NNET.size, NNET.decay=NNET.decay,
        FDA.formula=NULL, SVM.formula=NULL, SVME.formula=NULL,
        MAHAL.shape=MAHAL.shape)
    if(is.null(tests$evaluations$MAXENT.T)==F) {output["MAXENT",i+1] <- tests$evaluations$MAXENT.T@auc}
    if(is.null(tests$evaluations$GBM.T)==F) {output["GBM",i+1] <- tests$evaluations$GBM.T@auc} 
    if(is.null(tests$evaluations$GBMSTEP.T)==F) {output["GBMSTEP",i+1] <- tests$evaluations$GBMSTEP.T@auc} 
    if(is.null(tests$evaluations$RF.T)==F) {output["RF",i+1] <- tests$evaluations$RF.T@auc}
    if(is.null(tests$evaluations$GLM.T)==F) {output["GLM",i+1] <- tests$evaluations$GLM.T@auc} 
    if(is.null(tests$evaluations$GLMS.T)==F) {output["GLMSTEP",i+1] <- tests$evaluations$GLMS.T@auc}
    if(is.null(tests$evaluations$GAM.T)==F) {output["GAM",i+1] <- tests$evaluations$GAM.T@auc} 
    if(is.null(tests$evaluations$GAMS.T)==F) {output["GAMSTEP",i+1] <- tests$evaluations$GAMS.T@auc}
    if(is.null(tests$evaluations$MGCV.T)==F) {output["MGCV",i+1] <- tests$evaluations$MGCV.T@auc} 
    if(is.null(tests$evaluations$MGCVF.T)==F) {output["MGCVFIX",i+1] <- tests$evaluations$MGCVF.T@auc} 
    if(is.null(tests$evaluations$EARTH.T)==F) {output["EARTH",i+1] <- tests$evaluations$EARTH.T@auc} 
    if(is.null(tests$evaluations$RPART.T)==F) {output["RPART",i+1] <- tests$evaluations$RPART.T@auc}
    if(is.null(tests$evaluations$NNET.T)==F) {output["NNET",i+1] <- tests$evaluations$NNET.T@auc} 
    if(is.null(tests$evaluations$FDA.T)==F) {output["FDA",i+1] <- tests$evaluations$FDA.T@auc}
    if(is.null(tests$evaluations$SVM.T)==F) {output["SVM",i+1] <- tests$evaluations$SVM.T@auc}
    if(is.null(tests$evaluations$SVME.T)==F) {output["SVME",i+1] <- tests$evaluations$SVME.T@auc}
    if(is.null(tests$evaluations$BIOCLIM.T)==F) {output["BIOCLIM",i+1] <- tests$evaluations$BIOCLIM.T@auc}
    if(is.null(tests$evaluations$DOMAIN.T)==F) {output["DOMAIN",i+1] <- tests$evaluations$DOMAIN.T@auc}
    if(is.null(tests$evaluations$MAHAL.T)==F) {output["MAHAL",i+1] <- tests$evaluations$MAHAL.T@auc}
    if(is.null(tests$evaluations$ENSEMBLE.T)==F) {output["ENSEMBLE",i+1] <- tests$evaluations$ENSEMBLE.T@auc}

    }
    output <- 100*output
    if (difference == T) {
        for (i in 1:nv) {
            output[,i+1] <- output[,i+1] - output[,1]
        }
    }
    output <- output[order(output[,"all_vars"], decreasing=T),]
    cat(paste("\n", "Results (AUC as percentage)",  "\n\n", sep = ""))
    if (difference == T) {
        cat(paste("\n", "Note that positive differences indicate that the model without the variable", sep = ""))
        cat(paste("\n", "has higher AUC than the model with all the variables",  "\n\n", sep = ""))
    }
    print (output)
    cat(paste("\n\n"))
    return(list(table=output, call=match.call() ))
}


