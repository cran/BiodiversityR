`ensemble.test` <- function(
    x=NULL, p=NULL, a=NULL, an=1000, excludep=FALSE, ext=NULL, k=0, pt=NULL, at=NULL,
    TrainData=NULL, TestData=NULL,
    layer.drops=NULL, VIF=FALSE, COR=FALSE,
    PLOTS=TRUE, evaluations.keep=FALSE, models.keep=FALSE, 
    ENSEMBLE.decay=1, ENSEMBLE.best=0, ENSEMBLE.min=0.7,
    input.weights=NULL, models.list=NULL,
    MAXENT=1, GBM=1, GBMSTEP=1, RF=1, GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=1, MGCV=1, MGCVFIX=0,
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, BIOCLIM=1, DOMAIN=1, MAHAL=1, 
    GEODIST=0, 
    Yweights="BIOMOD", factors=NULL, dummy.vars=NULL,
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
    MAHAL.shape=1, MAHAL.OLD=NULL,
    GEODIST.file.name="Species001", RASTER.format="raster", GEODIST.OLD=NULL
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
    if (is.null(input.weights) == F) {
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
#
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
# check formulae and packages
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
# check if TestData is different from TrainData
    no.tests <- FALSE
    if (identical(TrainData, TestData) == T) {no.tests <- TRUE}
#
    if(is.null(factors)==F) {
        for (i in 1:length(factors)) {
            if (identical(levels(TrainData[,factors[i]]),levels(TestData[,factors[i]]))==F) {
                cat(paste("\n", "WARNING: differences in factor levels between calibration and evaluation data (variable ", factors[i], ")", "\n", sep = ""))
                cat(paste("same levels set for both data sets to avoid problems with evaluations for RF and SVM","\n", sep = ""))
                uniquelevels <- unique(c(levels(TrainData[,factors[i]]), levels(TestData[,factors[i]])))
                levels(TrainData[,factors[i]]) <- uniquelevels
                levels(TestData[,factors[i]]) <- uniquelevels
            } 
        }
    }
#
    if(models.keep==T) {
        models <- list(MAXENT=NULL, GBM=NULL, GBMSTEP=NULL, RF=NULL, GLM=NULL, 
            GLMSTEP=NULL, GAM=NULL, GAMSTEP=NULL, MGCV=NULL, MGCVFIX=NULL, EARTH=NULL, RPART=NULL, 
            NNET=NULL, FDA=NULL, SVM=NULL, SVME=NULL, BIOCLIM=NULL, DOMAIN=NULL, MAHAL=NULL, GEODIST=NULL,
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
    if (VIF == T) {
        if (! require(car)) {stop("Please install the car package")}
        GLM2.formula <- ensemble.formulae(TrainData, factors=factors)$RF.formula
        assign("GLM2.formula", GLM2.formula, envir=.BiodiversityR)
        GLM2.family=binomial(link="logit")
        assign("GLM2.family", GLM2.family, envir=.BiodiversityR)
        vifresult <- NULL
        tryCatch(vifresult <- vif(glm(formula=GLM2.formula, family=GLM2.family, data=TrainData, control=glm.control(maxit=maxit))),
            error= function(err) {print(paste("VIF evaluation failed"))},
                    silent=F)
        if (is.null(vifresult) == F) {
            cat(paste("\n", "Variance inflation (glm)", "\n", sep = ""))        
            print(vifresult)
        }else{
            cat(paste("\n", "NOTE: VIF evaluation failed", "\n", sep = "")) 
        }
    }
    if (COR == T) {
        TrainDataNum <- TrainData[, colnames(TrainData) != "pb"]
        if(is.null(factors)==F) {
            for (i in 1:length(factors)) {
                TrainDataNum <- TrainDataNum[, colnames(TrainDataNum) != factors[i]]            
            }
        }
        corresult <- cor(TrainDataNum)
        corresult <- round(100*corresult, digits=2)
        cat(paste("\n", "Correlation between numeric variables (as percentage)", "\n", sep = ""))        
        print(corresult)
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
        if (models.keep == T) {models$MAXENT.BackData <- MAXENT.BackData}
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
    weights <- as.numeric(array(dim=19, 0))
    names(weights) <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX", 
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL")
#
    if(evaluations.keep==T) {
        evaluations <- list(MAXENT.C=NULL, MAXENT.T=NULL, 
            GBM.trees=NULL, GBM.C=NULL, GBM.T=NULL, GBMSTEP.trees=NULL, GBMSTEP.C=NULL, GBMSTEP.T=NULL, 
            RF.C=NULL, RF.T=NULL, GLM.C=NULL, GLM.T=NULL, GLMS.C=NULL, GLMS.T=NULL, 
            GAM.C=NULL, GAM.T=NULL, GAMS.C=NULL, GAMS.T=NULL, MGCV.C=NULL, MGCV.T=NULL, MGCVF.C=NULL, MGCVF.T=NULL,
            EARTH.C=NULL, EARTH.T=NULL, RPART.C=NULL, RPART.T=NULL,
            NNET.C=NULL, NNET.T=NULL, FDA.C=NULL, FDA.T=NULL, SVM.C=NULL, SVM.T=NULL, SVME.C=NULL, SVME.T=NULL,
            BIOCLIM.C=NULL, BIOCLIM.T=NULL, DOMAIN.C=NULL, DOMAIN.T=NULL, MAHAL.C=NULL, MAHAL.T=NULL,
            GEODIST.C=NULL, GEODIST.T=NULL, ENSEMBLE.C=NULL, ENSEMBLE.T=NULL, STRATEGY.weights=NULL,
            TrainData=NULL, TestData=NULL)
    }else{
        evaluations <- NULL
    }
#
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
    cat(paste("\n", "Start of evaluations", "\n\n",sep = ""))
# Different modelling algorithms
    if (MAXENT > 0) {
        eval1 <- eval2 <- results <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        # Put the file 'maxent.jar' in the 'java' folder of dismo
        # the file 'maxent.jar' can be obtained from from http://www.cs.princeton.edu/~schapire/maxent/.
        jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
        if(is.null(MAXENT.OLD) == T) {
            tryCatch(results <- maxent(x=MAXENT.TrainData, p=MAXENT.pa, factors=factors, path=MAXENT.path),
                error= function(err) {print(paste("MAXENT calibration failed"))},
                silent=T)
        }else{ 
            results <- MAXENT.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "MAXENT calibration (tested with calibration data of other algorithms)","\n\n",sep = ""))
            TrainData[,"MAXENT"] <- dismo::predict(object=results, x=TrainData.vars)
            TrainPres <- TrainData[TrainData[,"pb"]==1,"MAXENT"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"MAXENT"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            weights["MAXENT"] <- max(c(eval1@auc, 0), na.rm=T)
#            TrainData[,"MAXENT"] <- trunc(1000*TrainData[,"MAXENT"])
            if (no.tests == F) {
                cat(paste("\n", "MAXENT evaluation","\n\n",sep = ""))
                TestData[,"MAXENT"] <- dismo::predict(object=results, x=TestData.vars)
                TestPres <- TestData[TestData[,"pb"]==1,"MAXENT"]
                TestAbs <- TestData[TestData[,"pb"]==0,"MAXENT"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("MAXENT evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["MAXENT"] <- max(c(eval2@auc, 0), na.rm=T)
#                    TestData[,"MAXENT"] <- trunc(1000*TestData[,"MAXENT"])
                    if(PLOTS==T) {
                        plot(eval2, "ROC") 
                        title(main="MAXENT", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: MAXENT evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$MAXENT.C <- eval1
                evaluations$MAXENT.T <- eval2
            }
            if(models.keep==T) {
                models$MAXENT <- results
            }
        }else{ cat(paste("\n", "WARNING: MAXENT calibration failed", "\n", "\n"))}
    }
    if (GBM > 0) {
#        require(gbm, quietly=T)        
        eval1 <- eval2 <- results <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(GBM.OLD) == T) {       
            tryCatch(results <- gbm(formula=GBM.formula, data=TrainData, weights=Yweights, distribution="bernoulli", 
                    interaction.depth=7, shrinkage=0.001, bag.fraction=0.5, train.fraction=1, 
                    n.trees=GBM.n.trees, verbose=F, cv.folds=5),
                error= function(err) {print(paste("GBM calibration failed"))},
                silent=T)
        }else{ 
            results <- GBM.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "GBM calibration","\n\n",sep = ""))
            TrainData[,"GBM"] <- predict(object=results, newdata=TrainData.vars, n.trees=results$n.trees, type="response")
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GBM"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GBM"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            weights["GBM"] <- max(c(eval1@auc, 0), na.rm=T)
#            TrainData[,"GBM"] <- trunc(1000*TrainData[,"GBM"])
            if (no.tests == F) {
                cat(paste("\n", "GBM evaluation","\n\n",sep = ""))
                TestData[,"GBM"] <- predict(object=results, newdata=TestData.vars, n.trees=results$n.trees, type="response")
                TestPres <- TestData[TestData[,"pb"]==1,"GBM"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GBM"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("GBM evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["GBM"] <- max(c(eval2@auc, 0), na.rm=T)
#                    TestData[,"GBM"] <- trunc(1000*TestData[,"GBM"])
                    if(PLOTS==T) {
                        plot(eval2, "ROC") 
                        title(main="BRT (gbm)", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: GBM evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$GBM.trees <- results$n.trees 
                evaluations$GBM.C <- eval1
                evaluations$GBM.T <- eval2
            }
            if(models.keep==T) {
                models$GBM <- results
                models$formulae$GBM.formula <- GBM.formula
            }
        }else{ cat(paste("\n", "WARNING: GBM calibration failed", "\n", "\n"))}
#       detach(package:gbm) 
    }
    if (GBMSTEP > 0) {
#        require(gbm, quietly=T)        
        eval1 <- eval2 <- results <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(GBMSTEP.OLD) == T) {       
            tryCatch(results <- gbm.step(data=TrainData, gbm.y=1, gbm.x=GBMSTEP.gbm.x, family="bernoulli",
                    site.weights=Yweights, tree.complexity=GBMSTEP.tree.complexity, learning.rate = GBMSTEP.learning.rate, 
                    bag.fraction=GBMSTEP.bag.fraction, step.size=GBMSTEP.step.size, verbose=F, silent=T, plot.main=F),
                error= function(err) {print(paste("stepwise GBM calibration failed"))},
                silent=T)
        }else{ 
            results <- GBMSTEP.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "stepwise GBM trees (target > 1000)", "\n", sep = ""))
            print(results$n.trees)
            cat(paste("\n", "stepwise GBM calibration","\n\n",sep = ""))
            TrainData[,"GBMSTEP"] <- predict(object=results, newdata=TrainData.vars, n.trees=results$n.trees, type="response")
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GBMSTEP"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GBMSTEP"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            weights["GBMSTEP"] <- max(c(eval1@auc, 0), na.rm=T)
#            TrainData[,"GBMSTEP"] <- trunc(1000*TrainData[,"GBMSTEP"])
            if (no.tests == F) {
                cat(paste("\n", "stepwise GBM evaluation","\n\n",sep = ""))
                TestData[,"GBMSTEP"] <- predict(object=results, newdata=TestData.vars, n.trees=results$n.trees, type="response")
                TestPres <- TestData[TestData[,"pb"]==1,"GBMSTEP"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GBMSTEP"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("stepwise GBM evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["GBMSTEP"] <- max(c(eval2@auc, 0), na.rm=T)
#                    TestData[,"GBMSTEP"] <- trunc(1000*TestData[,"GBMSTEP"])
                    if(PLOTS==T) {
                        plot(eval2, "ROC") 
                        title(main="BRT (gbm.step)", cex=2, adj=0, col.main="blue")
                   }
                }else{
                    cat(paste("\n", "WARNING: stepwise GBM evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$GBMSTEP.trees <- results$n.trees 
                evaluations$GBMSTEP.C <- eval1
                evaluations$GBMSTEP.T <- eval2
            }
            if(models.keep==T) {
                models$GBMSTEP <- results
            }
        }else{ cat(paste("\n", "WARNING: stepwise GBM calibration failed", "\n", "\n"))}
#       detach(package:gbm) 
    }
    if (RF > 0) {
#        require(randomForest, quietly=T)
        eval1 <- eval2 <- results <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(RF.OLD) == T) {        
            tryCatch(results <- randomForest(formula=RF.formula, ntree=RF.ntree, mtry=RF.mtry, data=TrainData, na.action=na.omit),
                error= function(err) {print(paste("random forest calibration failed"))},
                silent=T)
        }else{ 
            results <- RF.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "RF calibration","\n\n",sep = ""))
            TrainData[,"RF"] <- predict(object=results, newdata=TrainData.vars, type="response")
            TrainPres <- TrainData[TrainData[,"pb"]==1,"RF"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"RF"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            weights["RF"] <- max(c(eval1@auc, 0), na.rm=T)
#            TrainData[,"RF"] <- trunc(1000*TrainData[,"RF"])
            if (no.tests == F) {
                cat(paste("\n", "RF evaluation","\n\n",sep = ""))
                TestData[,"RF"] <- predict(object=results, newdata=TestData.vars, type="response")
                TestPres <- TestData[TestData[,"pb"]==1,"RF"]
                TestAbs <- TestData[TestData[,"pb"]==0,"RF"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("RF evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["RF"] <- max(c(eval2@auc, 0), na.rm=T)
#                    TestData[,"RF"] <- trunc(1000*TestData[,"RF"])
                    if(PLOTS==T) {
                        plot(eval2, "ROC") 
                        title(main="RF", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: RF evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$RF.C <- eval1
                evaluations$RF.T <- eval2
            }
            if(models.keep==T) {
                models$RF <- results
                models$formulae$RF.formula <- RF.formula
            }
        }else{ cat(paste("\n", "WARNING: random forest calibration failed", "\n", "\n"))}
#       detach(package:randomForest)
    }
    if (GLM > 0) {
        eval1 <- eval2 <- results <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(GLM.OLD) == T) { 
            tryCatch(results <- glm(formula=GLM.formula, family=GLM.family, data=TrainData, weights=Yweights, control=glm.control(maxit=maxit)),
                error= function(err) {print(paste("GLM calibration failed"))},
                silent=T)
        }else{ 
            results <- GLM.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "GLM calibration","\n\n",sep = ""))
            TrainData[,"GLM"] <- predict(object=results, newdata=TrainData.vars, type="response")
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GLM"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GLM"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            weights["GLM"] <- max(c(eval1@auc, 0), na.rm=T)
#            TrainData[,"GLM"] <- trunc(1000*TrainData[,"GLM"])
            if (no.tests == F) {
                cat(paste("\n", "GLM evaluation","\n\n",sep = ""))
                TestData[,"GLM"] <- predict(object=results, newdata=TestData.vars, type="response")
                TestPres <- TestData[TestData[,"pb"]==1,"GLM"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GLM"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("GLM evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["GLM"] <- max(c(eval2@auc, 0), na.rm=T)
#                    TestData[,"GLM"] <- trunc(1000*TestData[,"GLM"])
                    if(PLOTS==T) {
                        plot(eval2, "ROC") 
                        title(main="GLM", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: GLM evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$GLM.C <- eval1
                evaluations$GLM.T <- eval2
            }
            if(models.keep==T) {
                models$GLM <- results
                models$formulae$GLM.formula <- GLM.formula
            }
        }else{ cat(paste("\n", "WARNING: GLM calibration failed", "\n", "\n"))}
    }
    if (GLMSTEP > 0) {
#        require(MASS, quietly=T)
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(GLMSTEP.OLD) == T) {
            tryCatch(results <- glm(formula=STEP.formula, family=GLM.family, data=TrainData, weights=Yweights, control=glm.control(maxit=maxit)),
                error= function(err) {print(paste("first step of stepwise GLM calibration failed"))},
                silent=T)
            tryCatch(results2 <- stepAIC(results, scope=GLMSTEP.scope, direction="both", trace=F, steps=GLMSTEP.steps, k=GLMSTEP.k),
                error= function(err) {print(paste("stepwise GLM calibration failed"))},
                silent=T)
        }else{ 
            results2 <- GLMSTEP.OLD
        }
        if (is.null(results2) == F) {
            results <- results2
            cat(paste("\n", "stepwise GLM formula","\n\n",sep = ""))
            print(formula(results))
            cat(paste("\n", "stepwise GLM calibration","\n\n",sep = ""))
            TrainData[,"GLMSTEP"] <- predict(object=results, newdata=TrainData.vars, type="response")
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GLMSTEP"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GLMSTEP"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            weights["GLMSTEP"] <- max(c(eval1@auc, 0), na.rm=T)
#            TrainData[,"GLMSTEP"] <- trunc(1000*TrainData[,"GLMSTEP"])
            if (no.tests == F) {
                cat(paste("\n", "stepwise GLM evaluation","\n\n",sep = ""))
                TestData[,"GLMSTEP"] <- predict(object=results, newdata=TestData.vars, type="response")
                TestPres <- TestData[TestData[,"pb"]==1,"GLMSTEP"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GLMSTEP"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("stepwise GLM evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["GLMSTEP"] <- max(c(eval2@auc, 0), na.rm=T)
#                    TestData[,"GLMSTEP"] <- trunc(1000*TestData[,"GLMSTEP"])
                    if(PLOTS==T) {
                        plot(eval2, "ROC") 
                        title(main="STEPWISE GLM", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: stepwise GLM evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$GLMS.C <- eval1
                evaluations$GLMS.T <- eval2
            }
            if(models.keep==T) {
                models$GLMSTEP <- results
                models$formulae$STEP.formula <- STEP.formula
                models$formulae$GLMSTEP.scope <- GLMSTEP.scope
            }
        }else{ cat(paste("\n", "WARNING: stepwise GBM calibration failed", "\n", "\n"))}
#       detach(package:MASS)
    }
    if (GAM > 0) {
        cat(paste("\n\n"))
        require(gam, quietly=T)
        eval1 <- eval2 <- results <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(GAM.OLD) == T) {
            tryCatch(results <- gam(formula=GAM.formula, family=GAM.family, data=TrainData, weights=Yweights, control=gam.control(maxit=maxit, bf.maxit=50)),
                error= function(err) {print(paste("GAM calibration (gam package) failed"))},
                silent=T)
        }else{ 
            results <- GAM.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "GAM calibration (gam package)","\n\n",sep = ""))
            TrainData[,"GAM"] <- predict(object=results, newdata=TrainData.vars, type="response")
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GAM"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GAM"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            weights["GAM"] <- max(c(eval1@auc, 0), na.rm=T)
#            TrainData[,"GAM"] <- trunc(1000*TrainData[,"GAM"])
            if (no.tests == F) {
                cat(paste("\n", "GAM evaluation (gam package)","\n\n",sep = ""))
                TestData[,"GAM"] <- predict(object=results, newdata=TestData.vars, type="response")
                TestPres <- TestData[TestData[,"pb"]==1,"GAM"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GAM"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("GAM evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["GAM"] <- max(c(eval2@auc, 0), na.rm=T)
#                    TestData[,"GAM"] <- trunc(1000*TestData[,"GAM"])
                    if(PLOTS==T) {
                        plot(eval2, "ROC") 
                        title(main="GAM (gam)", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: GAM evaluation (gam package) failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$GAM.C <- eval1
                evaluations$GAM.T <- eval2
            }
            if(models.keep==T) {
                models$GAM <- results
                models$formulae$GAM.formula <- GAM.formula
            }
        }else{ cat(paste("\n", "WARNING: GAM calibration (gam package) failed", "\n", "\n"))}
# conflict with mgcv gam function
        detach(package:gam)
    }
    if (GAMSTEP > 0) {
        cat(paste("\n\n"))
        require(gam, quietly=T)
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
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
        }else{ 
            results2 <- GAMSTEP.OLD
        }
        if (is.null(results2) == F) {
            results <- results2
            cat(paste("\n", "stepwise GAM formula (gam package)","\n\n",sep = ""))
            print(formula(results))
            cat(paste("\n", "stepwise GAM calibration (gam package)","\n\n",sep = ""))
            TrainData[,"GAMSTEP"] <- predict(object=results, newdata=TrainData.vars, type="response")
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GAMSTEP"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GAMSTEP"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            weights["GAMSTEP"] <- max(c(eval1@auc, 0), na.rm=T)
#            TrainData[,"GAMSTEP"] <- trunc(1000*TrainData[,"GAMSTEP"])
            if (no.tests == F) {
                cat(paste("\n", "stepwise GAM evaluation (gam package)","\n\n","\n", sep = ""))
                TestData[,"GAMSTEP"] <- predict(object=results, newdata=TestData.vars, type="response")
                TestPres <- TestData[TestData[,"pb"]==1,"GAMSTEP"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GAMSTEP"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("stepwise GAM evaluation (gam package) failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["GAMSTEP"] <- max(c(eval2@auc, 0), na.rm=T)
#                    TestData[,"GAMSTEP"] <- trunc(1000*TestData[,"GAMSTEP"])
                    if(PLOTS==T) {
                        plot(eval2, "ROC") 
                        title(main="STEPWISE GAM (gam)", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: stepwise GAM evaluation (gam package) failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$GAMS.C <- eval1
                evaluations$GAMS.T <- eval2
            }
            if(models.keep==T) {
                models$GAMSTEP <- results
                models$formulae$STEP.formula <- STEP.formula
                models$formulae$GAMSTEP.scope <- GAMSTEP.scope
            }
        }else{ cat(paste("\n", "WARNING: stepwise GAM calibration (gam package) failed", "\n", "\n"))}
# conflict with mgcv gam function
        detach(package:gam)
    }
    if (MGCV > 0) {
        cat(paste("\n\n"))
        require(mgcv, quietly=T)
        eval1 <- eval2 <- results <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(MGCV.OLD) == T) {
            tryCatch(results <- gam(formula=MGCV.formula, family=GAM.family, data=TrainData, weights=Yweights, 
                        select=MGCV.select, control=gam.control(maxit=maxit)),
                error= function(err) {print(paste("GAM calibration (mgcv package) failed"))},
                silent=T)
        }else{ 
            results <- MGCV.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "GAM calibration (mgcv package)","\n\n",sep = ""))
            TrainData[,"MGCV"] <- predict.mgcv(object=results, newdata=TrainData.vars)
            TrainPres <- TrainData[TrainData[,"pb"]==1,"MGCV"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"MGCV"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            weights["MGCV"] <- max(c(eval1@auc, 0), na.rm=T)
#            TrainData[,"MGCV"] <- trunc(1000*TrainData[,"MGCV"])
            if (no.tests == F) {
                cat(paste("\n", "GAM evaluation (mgcv package)","\n\n",sep = ""))
                TestData[,"MGCV"] <- predict.mgcv(object=results, newdata=TestData.vars)
                TestPres <- TestData[TestData[,"pb"]==1,"MGCV"]
                TestAbs <- TestData[TestData[,"pb"]==0,"MGCV"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("GAM evaluation (mgcv package) failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["MGCV"] <- max(c(eval2@auc, 0), na.rm=T)
#                    TestData[,"MGCV"] <- trunc(1000*TestData[,"MGCV"])
                    if(PLOTS==T) {
                        plot(eval2, "ROC") 
                        title(main="GAM (mgcv)", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: GAM evaluation (mgcv package) failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$MGCV.C <- eval1
                evaluations$MGCV.T <- eval2
            }
            if(models.keep==T) {
                models$MGCV <- results
                models$formulae$MGCV.formula <- MGCV.formula
            }
        }else{ cat(paste("\n", "WARNING: GAM calibration (mgcv package) failed", "\n", "\n"))}
# conflict with gam package
        detach(package:mgcv)
    }
    if (MGCVFIX > 0) {
        cat(paste("\n\n"))
        require(mgcv, quietly=T)
        eval1 <- eval2 <- results <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(MGCVFIX.OLD) == T) {
            tryCatch(results <- gam(formula=MGCVFIX.formula, family=GAM.family, data=TrainData, weights=Yweights, select=FALSE, control=gam.control(maxit=maxit)),
                error= function(err) {print(paste("GAM calibration with fixed d.f. regression splines (mgcv package) failed"))},
                silent=T)
        }else{ 
            results <- MGCVFIX.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "GAM with fixed d.f. regression splines calibration (mgcv package)","\n\n",sep = ""))
            TrainData[,"MGCVFIX"] <- predict.mgcv(object=results, newdata=TrainData.vars)
            TrainPres <- TrainData[TrainData[,"pb"]==1,"MGCVFIX"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"MGCVFIX"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            weights["MGCVFIX"] <- max(c(eval1@auc, 0), na.rm=T)
#            TrainData[,"MGCVFIX"] <- trunc(1000*TrainData[,"MGCVFIX"])
            if (no.tests == F) {
                cat(paste("\n", "GAM with fixed d.f. regression splines evaluation (mgcv package)","\n\n",sep = ""))
                TestData[,"MGCVFIX"] <- predict.mgcv(object=results, newdata=TestData)
                TestPres <- TestData[TestData[,"pb"]==1,"MGCVFIX"]
                TestAbs <- TestData[TestData[,"pb"]==0,"MGCVFIX"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("GAMFIX evaluation (mgcv package) failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["MGCVFIX"] <- max(c(eval2@auc, 0), na.rm=T)
#                    TestData[,"MGCVFIX"] <- trunc(1000*TestData[,"MGCVFIX"])
                    if(PLOTS==T) {
                        plot(eval2, "ROC") 
                        title(main="GAM (mgcv)", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: GAM with fixed d.f. regression splines evaluation (mgcv package) failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$MGCVF.C <- eval1
                evaluations$MGCVF.T <- eval2
            }
            if(models.keep==T) {
                models$MGCVFIX <- results
                models$formulae$MGCVFIX.formula <- MGCVFIX.formula
            }
        }else{ cat(paste("\n", "WARNING: MGCVFIX calibration (mgcv package) failed", "\n", "\n"))}
# conflict with gam package
        detach(package:mgcv)
    }
    if (!is.null(factors) && EARTH > 0) {
        cat(paste("\n", "NOTE: MARS evaluation (earth package) with factors probably requires dummy variables", sep=""))
    } 
    if (EARTH > 0) {
#        require(earth, quietly=T)
        eval1 <- eval2 <- results <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(EARTH.OLD) == T) {
            tryCatch(results <- earth(formula=EARTH.formula, glm=EARTH.glm, data=TrainData, degree=2),
                error= function(err) {print(paste("MARS calibration (earth package) failed"))},
                silent=T)
        }else{ 
            results <- EARTH.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "MARS calibration (earth package)","\n\n",sep = ""))
            TrainData[,"EARTH"] <- predict.earth2(object=results, newdata=TrainData.vars)
            TrainPres <- TrainData[TrainData[,"pb"]==1,"EARTH"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"EARTH"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            weights["EARTH"] <- max(c(eval1@auc, 0), na.rm=T)
#            TrainData[,"EARTH"] <- trunc(1000*TrainData[,"EARTH"])
            if (no.tests == F) {
                cat(paste("\n", "MARS evaluation (earth package)","\n\n",sep = ""))
                TestData[,"EARTH"] <- predict.earth2(object=results, newdata=TestData.vars)
                TestPres <- TestData[TestData[,"pb"]==1,"EARTH"]
                TestAbs <- TestData[TestData[,"pb"]==0,"EARTH"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("MARS evaluation (earth package) failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["EARTH"] <- max(c(eval2@auc, 0), na.rm=T)
#                    TestData[,"EARTH"] <- trunc(1000*TestData[,"EARTH"])
                    if(PLOTS==T) {
                        plot(eval2, "ROC") 
                        title(main="MARS (earth)", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: MARS evaluation (earth package) failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$EARTH.C <- eval1
                evaluations$EARTH.T <- eval2
            }
            if(models.keep==T) {
                models$EARTH <- results
                models$formulae$EARTH.formula <- EARTH.formula
            }
        }else{ cat(paste("\n", "WARNING: MARS calibration (earth package) failed", "\n", "\n"))}
#       detach(package:earth)
    }
    if (RPART > 0) {
#        require(rpart, quietly=T)
        eval1 <- eval2 <- results <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(RPART.OLD) == T) {
            tryCatch(results <- rpart(formula=RPART.formula, data=TrainData, weights=Yweights,
                    control=rpart.control(xval=RPART.xval, minbucket=5, minsplit=5, cp=0.001, maxdepth=25)),
                error= function(err) {print(paste("RPART calibration failed"))},
                silent=T)
        }else{ 
            results <- RPART.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "RPART calibration","\n\n",sep = ""))
            TrainData[,"RPART"] <- predict(object=results, newdata=TrainData.vars, type="prob")[,2]
            TrainPres <- TrainData[TrainData[,"pb"]==1,"RPART"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"RPART"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            weights["RPART"] <- max(c(eval1@auc, 0), na.rm=T)
#            TrainData[,"RPART"] <- trunc(1000*TrainData[,"RPART"])
            if (no.tests == F) {
                cat(paste("\n", "RPART evaluation","\n\n",sep = ""))
                TestData[,"RPART"] <- predict(object=results, newdata=TestData.vars, type="prob")[,2]
                TestPres <- TestData[TestData[,"pb"]==1,"RPART"]
                TestAbs <- TestData[TestData[,"pb"]==0,"RPART"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("RPART evaluation failed"))},
                    silent=F)
                if (is.null(TestPres) == F && is.null(TestAbs) == F) {
                    eval2 <-  evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
                    weights["RPART"] <- max(c(eval2@auc, 0), na.rm=T)
#                    TestData[,"RPART"] <- trunc(1000*TestData[,"RPART"])
                    if(PLOTS==T) {
                        plot(eval2, "ROC") 
                        title(main="RPART", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: RPART evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$RPART.C <- eval1
                evaluations$RPART.T <- eval2
            }
            if(models.keep==T) {
                models$RPART <- results
                models$formulae$RPART.formula <- RPART.formula
            }
        }else{ cat(paste("\n", "WARNING: RPART calibration failed", "\n", "\n"))}
#       detach(package:rpart)
    }
    if (NNET > 0) {
#        require(nnet, quietly=T)
        eval1 <- eval2 <- results <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(NNET.OLD) == T) {
            tryCatch(results <- nnet(formula=NNET.formula, size=NNET.size, decay=NNET.decay, data=TrainData, weights=Yweights, 
                    rang=0.1, maxit=maxit, trace=F),
                error= function(err) {print(paste("ANN calibration (nnet package) failed"))},
                silent=T)
        }else{ 
            results <- NNET.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "ANN calibration (nnet package)","\n\n",sep = ""))
            TrainData[,"NNET"] <- predict.nnet2(object=results, newdata=TrainData.vars)
            TrainPres <- TrainData[TrainData[,"pb"]==1,"NNET"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"NNET"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            weights["NNET"] <- max(c(eval1@auc, 0), na.rm=T)
#            TrainData[,"NNET"] <- trunc(1000*TrainData[,"NNET"])
            if (no.tests == F) {
                cat(paste("\n", "ANN evaluation (nnet package)","\n\n",sep = ""))
                TestData[,"NNET"] <- predict.nnet2(object=results, newdata=TestData.vars)
                TestPres <- TestData[TestData[,"pb"]==1,"NNET"]
                TestAbs <- TestData[TestData[,"pb"]==0,"NNET"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("ANN evaluation (nnet package) failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["NNET"] <- max(c(eval2@auc, 0), na.rm=T)
#                    TestData[,"NNET"] <- trunc(1000*TestData[,"NNET"])
                    if(PLOTS==T) {
                        plot(eval2, "ROC") 
                        title(main="NNET", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: ANN evaluation (nnet package) failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$NNET.C <- eval1
                evaluations$NNET.T <- eval2
            }
            if(models.keep==T) {
                models$NNET <- results
                models$formulae$NNET.formula <- NNET.formula
            }
        }else{ cat(paste("\n", "WARNING: ANN calibration (nnet package) failed", "\n", "\n"))}
#       detach(package:nnet)
    }
    if (FDA > 0) {
#        require(mda, quietly=T)
        eval1 <- eval2 <- results <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(FDA.OLD) == T) {
            tryCatch(results <- fda(formula=FDA.formula, method=mars, data=TrainData, weights=Yweights),
                error= function(err) {print(paste("FDA calibration failed"))},
                silent=T)
        }else{ 
            results <- FDA.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "FDA calibration","\n\n",sep = ""))
            TrainData[,"FDA"] <- predict(object=results, newdata=TrainData.vars, type="posterior")[,2]
            TrainPres <- TrainData[TrainData[,"pb"]==1,"FDA"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"FDA"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            weights["FDA"] <- max(c(eval1@auc, 0), na.rm=T)
#            TrainData[,"FDA"] <- trunc(1000*TrainData[,"FDA"])
            if (no.tests == F) {
                cat(paste("\n", "FDA evaluation","\n\n",sep = ""))
                TestData[,"FDA"] <- predict(object=results, newdata=TestData.vars, type="posterior")[,2]
                TestPres <- TestData[TestData[,"pb"]==1,"FDA"]
                TestAbs <- TestData[TestData[,"pb"]==0,"FDA"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("FDA evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["FDA"] <- max(c(eval2@auc, 0), na.rm=T)
#                    TestData[,"FDA"] <- trunc(1000*TestData[,"FDA"])
                    if(PLOTS==T) {
                        plot(eval2, "ROC") 
                        title(main="FDA", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: FDA evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$FDA.C <- eval1
                evaluations$FDA.T <- eval2
            }
            if(models.keep==T) {
                models$FDA <- results
                models$formulae$FDA.formula <- FDA.formula
            }
        }else{ cat(paste("\n", "WARNING: FDA calibration failed", "\n", "\n"))}
#       detach(package:mda)
    }
    if (SVM > 0) {
        cat(paste("\n\n"))
#        require(kernlab, quietly=T)
        eval1 <- eval2 <- results <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(SVM.OLD) == T) {
            tryCatch(results <- ksvm(SVM.formula, data=TrainData, type="C-svc", prob.model=T),
                error= function(err) {print(paste("SVM calibration failed"))},
                silent=T)
        }else{ 
            results <- SVM.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "SVM calibration (kernlab package)","\n\n",sep = ""))
            TrainData[,"SVM"] <- kernlab::predict(object=results, newdata=TrainData.vars, type="probabilities")[,2]
            TrainPres <- TrainData[TrainData[,"pb"]==1,"SVM"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"SVM"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            weights["SVM"] <- max(c(eval1@auc, 0), na.rm=T)
#            TrainData[,"SVM"] <- trunc(1000*TrainData[,"SVM"])
            if (no.tests == F) {
                cat(paste("\n", "SVM evaluation (kernlab package)","\n\n",sep = ""))
                TestData[,"SVM"] <- kernlab::predict(object=results, newdata=TestData.vars, type="probabilities")[,2]
                TestPres <- TestData[TestData[,"pb"]==1,"SVM"]
                TestAbs <- TestData[TestData[,"pb"]==0,"SVM"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("SVM evaluation (kernlab package) failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["SVM"] <- max(c(eval2@auc, 0), na.rm=T)
#                    TestData[,"SVM"] <- trunc(1000*TestData[,"SVM"])
                    if(PLOTS==T) {
                        plot(eval2, "ROC") 
                        title(main="SVM (kernlab package)", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: SVM evaluation (kernlab package) failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$SVM.C <- eval1
                evaluations$SVM.T <- eval2
            }
            if(models.keep==T) {
                models$SVM <- results
                models$formulae$SVM.formula <- SVM.formula
            }
        }else{ cat(paste("\n", "WARNING: SVM calibration failed", "\n", "\n"))}
#       detach(package:kernlab)
    }
    if (SVME > 0) {
#        require(e1071, quietly=T)
        eval1 <- eval2 <- results <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(SVME.OLD) == T) {
            tryCatch(results <- svm(SVME.formula, data=TrainData, type="C-classification", kernel="polynomial", degree=3, probability=TRUE),
                error= function(err) {print(paste("SVM calibration (e1071 package) failed"))},
                silent=T)
        }else{ 
            results <- SVME.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "SVM calibration (e1071 package)","\n\n",sep = ""))
            TrainData[,"SVME"] <- predict.svme(model=results, newdata=TrainData.vars)
            TrainPres <- TrainData[TrainData[,"pb"]==1,"SVME"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"SVME"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            weights["SVME"] <- max(c(eval1@auc, 0), na.rm=T)
#            TrainData[,"SVME"] <- trunc(1000*TrainData[,"SVME"])
            if (no.tests == F) {
                cat(paste("\n", "SVM evaluation (e1071 package)","\n\n",sep = ""))
                TestData[,"SVME"] <- predict.svme(model=results, newdata=TestData.vars)
                TestPres <- TestData[TestData[,"pb"]==1,"SVME"]
                TestAbs <- TestData[TestData[,"pb"]==0,"SVME"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("SVM evaluation (e1071 package) failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["SVME"] <- max(c(eval2@auc, 0), na.rm=T)
#                    TestData[,"SVME"] <- trunc(1000*TestData[,"SVME"])
                    if(PLOTS==T) {
                        plot(eval2, "ROC") 
                        title(main="SVM (e1017 package)", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: SVM evaluation (e1071 package) failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$SVME.C <- eval1
                evaluations$SVME.T <- eval2
            }
            if(models.keep==T) {
                models$SVME <- results
                models$formulae$SVME.formula <- SVME.formula
            }
        }else{ cat(paste("\n", "WARNING: SVM calibration (e1071 package) failed", "\n", "\n"))}
#       detach(package:e1071)
    }
    if (BIOCLIM > 0 || DOMAIN > 0 || MAHAL > 0) {
        if(is.null(factors)==F) {
            for (i in 1:length(factors)) {
                TrainData.vars <- TrainData.vars[, which(colnames(TrainData.vars) != factors[i])]
                TrainData.pres <- TrainData.pres[, which(colnames(TrainData.pres) != factors[i])]
                TestData.vars <- TestData.vars[, which(colnames(TestData.vars) != factors[i])]
            }
        }
    }
    if (BIOCLIM > 0) {
        eval1 <- eval2 <- results <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(BIOCLIM.OLD) == T) {
            tryCatch(results <- bioclim(x=TrainData.pres),
                error= function(err) {print(paste("BIOCLIM calibration failed"))},
                silent=T)
        }else{ 
            results <- BIOCLIM.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "BIOCLIM calibration","\n\n",sep = ""))
            TrainData[,"BIOCLIM"] <- dismo::predict(object=results, x=TrainData.vars)
            TrainPres <- TrainData[TrainData[,"pb"]==1,"BIOCLIM"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"BIOCLIM"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            weights["BIOCLIM"] <- max(c(eval1@auc, 0), na.rm=T)
#            TrainData[,"BIOCLIM"] <- trunc(1000*TrainData[,"BIOCLIM"])
            if (no.tests == F) {
                cat(paste("\n", "BIOCLIM evaluation","\n\n",sep = ""))
                TestData[,"BIOCLIM"] <- dismo::predict(object=results, x=TestData.vars)
                TestPres <- TestData[TestData[,"pb"]==1,"BIOCLIM"]
                TestAbs <- TestData[TestData[,"pb"]==0,"BIOCLIM"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("BIOCLIM evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["BIOCLIM"] <- max(c(eval2@auc, 0), na.rm=T)
#                    TestData[,"BIOCLIM"] <- trunc(1000*TestData[,"BIOCLIM"])
                    if(PLOTS==T) {
                        plot(eval2, "ROC") 
                        title(main="BIOCLIM", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: BIOCLIM evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$BIOCLIM.C <- eval1
                evaluations$BIOCLIM.T <- eval2
            }
            if(models.keep==T) {
                models$BIOCLIM <- results
            }
        }else{ cat(paste("\n", "WARNING: BIOCLIM calibration failed", "\n", "\n"))}
    }
    if (DOMAIN > 0) {
        eval1 <- eval2 <- results <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(DOMAIN.OLD) == T) {
            tryCatch(results <- domain(x=TrainData.pres),
                error= function(err) {print(paste("DOMAIN calibration failed"))},
                silent=T)
        }else{ 
            results <- DOMAIN.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "DOMAIN calibration","\n\n",sep = ""))
            TrainData[,"DOMAIN"] <- dismo::predict(object=results, x=TrainData.vars)
            TrainPres <- TrainData[TrainData[,"pb"]==1,"DOMAIN"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"DOMAIN"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            weights["DOMAIN"] <- max(c(eval1@auc, 0), na.rm=T)
#            TrainData[,"DOMAIN"] <- trunc(1000*TrainData[,"DOMAIN"])
            if (no.tests == F) {
                cat(paste("\n", "DOMAIN evaluation","\n\n",sep = ""))
                TestData[,"DOMAIN"] <- dismo::predict(object=results, x=TestData.vars)
                TestPres <- TestData[TestData[,"pb"]==1,"DOMAIN"]
                TestAbs <- TestData[TestData[,"pb"]==0,"DOMAIN"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("DOMAIN evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["DOMAIN"] <- max(c(eval2@auc, 0), na.rm=T)
#                    TestData[,"DOMAIN"] <- trunc(1000*TestData[,"DOMAIN"])
                    if(PLOTS==T) {
                        plot(eval2, "ROC") 
                        title(main="DOMAIN", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: DOMAIN evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$DOMAIN.C <- eval1
                evaluations$DOMAIN.T <- eval2
            }
            if(models.keep==T) {
                models$DOMAIN <- results
            }
        }else{ cat(paste("\n", "WARNING: DOMAIN calibration failed", "\n", "\n"))}
    }
    if (MAHAL > 0) {
        eval1 <- eval2 <- results <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(MAHAL.OLD) == T) {
            tryCatch(results <- mahal(x=TrainData.pres),
                error= function(err) {print(paste("Mahalanobis calibration failed"))},
                silent=T)
        }else{ 
            results <- MAHAL.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Mahalanobis calibration ","\n\n",sep = ""))
            TrainData[,"MAHAL"] <- predict.mahal(model=results, newdata=TrainData.vars, MAHAL.shape=MAHAL.shape)
            TrainPres <- TrainData[TrainData[,"pb"]==1,"MAHAL"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"MAHAL"]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            weights["MAHAL"] <- max(c(eval1@auc, 0), na.rm=T)
#            TrainData[,"MAHAL"] <- trunc(1000*TrainData[,"MAHAL"])
            if (no.tests == F) {
                cat(paste("\n", "Mahalanobis evaluation","\n\n",sep = ""))
                TestData[,"MAHAL"] <- predict.mahal(model=results, newdata=TestData.vars, MAHAL.shape=MAHAL.shape)
                TestPres <- TestData[TestData[,"pb"]==1,"MAHAL"]
                TestAbs <- TestData[TestData[,"pb"]==0,"MAHAL"]
                tryCatch(eval2 <- evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("MAHAL evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["MAHAL"] <- max(c(eval2@auc, 0), na.rm=T)
#                    TestData[,"MAHAL"] <- trunc(1000*TestData[,"MAHAL"])
                    if(PLOTS==T) {
                        plot(eval2, "ROC") 
                        title(main="MAHAL", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: Mahalanobis evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$MAHAL.C <- eval1
                evaluations$MAHAL.T <- eval2
            }
            if(models.keep==T) {
                models$MAHAL <- results
            }
        }else{ cat(paste("\n", "WARNING: Mahalanobis calibration failed", "\n", "\n"))}
    }
    if (GEODIST > 0) {
        eval1 <- eval2 <- results <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(GEODIST.OLD) == T) {
            tryCatch(results <- geoDist(p=p, lonlat=TRUE),
                error= function(err) {print(paste("GEODIST calibration failed"))},
                silent=F)
        }else{ 
            results <- GEODIST.OLD
        }
        if (is.null(results) == F) {
            NAmask <- crop(x[[1]], x[[1]])
            fullname <- paste("models/", GEODIST.file.name, "_GEO", sep="")
            pgeo <- dismo::predict(object=results, x=NAmask, mask=TRUE, ext=ext, 
                 filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format)
            cat(paste("\n", "GEODIST calibration","\n\n",sep = ""))
            pres_geo <- extract(pgeo, p)
            abs_geo <- extract(pgeo, a)
            eval1 <- evaluate(p=pres_geo, a=abs_geo, tr=quantile(pgeo, na.rm=T, probs=c(1:50/50), names=F))
            print(eval1)
            if (no.tests == F) {
                cat(paste("\n", "GEODIST evaluation","\n\n","\n", sep = ""))
                prest_geo <- extract(pgeo, pt)
                abst_geo <- extract(pgeo, at)
                tryCatch(eval2 <- evaluate(p=prest_geo, a=abst_geo, tr=quantile(pgeo, na.rm=T, probs=c(1:50/50), names=F)),
                    error= function(err) {print(paste("GEODIST evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    if(PLOTS==T) {
                        plot(eval2, "ROC") 
                        title(main="GEODIST", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: GEODIST evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$GEODIST.C <- eval1
                evaluations$GEODIST.T <- eval2
            }
            if(models.keep==T) {
                models$GEODIST <- results
            }
        }else{ cat(paste("\n", "WARNING: GEODIST calibration failed", "\n", "\n"))}
    }
    ENSEMBLE.tune <- FALSE
    if((length(ENSEMBLE.decay) == 1) && (length(ENSEMBLE.best) == 1) && (length(ENSEMBLE.min) == 1)) {
        ws <- ensemble.weights(weights, decay=ENSEMBLE.decay, 
            best=ENSEMBLE.best, min.weight=ENSEMBLE.min, scale=TRUE, multiply=TRUE)
#        names(ws) <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX", 
#            "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL")
        cat(paste("\n", "Weights for ensemble forecasting", sep = ""))
        cat(paste("\n", "(minimum AUC was ", ENSEMBLE.min, ")", "\n\n", sep = ""))
        print(ws)
        if(evaluations.keep == T) {evaluations$ensemble.weights <- ws}
    }else{
# use different strategies for calculating the ensemble model
# similar to using test data for calculating input AUC, use test data for calculating best ensemble
# recalculating AUC does not require much computing time - initial calculations kept to spot problems for specific algorithms
        strategy.results <- ensemble.strategy(TrainData=TrainData, TestData=TestData, 
            ENSEMBLE.decay=ENSEMBLE.decay, ENSEMBLE.best=ENSEMBLE.best, ENSEMBLE.min=ENSEMBLE.min)
        ws <- strategy.results$weights
        ENSEMBLE.tune <- TRUE
#        names(ws) <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX", 
#            "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL")
        cat(paste("\n", "Weights for ensemble forecasting", "\n\n", sep = ""))
        print(ws)
        if(evaluations.keep == T) {evaluations$STRATEGY.weights <- ws}
    }
# do not return ensemble tests for only one model  
    if((sum(ws > 0, na.rm=T) > 1) || (ENSEMBLE.tune == T)) {
        TrainData[,"ENSEMBLE"] <- ws["MAXENT"]*TrainData[,"MAXENT"] + ws["GBM"]*TrainData[,"GBM"] +
            ws["GBMSTEP"]*TrainData[,"GBMSTEP"] + ws["RF"]*TrainData[,"RF"] + ws["GLM"]*TrainData[,"GLM"] +
            ws["GLMSTEP"]*TrainData[,"GLMSTEP"] + ws["GAM"]*TrainData[,"GAM"] + ws["GAMSTEP"]*TrainData[,"GAMSTEP"] +
            ws["MGCV"]*TrainData[,"MGCV"] + ws["MGCVFIX"]*TrainData[,"MGCVFIX"] + ws["EARTH"]*TrainData[,"EARTH"] +
            ws["RPART"]*TrainData[,"RPART"] + ws["NNET"]*TrainData[,"NNET"] + ws["FDA"]*TrainData[,"FDA"] +
            ws["SVM"]*TrainData[,"SVM"] + ws["SVME"]*TrainData[,"SVME"] + ws["BIOCLIM"]*TrainData[,"BIOCLIM"] +
            ws["DOMAIN"]*TrainData[,"DOMAIN"] + ws["MAHAL"]*TrainData[,"MAHAL"]
#        TrainData[,"ENSEMBLE"] <- trunc(TrainData[,"ENSEMBLE"])
        TestData[,"ENSEMBLE"] <- ws["MAXENT"]*TestData[,"MAXENT"] + ws["GBM"]*TestData[,"GBM"] +
            ws["GBMSTEP"]*TestData[,"GBMSTEP"] + ws["RF"]*TestData[,"RF"] + ws["GLM"]*TestData[,"GLM"] +
            ws["GLMSTEP"]*TestData[,"GLMSTEP"] + ws["GAM"]*TestData[,"GAM"] + ws["GAMSTEP"]*TestData[,"GAMSTEP"] +
            ws["MGCV"]*TestData[,"MGCV"] + ws["MGCVFIX"]*TestData[,"MGCVFIX"] + ws["EARTH"]*TestData[,"EARTH"] +
            ws["RPART"]*TestData[,"RPART"] + ws["NNET"]*TestData[,"NNET"] + ws["FDA"]*TestData[,"FDA"] +
            ws["SVM"]*TestData[,"SVM"] + ws["SVME"]*TestData[,"SVME"] + ws["BIOCLIM"]*TestData[,"BIOCLIM"] +
            ws["DOMAIN"]*TestData[,"DOMAIN"] + ws["MAHAL"]*TestData[,"MAHAL"]
#        TestData[,"ENSEMBLE"] <- trunc(TestData[,"ENSEMBLE"])
        eval1 <- eval2 <- NULL
        cat(paste("\n", "Ensemble evaluation with calibration data", "\n\n", sep = ""))
        TrainPres <- as.numeric(TrainData[TrainData[,"pb"]==1,"ENSEMBLE"])
        TrainAbs <- as.numeric(TrainData[TrainData[,"pb"]==0,"ENSEMBLE"])
        eval1 <- evaluate(p=TrainPres, a=TrainAbs)
        print(eval1)
        if (no.tests == F) {
            cat(paste("\n", "Ensemble evaluation with testing data", "\n\n", sep = ""))
            TestPres <- as.numeric(TestData[TestData[,"pb"]==1,"ENSEMBLE"])
            TestAbs <- as.numeric(TestData[TestData[,"pb"]==0,"ENSEMBLE"])
            eval2 <- evaluate(p=TestPres, a=TestAbs)
            print(eval2)
            if(PLOTS==T) {
                plot(eval2, "ROC") 
                title(main="ENSEMBLE", cex=2, adj=0, col.main="blue")
            }
        }
        if(evaluations.keep==T) {
            evaluations$ENSEMBLE.C <- eval1
            evaluations$ENSEMBLE.T <- eval2
            evaluations$TrainData <- TrainData
            evaluations$TestData <- TestData
        }
    }
    cat(paste("\n", "End of evaluations", "\n\n",sep = ""))
    remove(Yweights, envir=.BiodiversityR)
    remove(TrainData, envir=.BiodiversityR)
    remove(TestData, envir=.BiodiversityR)
    if (evaluations.keep==T || models.keep==T) {return(list(evaluations=evaluations, models=models, call=match.call() ))}
}

