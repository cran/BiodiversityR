`ensemble.test` <- function(
    x, p, a=NULL, an=1000, ext=NULL, k=0, pt=NULL, at=NULL,
    layer.drops=NULL, VIF=FALSE, COR=FALSE,
    PLOTS=TRUE, evaluations.keep=FALSE, models.keep=FALSE, 
    input.weights=NULL,
    MAXENT=1, GBM=1, GBMSTEP=1, RF=1, GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=1, MGCV=1, MGCVFIX=0,
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, BIOCLIM=1, DOMAIN=1, MAHAL=1, 
    GEODIST=0, GEODIST.file.name="Species001", RASTER.format="raster", 
    Yweights="BIOMOD", factors=NULL, dummy.vars=NULL,
    formulae.defaults=TRUE, maxit=100,
    GBM.formula=NULL, GBM.n.trees=2000,
    GBMSTEP.gbm.x=2:(1+nlayers(x)), GBMSTEP.tree.complexity=5, GBMSTEP.learning.rate=0.005, 
    GBMSTEP.bag.fraction=0.5, GBMSTEP.step.size=100,
    RF.formula=NULL, RF.ntree=750, RF.mtry=log(nlayers(x)), 
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
    SVM.formula=NULL, SVME.formula=NULL
)
{
    .BiodiversityR <- new.env()
    if (! require(dismo)) {stop("Please install the dismo package")}
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
    if (is.null(input.weights)==F) {
# use the last column in case output from the ensemble.test.splits function is used
        if (length(dim(input.weights)) == 2) {
            input.weights <- input.weights[,"MEAN"]
            input.weights <- input.weights - 50
            input.weights[input.weights < 0] <- 0
        }
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
    if (MAXENT > 0) {
	jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
	if (!file.exists(jar)) {stop('maxent program is missing: ', jar, '\nPlease download it here: http://www.cs.princeton.edu/~schapire/maxent/')}
    }
    if (formulae.defaults == T) {
        formulae <- ensemble.formulae(x, factors=factors, dummy.vars=dummy.vars)
    }
    if (GBM > 0) {
        if (! require(gbm)) {stop("Please install the gbm package")}
        if (is.null(GBM.formula) == T && formulae.defaults == T) {GBM.formula <- formulae$GBM.formula}
        if (is.null(GBM.formula) == T) {stop("Please provide the GBM.formula (hint: use ensemble.formulae function)")}
        environment(GBM.formula) <- .BiodiversityR
    }
    if (GBMSTEP > 0) {
        if (! require(gbm)) {stop("Please install the gbm package")}
        if (is.null(GBM.formula) == T && formulae.defaults == T) {GBM.formula <- formulae$GBM.formula}
        if (is.null(GBM.formula) == T) {stop("Please provide the GBM.formula (hint: use ensemble.formulae function)")}
        environment(GBM.formula) <- .BiodiversityR
    }
    if (RF > 0) {
        if (! require(randomForest)) {stop("Please install the randomForest package")}
        if (is.null(RF.formula) == T && formulae.defaults == T) {RF.formula <- formulae$RF.formula}
        if (is.null(RF.formula) == T) {stop("Please provide the RF.formula (hint: use ensemble.formulae function)")}
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
    }
    if (MGCVFIX > 0) {
        if (! require(mgcv)) {stop("Please install the mgcv package")}
        if (is.null(MGCVFIX.formula) == T && formulae.defaults == T) {MGCVFIX.formula <- formulae$MGCVFIX.formula}
        if (is.null(MGCVFIX.formula) == T) {stop("Please provide the MGCVFIX.formula (hint: use ensemble.formulae function)")}
        environment(MGCVFIX.formula) <- .BiodiversityR
        assign("GAM.family", GAM.family, envir=.BiodiversityR)
        detach(package:mgcv)   
    }
    if (EARTH > 0) {
        if (! require(earth)) {stop("Please install the earth package")}
        if (is.null(EARTH.formula) == T && formulae.defaults == T) {EARTH.formula <- formulae$EARTH.formula}
        if (is.null(EARTH.formula) == T) {stop("Please provide the EARTH.formula (hint: use ensemble.formulae function)")}
        environment(EARTH.formula) <- .BiodiversityR
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
    if (is.null(a)==T) {
        a <- randomPoints(x, n=an, ext=ext)
    }
    if (is.null(pt)==T && k < 2) {
        pt <- p
    }
    if (is.null(pt)==T && k > 1) {
        groupp <- kfold(p, k=k)
        pc <- p[groupp != 1,]
        pt <- p[groupp == 1,]
        p <- pc
    }
    if (is.null(at)==T && k < 2) {
        at <- a
    }
    if (is.null(at)==T && k > 1) {
        groupa <- kfold(a, k=k)
        ac <- a[groupa != 1,]
        at <- a[groupa == 1,]
        a <- ac
    }
    if(is.null(factors)==F) {
        for (i in 1:length(factors)) {
            j <- which(names(x) == factors[i])
            x[[j]] <- raster::as.factor(x[[j]])
#            levels1 <- as.numeric(levels(as.factor(predictors[[j]]))[[1]]$ID)
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
    if(is.null(factors)==F) {
        for (i in 1:length(factors)) {
            if (identical(levels(TrainData[,factors[i]]),levels(TestData[,factors[i]]))==F) {
                missinglevels <- T
                cat(paste("\n", "WARNING: differences in factor levels between calibration and evaluation data (variable ", factors[i], ")", "\n", sep = ""))
                cat(paste("same levels set for both data sets to avoid problems with evaluations for RF and SVM","\n", sep = ""))
                uniquelevels <- unique(c(levels(TrainData[,factors[i]]), levels(TestData[,factors[i]])))
                levels(TrainData[,factors[i]]) <- uniquelevels
                levels(TestData[,factors[i]]) <- uniquelevels
            } 
        }
    }
    assign("TestData", TestData, envir=.BiodiversityR)
    assign("TrainData", TrainData, envir=.BiodiversityR)
    if (VIF == T) {
        if (! require(car)) {stop("Please install the car package")}
        GLM2.formula <- ensemble.formulae(x, factors=factors)$RF.formula
        assign("GLM2.formula", GLM2.formula, envir=.BiodiversityR)
        GLM2.family=binomial(link="logit")
        assign("GLM2.family", GLM2.family, envir=.BiodiversityR)
        vifresult <- NULL
        tryCatch(vifresult <- vif(glm(formula=GLM2.formula, family=GLM2.family, data=TrainData, control=glm.control(maxit=maxit))),
            error= function(err) {print(paste("VIF evaluation failed"))},
                    silent=T)
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
        corresult <- round(100*corresult, digits=1)
        cat(paste("\n", "Correlation between numeric variables (as percentage)", "\n", sep = ""))        
        print(corresult)
    }
    if(evaluations.keep==T) {
        evaluations <- list(p=p, a=a, pt=pt, at=at, MAXENT.C=NULL, MAXENT.T=NULL, 
            GBM.trees=NULL, GBM.C=NULL, GBM.T=NULL, GBMSTEP.trees=NULL, GBMSTEP.C=NULL, GBMSTEP.T=NULL, 
            RF.C=NULL, RF.T=NULL, GLM.C=NULL, GLM.T=NULL, GLMS.C=NULL, GLMS.T=NULL, 
            GAM.C=NULL, GAM.T=NULL, GAMS.C=NULL, GAMS.T=NULL, MGCV.C=NULL, MGCV.T=NULL, MGCVF.C=NULL, MGCVF.T=NULL,
            EARTH.C=NULL, EARTH.T=NULL, RPART.C=NULL, RPART.T=NULL,
            NNET.C=NULL, NNET.T=NULL, FDA.C=NULL, FDA.T=NULL, SVM.C=NULL, SVM.T=NULL, SVME.C=NULL, SVME.T=NULL,
            BIOCLIM.C=NULL, BIOCLIM.T=NULL, DOMAIN.C=NULL, DOMAIN.T=NULL, MAHAL.C=NULL, MAHAL.T=NULL,
            GEODIST.C=NULL, GEODIST.T=NULL)
    }
    if(models.keep==T) {
        models <- list(MAXENT=NULL, GBM=NULL, GBMSTEP=NULL, RF=NULL, GLM=NULL, 
            GLMSTEP=NULL, GAM=NULL, GAMSTEP=NULL, MGCV=NULL, MGCVFIX=NULL, EARTH=NULL, RPART=NULL, 
            NNET=NULL, FDA=NULL, SVM=NULL, SVME=NULL, BIOCLIM=NULL, DOMAIN=NULL, MAHAL=NULL, GEODIST=NULL)
    }
    Yweights1 <- Yweights
    if (Yweights == "BIOMOD") {
        #have equal weight of presence vs. background
        Yweights1 <- numeric(length = nrow(TrainData))
        pres <- nrow(p)
        abs <- nrow(a)
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
        eval1 <- eval2 <- results <- NULL
        # Put the file 'maxent.jar' in the 'java' folder of dismo
        # the file 'maxent.jar' can be obtained from from http://www.cs.princeton.edu/~schapire/maxent/.
        jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
        tryCatch(results <- maxent(x=x, p=p, a=a, factors=factors),
            error= function(err) {print(paste("MAXENT calibration failed"))},
            silent=T)
        if (is.null(results) == F) {
            cat(paste("\n", "MAXENT calibration","\n\n",sep = ""))
            eval1 <- evaluate(p=p, a=a, model=results, x=x)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "MAXENT evaluation","\n\n",sep = ""))
                tryCatch(eval2 <- evaluate(p=pt, a=at, model=results, x=x),
                    error= function(err) {print(paste("MAXENT evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
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
        eval1 <- eval2 <- results <- NULL
        tryCatch(results <- gbm(formula=GBM.formula, data=TrainData, weights=Yweights, distribution="bernoulli", 
            interaction.depth=7, shrinkage=0.001, bag.fraction=0.5, train.fraction=1, 
            n.trees=GBM.n.trees, verbose=F, cv.folds=5),
            error= function(err) {print(paste("GBM calibration failed"))},
            silent=T)
        if (is.null(results) == F) {
            cat(paste("\n", "GBM calibration","\n\n",sep = ""))
            eval1 <- evaluate(p=p, a=a, model=results, x=x, n.trees=results$n.trees, type="response")
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "GBM evaluation","\n\n",sep = ""))
                tryCatch(eval2 <- evaluate(p=pt, a=at, model=results, x=x, n.trees=results$n.trees, type="response"),
                    error= function(err) {print(paste("GBM evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
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
            }
        }else{ cat(paste("\n", "WARNING: MAXENT calibration failed", "\n", "\n"))}
#       detach(package:gbm) 
    }
    if (GBMSTEP > 0) {
#        require(gbm, quietly=T)        
        eval1 <- eval2 <- results <- NULL
        tryCatch(results <- gbm.step(data=TrainData, gbm.y=1, gbm.x=GBMSTEP.gbm.x, family="bernoulli",
            site.weights=Yweights, tree.complexity=GBMSTEP.tree.complexity, learning.rate = GBMSTEP.learning.rate, 
            bag.fraction=GBMSTEP.bag.fraction, step.size=GBMSTEP.step.size, verbose=F, silent=T, plot.main=F),
            error= function(err) {print(paste("stepwise GBM calibration failed"))},
            silent=T)
        if (is.null(results) == F) {
            cat(paste("\n", "stepwise GBM trees (target > 1000)", "\n", sep = ""))
            print(results$n.trees)
            cat(paste("\n", "stepwise GBM calibration","\n\n",sep = ""))
            eval1 <- evaluate(p=p, a=a, model=results, x=x, n.trees=results$n.trees, type="response")
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "stepwise GBM evaluation","\n\n",sep = ""))
                tryCatch(eval2 <- evaluate(p=pt, a=at, model=results, x=x, n.trees=results$n.trees, type="response"),
                    error= function(err) {print(paste("stepwise GBM evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
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
        eval1 <- eval2 <- results <- TestPres <- TestAbs <- NULL
        tryCatch(results <- randomForest(formula=RF.formula, ntree=RF.ntree, mtry=RF.mtry, data=TrainData, na.action=na.omit),
            error= function(err) {print(paste("random forest calibration failed"))},
            silent=T)
        if (is.null(results) == F) {
            cat(paste("\n", "RF calibration","\n\n",sep = ""))
#           treat factors as factors during evaluations
            TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="response")
            TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="response")
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "RF evaluation","\n\n",sep = ""))
                tryCatch(TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="response"),
                    error= function(err) {print(paste("RF evaluation (step 1) failed"))},
                    silent=T)
                tryCatch(TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="response"),
                    error= function(err) {print(paste("RF evaluation (step 2) failed"))},
                    silent=T)
                if (is.null(TestPres) == F && is.null(TestAbs) == F) {
                    eval2 <-  evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
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
            }
        }else{ cat(paste("\n", "WARNING: random forest calibration failed", "\n", "\n"))}
#       detach(package:randomForest)
    }
    if (GLM > 0) {
        eval1 <- eval2 <- results <- TestPres <- TestAbs <- NULL
        tryCatch(results <- glm(formula=GLM.formula, family=GLM.family, data=TrainData, weights=Yweights, control=glm.control(maxit=maxit)),
            error= function(err) {print(paste("GLM calibration failed"))},
            silent=T)
        if (is.null(results) == F) {
            cat(paste("\n", "GLM calibration","\n\n",sep = ""))
#           treat factors as factors during evaluations
            TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="response")
            TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="response")
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "GLM evaluation","\n\n",sep = ""))
                tryCatch(TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="response"),
                    error= function(err) {print(paste("GLM evaluation (step 1) failed"))},
                    silent=T)
                tryCatch(TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="response"),
                    error= function(err) {print(paste("GLM evaluation (step 2) failed"))},
                    silent=T)
                if (is.null(TestPres) == F && is.null(TestAbs) == F) {
                    eval2 <-  evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
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
            }
        }else{ cat(paste("\n", "WARNING: GLM calibration failed", "\n", "\n"))}
    }
    if (GLMSTEP > 0) {
#        require(MASS, quietly=T)
        eval1 <- eval2 <- results <- results2 <- TestPres <- TestAbs <- NULL
        tryCatch(results <- glm(formula=STEP.formula, family=GLM.family, data=TrainData, weights=Yweights, control=glm.control(maxit=maxit)),
            error= function(err) {print(paste("first step of stepwise GLM calibration failed"))},
            silent=T)
        tryCatch(results2 <- stepAIC(results, scope=GLMSTEP.scope, direction="both", trace=F, steps=GLMSTEP.steps, k=GLMSTEP.k),
            error= function(err) {print(paste("stepwise GLM calibration failed"))},
            silent=T)
        if (is.null(results2) == F) {
            results <- results2
            cat(paste("\n", "stepwise GLM formula","\n\n",sep = ""))
            print(formula(results))
            cat(paste("\n", "stepwise GLM calibration","\n\n",sep = ""))
            TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="response")
            TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="response")
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "stepwise GLM evaluation","\n\n",sep = ""))
                tryCatch(TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="response"),
                    error= function(err) {print(paste("stepwise GLM evaluation (step 1) failed"))},
                    silent=T)
                tryCatch(TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="response"),
                    error= function(err) {print(paste("stepwise GLM evaluation (step 2) failed"))},
                    silent=T)
                if (is.null(TestPres) == F && is.null(TestAbs) == F) {
                    eval2 <-  evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
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
            }
        }else{ cat(paste("\n", "WARNING: stepwise GBM calibration failed", "\n", "\n"))}
#       detach(package:MASS)
    }
    if (GAM > 0) {
        require(gam, quietly=T)
        eval1 <- eval2 <- results <- TestPres <- TestAbs <- NULL
        tryCatch(results <- gam(formula=GAM.formula, family=GAM.family, data=TrainData, weights=Yweights, control=gam.control(maxit=maxit, bf.maxit=50)),
            error= function(err) {print(paste("GAM calibration (gam package) failed"))},
            silent=T)
        if (is.null(results) == F) {
            cat(paste("\n", "GAM calibration (gam package)","\n\n",sep = ""))
            TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="response")
            TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="response")
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "GAM evaluation (gam package)","\n\n",sep = ""))
                tryCatch(TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="response"),
                    error= function(err) {print(paste("GAM evaluation (step 1, gam package) failed"))},
                    silent=T)
                tryCatch(TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="response"),
                    error= function(err) {print(paste("GAM evaluation (step 2, gam package) failed"))},
                    silent=T)
                if (is.null(TestPres) == F && is.null(TestAbs) == F) {
                    eval2 <-  evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
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
            }
        }else{ cat(paste("\n", "WARNING: GAM calibration (gam package) failed", "\n", "\n"))}
# conflict with mgcv gam function
        detach(package:gam)
    }
    if (GAMSTEP > 0) {
        require(gam, quietly=T)
        eval1 <- eval2 <- results <- results2 <- TestPres <- TestAbs <- NULL
        tryCatch(results <- gam(formula=STEP.formula, family=GAM.family, data=TrainData, weights=Yweights, control=gam.control(maxit=maxit, bf.maxit=50)), 
            error= function(err) {print(paste("first step of stepwise GAM calibration (gam package) failed"))},
            silent=T)
        assign("TrainData", TrainData, pos=GAMSTEP.pos)
        assign("GAM.family", GAM.family, pos=GAMSTEP.pos)
        assign("maxit", maxit, pos=GAMSTEP.pos)      
        tryCatch(results2 <- step.gam(results, scope=GAMSTEP.scope, direction="both", trace=F, steps=GAMSTEP.steps), 
            error= function(err) {print(paste("stepwise GAM calibration (gam package) failed"))},
            silent=F)
        remove(TrainData, pos=GAMSTEP.pos)
        remove(GAM.family, pos=GAMSTEP.pos)
        remove(maxit, pos=GAMSTEP.pos)
        if (is.null(results2) == F) {
            results <- results2
            cat(paste("\n", "stepwise GAM formula (gam package)","\n\n",sep = ""))
            print(formula(results))
            cat(paste("\n", "stepwise GAM calibration (gam package)","\n\n",sep = ""))
            TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="response")
            TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="response")
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "stepwise GAM evaluation (gam package)","\n\n","\n", sep = ""))
                tryCatch(TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="response"),
                    error= function(err) {print(paste("stepwise GAM evaluation (step 1, gam package) failed"))},
                    silent=T)
                tryCatch(TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="response"),
                    error= function(err) {print(paste("stepwise GAM evaluation (step 2, gam package) failed"))},
                    silent=T)
                if (is.null(TestPres) == F && is.null(TestAbs) == F) {
                    eval2 <-  evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
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
            }
        }else{ cat(paste("\n", "WARNING: stepwise GAM calibration (gam package) failed", "\n", "\n"))}
# conflict with mgcv gam function
        detach(package:gam)
    }
    if (MGCV > 0) {
        require(mgcv, quietly=T)
        eval1 <- eval2 <- results <- NULL
        tryCatch(results <- gam(formula=MGCV.formula, family=GAM.family, data=TrainData, weights=Yweights, na.action=na.omit, 
                select=MGCV.select, control=gam.control(maxit=maxit)),
            error= function(err) {print(paste("GAM calibration (mgcv package) failed"))},
            silent=T)
        if (is.null(results) == F) {
            cat(paste("\n", "GAM calibration (mgcv package)","\n\n",sep = ""))
            eval1 <- evaluate(p=p, a=a, model=results, x=x, type="response", na.action=na.omit)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "GAM evaluation (mgcv package)","\n\n",sep = ""))
                tryCatch(eval2 <- evaluate(p=pt, a=at, model=results, x=x, type="response", na.action=na.omit),
                    error= function(err) {print(paste("GAM evaluation (mgcv package) failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
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
            }
        }else{ cat(paste("\n", "WARNING: GAM calibration (mgcv package) failed", "\n", "\n"))}
# conflict with gam package
        detach(package:mgcv)
    }
    if (MGCVFIX > 0) {
        require(mgcv, quietly=T)
        eval1 <- eval2 <- results <- NULL
        tryCatch(results <- gam(formula=MGCVFIX.formula, family=GAM.family, data=TrainData, weights=Yweights, na.action=na.omit, select=FALSE, control=gam.control(maxit=maxit)),
            error= function(err) {print(paste("MGCVFIX calibration (mgcv package) failed"))},
            silent=T)
        if (is.null(results) == F) {
            cat(paste("\n", "GAM with fixed d.f. regression splines calibration (mgcv package)","\n\n",sep = ""))
            eval1 <- evaluate(p=p, a=a, model=results, x=x, type="response", na.action=na.omit)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "GAM with fixed d.f. regression splines evaluation (mgcv package)","\n\n",sep = ""))
                tryCatch(eval2 <- evaluate(p=pt, a=at, model=results, x=x, type="response", na.action=na.omit),
                    error= function(err) {print(paste("GAMFIX evaluation (mgcv package) failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
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
        eval1 <- eval2 <- results <- NULL
        tryCatch(results <- earth(formula=EARTH.formula, glm=EARTH.glm, data=TrainData, degree=2),
            error= function(err) {print(paste("MARS calibration (earth package) failed"))},
            silent=T)
        if (is.null(results) == F) {
            cat(paste("\n", "MARS calibration (earth package)","\n\n",sep = ""))
            eval1 <- evaluate(p=p, a=a, model=results, x=x, type="response")
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "MARS evaluation (earth package)","\n\n",sep = ""))
                tryCatch(eval2 <- evaluate(p=pt, a=at, model=results, x=x, type="response"),
                    error= function(err) {print(paste("MARS evaluation (earth package) failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
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
            }
        }else{ cat(paste("\n", "WARNING: MARS calibration (earth package) failed", "\n", "\n"))}
#       detach(package:earth)
    }
    if (RPART > 0) {
#        require(rpart, quietly=T)
        eval1 <- eval2 <- results <- TestPres <- TestAbs <- NULL
        tryCatch(results <- rpart(formula=RPART.formula, data=TrainData, weights=Yweights,
            control=rpart.control(xval=RPART.xval, minbucket=5, minsplit=5, cp=0.001, maxdepth=25)),
            error= function(err) {print(paste("RPART calibration failed"))},
            silent=T)
        assign("results", results, envir=.BiodiversityR)
        if (is.null(results) == F) {
            cat(paste("\n", "RPART calibration","\n\n",sep = "")) 
            TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="prob")[,2]
            TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="prob")[,2]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "RPART evaluation","\n\n",sep = ""))
                tryCatch(TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="prob")[,2],
                    error= function(err) {print(paste("RPART evaluation (step 1) failed"))},
                    silent=T)
                tryCatch(TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="prob")[,2],
                    error= function(err) {print(paste("RPART evaluation (step 2) failed"))},
                    silent=T)
                if (is.null(TestPres) == F && is.null(TestAbs) == F) {
                    eval2 <-  evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
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
            }
        }else{ cat(paste("\n", "WARNING: RPART calibration failed", "\n", "\n"))}
#       detach(package:rpart)
    }
    if (NNET > 0) {
#        require(nnet, quietly=T)
        eval1 <- eval2 <- results <- NULL
        tryCatch(results <- nnet(formula=NNET.formula, size=NNET.size, decay=NNET.decay, data=TrainData, weights=Yweights, 
            rang=0.1, maxit=maxit, trace=F),
            error= function(err) {print(paste("ANN calibration (nnet package) failed"))},
            silent=T)
        if (is.null(results) == F) {
            cat(paste("\n", "ANN calibration (nnet package)","\n\n",sep = ""))
            eval1 <- evaluate(p=p, a=a, model=results, x=x, type="raw")
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "ANN evaluation (nnet package)","\n\n",sep = ""))
                tryCatch(eval2 <- evaluate(p=pt, a=at, model=results, x=x, type="raw"),
                    error= function(err) {print(paste("ANN evaluation (nnet package) failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
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
            }
        }else{ cat(paste("\n", "WARNING: ANN calibration (nnet package) failed", "\n", "\n"))}
#       detach(package:nnet)
    }
    if (FDA > 0) {
#        require(mda, quietly=T)
        eval1 <- eval2 <- results <- TestPres <- TestAbs <- NULL
        tryCatch(results <- fda(formula=FDA.formula, method=mars, data=TrainData, weights=Yweights),
            error= function(err) {print(paste("FDA calibration failed"))},
            silent=T)
        if (is.null(results) == F) {
            cat(paste("\n", "FDA calibration","\n\n",sep = "")) 
            TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="posterior")[,2]
            TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="posterior")[,2]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "FDA evaluation","\n\n",sep = ""))
                tryCatch(TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="posterior")[,2],
                    error= function(err) {print(paste("FDA evaluation (step 1) failed"))},
                    silent=T)
                tryCatch(TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="posterior")[,2],
                    error= function(err) {print(paste("FDA evaluation (step 2) failed"))},
                    silent=T)
                if (is.null(TestPres) == F && is.null(TestAbs) == F) {
                    eval2 <-  evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
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
            }
        }else{ cat(paste("\n", "WARNING: FDA calibration failed", "\n", "\n"))}
#       detach(package:mda)
    }
    if (SVM > 0) {
#        require(kernlab, quietly=T)
        eval1 <- eval2 <- results <- TestPres <- TestAbs <- NULL
        tryCatch(results <- ksvm(SVM.formula, data=TrainData, type="C-svc", prob.model=T),
            error= function(err) {print(paste("SVM calibration failed"))},
            silent=T)
        if (is.null(results) == F) {
            cat(paste("\n", "SVM calibration (kernlab package)","\n\n",sep = ""))
            TrainPres <- kernlab::predict(results, newdata=TrainData[TrainData[,"pb"]==1,], type="probabilities")[,2]
            TrainAbs <- kernlab::predict(results, newdata=TrainData[TrainData[,"pb"]==0,], type="probabilities")[,2]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "SVM evaluation (kernlab package)","\n\n",sep = ""))
                tryCatch(TestPres <- kernlab::predict(results, newdata=TestData[TestData[,"pb"]==1,], type="probabilities")[,2],
                    error= function(err) {print(paste("SVM evaluation (step 1, kernlab package) failed"))},
                    silent=T)
                tryCatch(TestAbs <- kernlab::predict(results, newdata=TestData[TestData[,"pb"]==0,], type="probabilities")[,2],
                    error= function(err) {print(paste("SVM evaluation (step 2, kernlab package) failed"))},
                    silent=T)
                if (is.null(TestPres) == F && is.null(TestAbs) == F) {
                    eval2 <-  evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
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
            }
        }else{ cat(paste("\n", "WARNING: SVM calibration failed", "\n", "\n"))}
#       detach(package:kernlab)
    }
    if (SVME > 0) {
#        require(e1071, quietly=T)
        eval1 <- eval2 <- results <- TestPres <- TestAbs <- NULL
        tryCatch(results <- svm(SVME.formula, data=TrainData, type="C-classification", kernel="polynomial", degree=3, probability=T),
            error= function(err) {print(paste("SVM calibration (e1071 package) failed"))},
            silent=T)
        if (is.null(results) == F) {
            cat(paste("\n", "SVM calibration (e1071 package)","\n\n",sep = ""))
            TrainPres <- predict.svme(results, newdata=TrainData[TrainData[,"pb"]==1,])
            TrainAbs <- predict.svme(results, newdata=TrainData[TrainData[,"pb"]==0,])
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "SVM evaluation (e1071 package)","\n\n",sep = ""))
                tryCatch(TestPres <- predict.svme(results, newdata=TestData[TestData[,"pb"]==1,]),
                    error= function(err) {print(paste("SVM evaluation (step 1, e1071 package) failed"))},
                    silent=T)
                tryCatch(TestAbs <- predict.svme(results, newdata=TestData[TestData[,"pb"]==0,]),
                    error= function(err) {print(paste("SVM evaluation (step 2, e1071 package) failed"))},
                    silent=T)
                if (is.null(TestPres) == F && is.null(TestAbs) == F) {
                    eval2 <-  evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
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
            }
        }else{ cat(paste("\n", "WARNING: SVM calibration (e1071 package) failed", "\n", "\n"))}
#       detach(package:e1071)
    }
    if (BIOCLIM > 0 || DOMAIN > 0 || MAHAL > 0) {
        if(is.null(factors)==F) {
            for (i in 1:length(factors)) {
                x <- dropLayer(x, which(names(x) == factors[i]))              
            }
        }
    }
    if (BIOCLIM > 0) {
        eval1 <- eval2 <- results <- NULL
        tryCatch(results <- bioclim(x=x, p=p),
            error= function(err) {print(paste("BIOCLIM calibration failed"))},
            silent=T)
        if (is.null(results) == F) {
            cat(paste("\n", "BIOCLIM calibration","\n\n",sep = ""))
            eval1 <- evaluate(p=p, a=a, model=results, x=x)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "BIOCLIM evaluation","\n\n",sep = ""))
                tryCatch(eval2 <- evaluate(p=pt, a=at, model=results, x=x),
                    error= function(err) {print(paste("BIOCLIM evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
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
        eval1 <- eval2 <- results <- NULL
        tryCatch(results <- domain(x=x, p=p, factors=factors),
            error= function(err) {print(paste("DOMAIN calibration failed"))},
            silent=T)
        if (is.null(results) == F) {
            cat(paste("\n", "DOMAIN calibration","\n\n",sep = ""))
            eval1 <- evaluate(p=p, a=a, model=results, x=x)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "DOMAIN evaluation","\n\n",sep = ""))
                tryCatch(eval2 <- evaluate(p=pt, a=at, model=results, x=x),
                    error= function(err) {print(paste("DOMAIN evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
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
        eval1 <- eval2 <- results <- NULL
        tryCatch(results <- mahal(x=x, p=p),
            error= function(err) {print(paste("Mahalanobis calibration failed"))},
            silent=T)
        if (is.null(results) == F) {
            cat(paste("\n", "Mahalanobis calibration","\n\n",sep = ""))
            eval1 <- evaluate(p=p, a=a, model=results, x=x)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "Mahalanobis evaluation","\n\n",sep = ""))
                tryCatch(eval2 <- evaluate(p=pt, a=at, model=results, x=x),
                    error= function(err) {print(paste("Mahalanobis evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
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
        eval1 <- eval2 <- results <- NULL
        tryCatch(results <- geoDist(p=p, lonlat=TRUE),
            error= function(err) {print(paste("GEODIST calibration failed"))},
            silent=T)
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
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "GEODIST evaluation","\n\n","\n", sep = ""))
                prest_geo <- extract(pgeo, pt)
                abst_geo <- extract(pgeo, at)
                tryCatch(eval2 <- evaluate(p=prest_geo, a=abst_geo, tr=quantile(pgeo, na.rm=T, probs=c(1:50/50), names=F)),
                    error= function(err) {print(paste("GEODIST evaluation failed"))},
                    silent=T)
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
    cat(paste("\n", "End of evaluations", "\n\n",sep = "")) 
    remove(Yweights, envir=.BiodiversityR)
    remove(TrainData, envir=.BiodiversityR)
    remove(TestData, envir=.BiodiversityR)
    if (evaluations.keep==T && models.keep==F) {return(evaluations)}
    if (evaluations.keep==F && models.keep==T) {return(models)}
    if (evaluations.keep==T && models.keep==T) {return(list(evaluations=evaluations, models=models))}
}




