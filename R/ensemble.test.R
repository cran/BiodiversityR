`ensemble.test` <- function(
    x, p, a=NULL, an=1000, ext=NULL, k=0, pt=NULL, at=NULL,
    PLOTS=TRUE, evaluations.keep=FALSE,
    MAXENT=1, GBM=1, GBMSTEP=1, RF=1, GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=1, 
    MGCV=1, EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, BIOCLIM=1, DOMAIN=1, MAHAL=1,   
    Yweights="BIOMOD", factors=NULL,
    formulae.defaults=TRUE,
    GBM.formula=NULL, GBM.n.trees=3000,
    GBMSTEP.gbm.x=2:(1+nlayers(x)), GBMSTEP.tree.complexity=5, GBMSTEP.learning.rate=0.01, GBMSTEP.bag.fraction=0.5,
    RF.formula=NULL, RF.ntree=750, RF.mtry=log(nlayers(x)), 
    GLM.formula=NULL, GLM.family=binomial(link="logit"), 
    GLMSTEP.steps=1000, STEP.formula=NULL, GLMSTEP.scope=NULL, GLMSTEP.k=2,
    GAM.formula=NULL, GAM.family=binomial(link="logit"), 
    GAMSTEP.steps=1000, GAMSTEP.scope=NULL,
    MGCV.formula=NULL,
    EARTH.formula=NULL, EARTH.glm=list(family=binomial(link="logit"), maxit=1000),
    RPART.formula=NULL, RPART.xval=50, 
    NNET.formula=NULL, NNET.size=8, NNET.decay=0.01,
    FDA.formula=NULL,
    SVM.formula=NULL
)
{
    .BiodiversityR <- new.env()
    if (! require(dismo)) {stop("Please install the dismo package")}
    if (MAXENT > 0) {
	jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
	if (!file.exists(jar)) {stop('maxent program is missing: ', jar, '\nPlease download it here: http://www.cs.princeton.edu/~schapire/maxent/')}
    }
    if (formulae.defaults == T) {
        formulae <- ensemble.formulae(x, factors=factors)
    }
    if (GBM > 0) {
        if (! require(gbm)) {stop("Please install the gbm package")}
        if (is.null(GBM.formula) == T && formulae.defaults == T) {GBM.formula <- formulae$GBM.formula}
        if (is.null(GBM.formula) == T) {stop("Please provide the GBM.formula (hint: use ensemble.formulae function)")}
        environment(GBM.formula) <- .BiodiversityR
    }
    if (GBMSTEP > 0) {if (! require(gbm)) {stop("Please install the gbm package")}}
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
    if (is.null(a)==T) {a <- randomPoints(x, n=an, ext=ext)}
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
#    if(is.null(factors)==F) {
#        for (i in 1:length(factors)) {
#            j <- which(names(x) == factors[i])
#            x[[j]] <- as.factor(x[[j]])
#        }
#    }
    TrainData <- prepareData(x, p, b=a, factors=factors, xy=FALSE)
    if(any(is.na(TrainData[TrainData[,"pb"]==1,]))) {
        cat(paste("\n", "WARNING: presence locations with missing data removed from calibration data","\n", "\n", sep = ""))
    }
    TrainValid <- complete.cases(TrainData[TrainData[,"pb"]==1,])
    p <- p[TrainValid,]
    if(any(is.na(TrainData[TrainData[,"pb"]==0,]))) {
        cat(paste("\n", "WARNING: background locations with missing data removed from calibration data","\n", "\n", sep = ""))
    }
    TrainValid <- complete.cases(TrainData[TrainData[,"pb"]==0,])
    a <- a[TrainValid,]
    TrainData <- prepareData(x, p, b=a, factors=factors, xy=FALSE)
    assign("TrainData", TrainData, envir=.BiodiversityR)
    TestData <- prepareData(x, p=pt, b=at, factors=factors, xy=FALSE)
    if(any(is.na(TestData[TestData[,"pb"]==1,]))) {
        cat(paste("\n", "WARNING: presence locations with missing data removed from evaluation data","\n", "\n", sep = ""))
    }
    TestValid <- complete.cases(TestData[TestData[,"pb"]==1,])
    pt <- pt[TestValid,]
    if(any(is.na(TestData[TestData[,"pb"]==0,]))) {
        cat(paste("\n", "WARNING: background locations with missing data removed from evaluation data","\n", "\n", sep = ""))
    }
    TestValid <- complete.cases(TestData[TestData[,"pb"]==0,])
    at <- at[TestValid,]
    TestData <- prepareData(x, p=pt, b=at, factors=factors, xy=FALSE)
    assign("TestData", TestData, envir=.BiodiversityR)
    missinglevels <- F
    if(is.null(factors)==F) {
        for (i in 1:length(factors)) {
            if (identical(levels(TrainData[,factors[i]]),levels(TestData[,factors[i]]))==F) {
                missinglevels <- T
                cat(paste("\n", "WARNING: differences in factor levels between calibration and evaluation data", "\n", sep = ""))
                cat(paste("random forest and support vector machine evaluations not done for this reason","\n", sep = ""))
            }                 
        }
    }
    if(evaluations.keep==T) {
        evaluations <- list(MAXENT.C=NULL, MAXENT.T=NULL, GBM.trees=NULL, GBM.C=NULL, GBM.T=NULL,
            GBMSTEP.trees=NULL, GBMSTEP.C=NULL, GBMSTEP.T=NULL, 
            RF.C=NULL, RF.T=NULL, GLM.C=NULL, GLM.T=NULL, GLMS.C=NULL, GLMS.T=NULL, 
            GAM.C=NULL, GAM.T=NULL, GAMS.C=NULL, GAMS.T=NULL,
            MGCV.C=NULL, MGCV.T=NULL, EARTH.C=NULL, EARTH.T=NULL, RPART.C=NULL, RPART.T=NULL,
            NNET.C=NULL, NNET.T=NULL, FDA.C=NULL, FDA.T=NULL, SVM.C=NULL, SVM.T=NULL,
            BIOCLIM.C=NULL, BIOCLIM.T=NULL, DOMAIN.C=NULL, DOMAIN.T=NULL, MAHAL.C=NULL, MAHAL.T=NULL)
    }
    results <- NULL
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
    if (MAXENT > 0) {
        eval1 <- eval2 <- results <- NULL
        # Put the file 'maxent.jar' in the 'java' folder of dismo
        # the file 'maxent.jar' can be obtained from from http://www.cs.princeton.edu/~schapire/maxent/.
        jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
        results <- maxent(x=x, p=p, a=a, factors=factors)
        cat(paste("\n", "MAXENT calibration","\n", "\n", sep = ""))
        eval1 <- evaluate(p=p, a=a, model=results, x=x)
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "MAXENT evaluation","\n", "\n", "\n", sep = ""))
            eval2 <- evaluate(p=pt, a=at, model=results, x=x)
            print(eval2)
            if(PLOTS==T) {
                plot(eval2, "ROC") 
                title(main="MAXENT", cex=2, adj=0, col.main="blue")
            }
        }
        if(evaluations.keep ==T) {
            evaluations$MAXENT.C <- eval1
            evaluations$MAXENT.T <- eval2
        }
    }
    if (GBM > 0) {
#        require(gbm, quietly=T)        
        eval1 <- eval2 <- results <- NULL
        results <- gbm(formula=GBM.formula, data=TrainData, weights=Yweights, distribution = "bernoulli", 
            interaction.depth = 7, shrinkage = 0.001, bag.fraction = 0.5, train.fraction = 1, 
            n.trees = GBM.n.trees, verbose = F, cv.folds = 5)
        cat(paste("\n", "GBM calibration","\n", "\n", sep = ""))
        eval1 <- evaluate(p=p, a=a, model=results, x=x, n.trees=results$n.trees, type="response")
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "GBM evaluation","\n", "\n", sep = ""))
            eval2 <- evaluate(p=pt, a=at, model=results, x=x, n.trees=results$n.trees, type="response")
            print(eval2)
            if(PLOTS==T) {
                plot(eval2, "ROC") 
                title(main="BRT (gbm)", cex=2, adj=0, col.main="blue")
            }
        }
#       detach(package:gbm) 
        if(evaluations.keep ==T) {
            evaluations$GBM.trees <- results$n.trees 
            evaluations$GBM.C <- eval1
            evaluations$GBM.T <- eval2
        }
    }
    if (GBMSTEP > 0) {
#        require(gbm, quietly=T)        
        eval1 <- eval2 <- results <- NULL
        results <- gbm.step(data=TrainData, gbm.y=1, gbm.x=GBMSTEP.gbm.x, family="bernoulli",
            site.weights=Yweights, tree.complexity=GBMSTEP.tree.complexity, learning.rate = GBMSTEP.learning.rate, 
            bag.fraction=GBMSTEP.bag.fraction, verbose=F, silent=T, plot.main=F)
        cat(paste("\n", "stepwise GBM trees (target > 1000)", "\n", sep = ""))
        print(results$n.trees)
        cat(paste("\n", "stepwise GBM calibration","\n", "\n", sep = ""))
        eval1 <- evaluate(p=p, a=a, model=results, x=x, n.trees=results$n.trees, type="response")
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "stepwise GBM evaluation","\n", "\n", sep = ""))
            eval2 <- evaluate(p=pt, a=at, model=results, x=x, n.trees=results$n.trees, type="response")
            print(eval2)
            if(PLOTS==T) {
                plot(eval2, "ROC") 
                title(main="BRT (gbm.step)", cex=2, adj=0, col.main="blue")
            }
        }
#       detach(package:gbm) 
        if(evaluations.keep ==T) {
            evaluations$GBMSTEP.trees <- results$n.trees 
            evaluations$GBMSTEP.C <- eval1
            evaluations$GBMSTEP.T <- eval2
        }
    }
    if (RF > 0) {
#        require(randomForest, quietly=T)
        eval1 <- eval2 <- results <- NULL
        results <- randomForest(formula=RF.formula, ntree=RF.ntree, mtry=RF.mtry, data=TrainData, na.action=na.omit)
        cat(paste("\n", "RF calibration","\n", "\n", sep = ""))
#       treat factors as factors during evaluations
        TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="response")
        TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="response")
        eval1 <- evaluate(p=TrainPres, a=TrainAbs)
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            if (missinglevels==T) {
                cat(paste("\n", "random forest evaluations not done as factor levels are different","\n", "\n", sep = ""))
            }else{
                cat(paste("\n", "RF evaluation","\n", "\n", sep = ""))
                TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="response")
                TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="response")
                eval2 <-  evaluate(p=TestPres, a=TestAbs)
                print(eval2)
                if(PLOTS==T) {
                    plot(eval2, "ROC") 
                    title(main="RF", cex=2, adj=0, col.main="blue")
                }
            }
        }
#       detach(package:randomForest)
        if(evaluations.keep ==T) {
            evaluations$RF.C <- eval1
            evaluations$RF.T <- eval2
        }
    }
    if (GLM > 0) {
        eval1 <- eval2 <- results <- NULL
        results <- glm(formula=GLM.formula, family=GLM.family, data=TrainData, weights=Yweights, maxit=100)
        cat(paste("\n", "GLM calibration","\n", "\n", sep = ""))
#       treat factors as factors during evaluations
        TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="response")
        TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="response")
        eval1 <- evaluate(p=TrainPres, a=TrainAbs)
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "GLM evaluation","\n", "\n", sep = ""))
            TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="response")
            TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="response")
            eval2 <-  evaluate(p=TestPres, a=TestAbs)
            print(eval2)
            if(PLOTS==T) {
                plot(eval2, "ROC") 
                title(main="GLM", cex=2, adj=0, col.main="blue")
            }
        }
        if(evaluations.keep ==T) {
            evaluations$GLM.C <- eval1
            evaluations$GLM.T <- eval2
        }
    }
    if (GLMSTEP > 0) {
#        require(MASS, quietly=T)
        eval1 <- eval2 <- results <- NULL
        results <- glm(formula=STEP.formula, family=GLM.family, data=TrainData, weights=Yweights, maxit=100)
        results <- stepAIC(results, scope=GLMSTEP.scope, direction="both", trace=F, steps=GLMSTEP.steps, k=GLMSTEP.k, maxit=100)
        cat(paste("\n", "stepwise GLM calibration","\n", "\n", sep = ""))
        TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="response")
        TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="response")
        eval1 <- evaluate(p=TrainPres, a=TrainAbs)
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "stepwise GLM evaluation","\n", "\n", sep = ""))
            TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="response")
            TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="response")
            eval2 <-  evaluate(p=TestPres, a=TestAbs)
            print(eval2)
            if(PLOTS==T) {
                plot(eval2, "ROC") 
                title(main="STEPWISE GLM", cex=2, adj=0, col.main="blue")
            }
        }
#       detach(package:MASS)
        if(evaluations.keep ==T) {
            evaluations$GLMS.C <- eval1
            evaluations$GLMS.T <- eval2
        }
    }
    if (GAM > 0) {
        require(gam, quietly=T)
        eval1 <- eval2 <- results <- NULL
        results <- gam(formula=GAM.formula, family=GAM.family, data=TrainData, weights=Yweights, maxit=50, bf.maxit=50)
        cat(paste("\n", "GAM calibration (gam package)","\n", "\n", sep = ""))
        TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="response")
        TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="response")
        eval1 <- evaluate(p=TrainPres, a=TrainAbs)
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "GAM evaluation (gam package)","\n", "\n", sep = ""))
            TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="response")
            TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="response")
            eval2 <-  evaluate(p=TestPres, a=TestAbs)
            print(eval2)
            if(PLOTS==T) {
                plot(eval2, "ROC") 
                title(main="GAM (gam)", cex=2, adj=0, col.main="blue")
            }
        }
# conflict with mgcv gam function
        detach(package:gam)
        if(evaluations.keep ==T) {
            evaluations$GAM.C <- eval1
            evaluations$GAM.T <- eval2
        }
    }

    if (GAMSTEP > 0) {
        cat(paste("\n", "NOTE: stepwise GAM seems to require assignments of family and data to pos 1", "\n", "\n", sep=""))
        GAMSTEP <- 0
    } 
    if (GAMSTEP > 0) {
        require(gam, quietly=T)
        eval1 <- eval2 <- results <- NULL
        results <- gam(formula=STEP.formula, family=GAM.family, data=TrainData, weights=Yweights, maxit=50, bf.maxit=50)
        results <- step.gam(results, scope=GAMSTEP.scope, direction="both", trace=F, steps=GAMSTEP.steps, maxit=50, bf.maxit=50)
        cat(paste("\n", "stepwise GAM calibration","\n", "\n", sep = ""))
        TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="response")
        TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="response")
        eval1 <- evaluate(p=TrainPres, a=TrainAbs)
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "stepwise GAM evaluation","\n", "\n", "\n", sep = ""))
            TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="response")
            TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="response")
            eval2 <-  evaluate(p=TestPres, a=TestAbs)
            print(eval2)
            if(PLOTS==T) {
                plot(eval2, "ROC") 
                title(main="STEPWISE GAM (gam)", cex=2, adj=0, col.main="blue")
            }
        }
# conflict with mgcv gam function
        detach(package:gam)
        if(evaluations.keep ==T) {
            evaluations$GAMS.C <- eval1
            evaluations$GAMS.T <- eval2
        }
    }
    if (MGCV > 0) {
        require(mgcv, quietly=T)
        eval1 <- eval2 <- results <- NULL
        results <- gam(formula=MGCV.formula, family=GAM.family, data=TrainData, weights=Yweights, na.action=na.omit)
        cat(paste("\n", "GAM calibration (mgcv package)","\n", "\n", sep = ""))
        eval1 <- evaluate(p=p, a=a, model=results, x=x, type="response", na.action=na.omit)
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "GAM evaluation (mgcv package)","\n", "\n", sep = ""))
            eval2 <- evaluate(p=pt, a=at, model=results, x=x, type="response", na.action=na.omit)
            print(eval2)
            if(PLOTS==T) {
                plot(eval2, "ROC") 
                title(main="GAM (mgcv)", cex=2, adj=0, col.main="blue")
            }
        }
# conflict with gam package
        detach(package:mgcv)
        if(evaluations.keep ==T) {
            evaluations$MGCV.C <- eval1
            evaluations$MGCV.T <- eval2
        }
    }
    if (!is.null(factors) && EARTH > 0) {
        cat(paste("\n", "NOTE: MARS evaluation (earth package) with factors probably requires dummy variables", sep=""))
    } 
    if (EARTH > 0) {
#        require(earth, quietly=T)
        eval1 <- eval2 <- results <- NULL
        results <- earth(formula=EARTH.formula, glm=EARTH.glm, data=TrainData, degree=2)
        cat(paste("\n", "MARS calibration (earth package)","\n", "\n", sep = ""))
        eval1 <- evaluate(p=p, a=a, model=results, x=x, type="response")
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "MARS evaluation (earth package)","\n", "\n", sep = ""))
            eval2 <- evaluate(p=pt, a=at, model=results, x=x, type="response")
            print(eval2)
            if(PLOTS==T) {
                plot(eval2, "ROC") 
                title(main="MARS (earth)", cex=2, adj=0, col.main="blue")
            }
        }
#       detach(package:earth)
        if(evaluations.keep ==T) {
            evaluations$EARTH.C <- eval1
            evaluations$EARTH.T <- eval2
        }
    }
    if (RPART > 0) {
#        require(rpart, quietly=T)
        eval1 <- eval2 <- results <- NULL
        results <- rpart(formula=RPART.formula, data=TrainData, weights=Yweights,
            control=rpart.control(xval=RPART.xval, minbucket=5, minsplit=5, cp=0.001, maxdepth=25))
        cat(paste("\n", "RPART calibration","\n", "\n", sep = "")) 
        TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="prob")[,2]
        TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="prob")[,2]
        eval1 <- evaluate(p=TrainPres, a=TrainAbs)
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "RPART evaluation","\n", "\n", sep = ""))
            TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="prob")[,2]
            TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="prob")[,2]
            eval2 <-  evaluate(p=TestPres, a=TestAbs)
            print(eval2)
            if(PLOTS==T) {
                plot(eval2, "ROC") 
                title(main="RPART", cex=2, adj=0, col.main="blue")
            }
        }
#       detach(package:rpart)
        if(evaluations.keep ==T) {
            evaluations$RPART.C <- eval1
            evaluations$RPART.T <- eval2
        }
    }
    if (NNET > 0) {
#        require(nnet, quietly=T)
        eval1 <- eval2 <- results <- NULL
        results <- nnet(formula=NNET.formula, size=NNET.size, decay=NNET.decay, data=TrainData, weights=Yweights, 
            rang=0.1, maxit=200, trace=F)
        cat(paste("\n", "ANN calibration (nnet package)","\n", "\n", sep = ""))
        eval1 <- evaluate(p=p, a=a, model=results, x=x, type="raw")
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "ANN evaluation (nnet package)","\n", "\n", sep = ""))
            eval2 <- evaluate(p=pt, a=at, model=results, x=x, type="raw")
            print(eval2)
            if(PLOTS==T) {
                plot(eval2, "ROC") 
                title(main="NNET", cex=2, adj=0, col.main="blue")
            }
        }
#       detach(package:nnet)
        if(evaluations.keep ==T) {
            evaluations$NNET.C <- eval1
            evaluations$NNET.T <- eval2
        }
    }
    if (FDA > 0) {
#        require(mda, quietly=T)
        eval1 <- eval2 <- results <- NULL
        results <- fda(formula=FDA.formula, method=mars, data=TrainData, weights=Yweights)
        cat(paste("\n", "FDA calibration","\n", "\n", sep = "")) 
        TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="posterior")[,2]
        TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="posterior")[,2]
        eval1 <- evaluate(p=TrainPres, a=TrainAbs)
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "FDA evaluation","\n", "\n", sep = ""))
            TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="posterior")[,2]
            TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="posterior")[,2]
            eval2 <-  evaluate(p=TestPres, a=TestAbs)
            print(eval2)
            if(PLOTS==T) {
                plot(eval2, "ROC") 
                title(main="FDA", cex=2, adj=0, col.main="blue")
            }
        }
#       detach(package:mda)
        if(evaluations.keep ==T) {
            evaluations$FDA.C <- eval1
            evaluations$FDA.T <- eval2
        }
    }
    if (SVM > 0) {
#        require(kernlab, quietly=T)
        eval1 <- eval2 <- results <- NULL
        results <- ksvm(SVM.formula, data=TrainData, type="C-svc", prob.model=T)
        cat(paste("\n", "SVM calibration (kernlab package)","\n", "\n", sep = ""))
        TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="probabilities")[,2]
        TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="probabilities")[,2]
        eval1 <- evaluate(p=TrainPres, a=TrainAbs)
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            if (missinglevels==T) {
                cat(paste("\n", "SVM evaluations not done as factor levels are different","\n", "\n", sep = ""))
            }else{
                cat(paste("\n", "SVM evaluation (kernlab package)","\n", "\n", sep = ""))
                TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="probabilities")[,2]
                TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="probabilities")[,2]
                eval2 <-  evaluate(p=TestPres, a=TestAbs)
                print(eval2)
                if(PLOTS==T) {
                    plot(eval2, "ROC") 
                    title(main="SVM", cex=2, adj=0, col.main="blue")
                }
            }
        }
#       detach(package:kernlab)
        if(evaluations.keep ==T) {
            evaluations$SVM.C <- eval1
            evaluations$SVM.T <- eval2
        }
    }
    if (BIOCLIM > 0) {
        eval1 <- eval2 <- results <- NULL
        results <- bioclim(x=x, p=p)
        cat(paste("\n", "BIOCLIM calibration","\n", "\n", sep = ""))
        eval1 <- evaluate(p=p, a=a, model=results, x=x)
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "BIOCLIM evaluation","\n", "\n", sep = ""))
            eval2 <- evaluate(p=pt, a=at, model=results, x=x)
            print(eval2)
            if(PLOTS==T) {
                plot(eval2, "ROC") 
                title(main="BIOCLIM", cex=2, adj=0, col.main="blue")
            }
        }
        if(evaluations.keep ==T) {
            evaluations$BIOCLIM.C <- eval1
            evaluations$BIOCLIM.T <- eval2
        }
    }
    if (DOMAIN > 0) {
        eval1 <- eval2 <- results <- NULL
        results <- domain(x=x, p=p, factors=factors)
        cat(paste("\n", "DOMAIN calibration","\n", "\n", sep = ""))
        eval1 <- evaluate(p=p, a=a, model=results, x=x)
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "DOMAIN evaluation","\n", "\n", sep = ""))
            eval2 <- evaluate(p=pt, a=at, model=results, x=x)
            print(eval2)
            if(PLOTS==T) {
                plot(eval2, "ROC") 
                title(main="DOMAIN", cex=2, adj=0, col.main="blue")
            }
        }
        if(evaluations.keep ==T) {
            evaluations$DOMAIN.C <- eval1
            evaluations$DOMAIN.T <- eval2
        }
    }
    if (MAHAL > 0) {
        eval1 <- eval2 <- results <- NULL
        results <- mahal(x=x, p=p)
        cat(paste("\n", "Mahalanobis calibration","\n", "\n", sep = ""))
        eval1 <- evaluate(p=p, a=a, model=results, x=x)
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "Mahalanobis evaluation","\n", "\n", sep = ""))
            eval2 <- evaluate(p=pt, a=at, model=results, x=x)
            print(eval2)
            if(PLOTS==T) {
                plot(eval2, "ROC") 
                title(main="MAHAL", cex=2, adj=0, col.main="blue")
            }
        }
        if(evaluations.keep ==T) {
            evaluations$MAHAL.C <- eval1
            evaluations$MAHAL.T <- eval2
        }
    }
    cat(paste("\n","end of evaluations", "\n", "\n", sep = ""))
    remove(Yweights, envir=.BiodiversityR)
    remove(TrainData, envir=.BiodiversityR)
    remove(TestData, envir=.BiodiversityR)
    if (evaluations.keep==T) {return(evaluations)}
}




