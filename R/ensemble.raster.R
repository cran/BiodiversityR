`ensemble.raster` <- function(
    x, p, a=NULL, an=1000, ext=NULL, k=0, pt=NULL, at=NULL, xn=x, models.keep=FALSE, 
    RASTER.file.name="Species001", RASTER.format="raster", RASTER.datatype="INT2S", RASTER.NAflag=-32767,
    RASTER.models.overwrite=TRUE,
    ENSEMBLE.decay=1, ENSEMBLE.multiply=TRUE, ENSEMBLE.best=0,
    MAXENT=1, GBM=1, GBMSTEP=1, RF=1, GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=1, MGCV=1, 
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, BIOCLIM=1, DOMAIN=1, MAHAL=1,  
    Yweights="BIOMOD", factors=NULL,
    evaluation.strip=TRUE, 
    formulae.defaults=TRUE, 
    MAXENT.OLD=NULL,
    GBM.formula=NULL, GBM.n.trees=3000, GBM.OLD=NULL,    
    GBMSTEP.gbm.x=2:(1+nlayers(x)), GBMSTEP.tree.complexity=5, GBMSTEP.learning.rate=0.01, GBMSTEP.bag.fraction=0.5, GBMSTEP.OLD=NULL,
    RF.formula=NULL, RF.ntree=750, RF.mtry=log(nlayers(x)), RF.OLD=NULL,
    GLM.formula=NULL, GLM.family=binomial(link="logit"), GLM.OLD=NULL,
    GLMSTEP.steps=1000, STEP.formula=NULL, GLMSTEP.scope=NULL, GLMSTEP.k=2, GLMSTEP.OLD=NULL,
    GAM.formula=NULL, GAM.family=binomial(link="logit"), GAM.OLD=NULL, 
    GAMSTEP.steps=1000, GAMSTEP.scope=NULL, GAMSTEP.OLD=NULL,
    MGCV.formula=NULL, MGCV.OLD=NULL,
    EARTH.formula=NULL, EARTH.glm=list(family=binomial(link="logit")), EARTH.OLD=NULL,
    RPART.formula=NULL, RPART.xval=50, RPART.OLD=NULL,
    NNET.formula=NULL, NNET.size=8, NNET.decay=0.01, NNET.OLD=NULL,
    FDA.formula=NULL, FDA.OLD=NULL, 
    SVM.formula=NULL, SVM.OLD=NULL,
    BIOCLIM.OLD=NULL, DOMAIN.OLD=NULL, MAHAL.OLD=NULL    
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
    weights <- as.numeric(c(MAXENT, GBM, GBMSTEP, RF, GLM, GLMSTEP, GAM, GAMSTEP,
        MGCV, EARTH, RPART, NNET, FDA, SVM, BIOCLIM, DOMAIN, MAHAL))
    ws <- ensemble.weights(weights, decay=ENSEMBLE.decay, multiply=ENSEMBLE.multiply,
        best=ENSEMBLE.best, scale=TRUE)
    ws <- round(ws, digits=4)  
    cat(paste("\n", "Weights for ensemble forecasting", "\n", sep = ""))
    cat(paste("\n", "MAXENT=", ws[1], ", GBM=", ws[2], ", GBMSTEP=", ws[3], ", RF=", ws[4], 
            ", GLM=", ws[5], ", GLMSTEP=", ws[6],  
            ", GAM=", ws[7], ", GAMSTEP=", ws[8], ",", "\n",
            "MGCV=", ws[9], ", EARTH=", ws[10], ", RPART=", ws[11],
            ", NNET=", ws[12], ", FDA=", ws[13], ", SVM=", ws[14],", BIOCLIM=", ws[15],
            ", DOMAIN=", ws[16], ", MAHAL=", ws[17], "\n", "\n", sep=""))
#    if(is.null(factors)==F) {
#        for (i in 1:length(factors)) {
#            j <- which(names(x) == factors[i])
#            x[[j]] <- as.factor(x[[j]])
#            j <- which(names(xn) == factors[i])
#            xn[[j]] <- as.factor(xn[[j]])
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
    if(evaluation.strip==T) {
        ResponseData <- evaluation.strip.data(x, ext=ext, vars=names(x), factors=factors)
        assign("ResponseData", ResponseData, envir=.BiodiversityR)
    }
    if(models.keep==T) {
        ensemble.models <- list(evaluation.strip=ResponseData, MAXENT=NULL, GBM=NULL, GBMSTEP=NULL, RF=NULL, GLM=NULL, 
        GLMSTEP=NULL, GAM=NULL, GAMSTEP=NULL, MGCV=NULL, EARTH=NULL, RPART=NULL, 
        NNET=NULL, FDA=NULL, SVM=NULL, BIOCLIM=NULL, DOMAIN=NULL, MAHAL=NULL)
    }else{
        ensemble.models <- list(evaluation.strip=ResponseData)
    }
    dir.create("models", showWarnings = F)
    dir.create("ensembles", showWarnings = F)
    rasterfull <- paste("ensembles/", RASTER.file.name, "_ENSEMBLE", sep="")
    rastercount <- paste("ensembles/", RASTER.file.name, "_ENSEMBLE_COUNT", sep="")
    rasterpresence <- paste("ensembles/", RASTER.file.name, "_ENSEMBLE_PRESENCE", sep="")
    if (RASTER.models.overwrite==T) {RASTER.file.name <- "working"}
    ensemble <- raster(xn[[1]])
    ensemble[,] <- 0
    if (is.null(ext) == F) {ensemble <- crop(ensemble, y=ext, snap="in")}
#    writeRaster(x=ensemble, filename=rasterfull, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag, NAflag=RASTER.NAflag)
    enscount <- ensemble
    enspresence <- ensemble
#    writeRaster(x=enscount, filename=rastercount, progress='window', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
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
        # Put the file 'maxent.jar' in the 'java' folder of dismo
        # the file 'maxent.jar' can be obtained from from http://www.cs.princeton.edu/~schapire/maxent/.
        jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
        eval1 <- eval2 <- results <- NULL
        if(is.null(MAXENT.OLD) == T) {
            results <- maxent(x=x, p=p, a=a, factors=factors)
        }else{ results <- MAXENT.OLD}
        cat(paste("\n", "MAXENT calibration","\n", "\n", sep = ""))
        eval1 <- evaluate(p=p, a=a, model=results, x=x)
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "MAXENT evaluation","\n", "\n", "\n", sep = ""))
            eval2 <- evaluate(p=pt, a=at, model=results, x=x)
            print(eval2)
        }
        fullname <- paste("models/", RASTER.file.name, "_MAXENT", sep="")
        pmaxent <- predict(xn, results, ext=ext, 
             filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format)
        pmaxent <- trunc(1000*pmaxent)
        writeRaster(x=pmaxent, filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        if (models.keep==T) {ensemble.models$MAXENT <- results}
        ensemble <- ensemble + ws[1] * pmaxent
        writeRaster(x=ensemble, filename=rasterfull, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        threshold <- eval1@t[which.max(eval1@TPR + eval1@TNR)] * 1000        
        pmaxent <- pmaxent > threshold
        enscount <- enscount + pmaxent
        writeRaster(x=enscount, filename=rastercount, progress='window', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
        if (evaluation.strip == T) {
            ensemble.models$evaluation.strip$MAXENT <- as.numeric(predict(results, ResponseData))
            ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + ws[1] * ensemble.models$evaluation.strip$MAXENT            
        }    
    }
    if (GBM > 0) {
#        require(gbm, quietly=T) 
        eval1 <- eval2 <- results <- NULL
        if(is.null(GBM.OLD) == T) {       
        results <- gbm(formula=GBM.formula, data=TrainData, weights=Yweights, distribution = "bernoulli", 
            interaction.depth = 7, shrinkage = 0.001, bag.fraction = 0.5, train.fraction = 1, 
            n.trees = GBM.n.trees, verbose = F, cv.folds = 5)
        }else{ results <- GBM.OLD}
        cat(paste("\n", "GBM calibration","\n", "\n", sep = ""))
        eval1 <- evaluate(p=p, a=a, model=results, x=x, n.trees=results$n.trees, type="response")
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "GBM evaluation","\n", "\n", sep = ""))
            eval2 <- evaluate(p=pt, a=at, model=results, x=x, n.trees=results$n.trees, type="response")
            print(eval2)
        }
        fullname <- paste("models/", RASTER.file.name, "_GBM", sep="")
        pgbm <- predict(xn, results, ext=ext, n.trees=results$n.trees, type="response", 
             filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format)
        pgbm <- trunc(1000*pgbm)
        writeRaster(x=pgbm, filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        if (models.keep==T) {ensemble.models$GBM <- results}
        ensemble <- ensemble + ws[2] * pgbm
        writeRaster(x=ensemble, filename=rasterfull, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        threshold <- eval1@t[which.max(eval1@TPR + eval1@TNR)] * 1000        
        pgbm <- pgbm > threshold
        enscount <- enscount + pgbm
        writeRaster(x=enscount, filename=rastercount, progress='window', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
        if (evaluation.strip == T) {
            ensemble.models$evaluation.strip$GBM <- as.numeric(predict(results, newdata=ResponseData, n.trees=results$n.trees, type="response"))
            ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + ws[2] * ensemble.models$evaluation.strip$GBM
        }    
#        detach(package:gbm)
    }
    if (GBMSTEP > 0) {
#        require(gbm, quietly=T) 
        eval1 <- eval2 <- results <- NULL
        if(is.null(GBMSTEP.OLD) == T) {       
        results <- gbm.step(data=TrainData, gbm.y=1, gbm.x=GBMSTEP.gbm.x, family="bernoulli",
            site.weights=Yweights, tree.complexity=GBMSTEP.tree.complexity, learning.rate = GBMSTEP.learning.rate, 
            bag.fraction=GBMSTEP.bag.fraction, verbose=F, silent=T, plot.main=F)
        }else{ results <- GBMSTEP.OLD}
        cat(paste("\n", "stepwise GBM trees (target > 1000)", "\n", sep = ""))
        print(results$n.trees)
        cat(paste("\n", "stepwise GBM calibration","\n", "\n", sep = ""))
        eval1 <- evaluate(p=p, a=a, model=results, x=x, n.trees=results$n.trees, type="response")
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "stepwise GBM evaluation","\n", "\n", sep = ""))
            eval2 <- evaluate(p=pt, a=at, model=results, x=x, n.trees=results$n.trees, type="response")
            print(eval2)
        }
        fullname <- paste("models/", RASTER.file.name, "_GBMSTEP", sep="")
        pgbms <- predict(xn, results, ext=ext, n.trees=results$n.trees, type="response", 
             filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format)
        pgbms <- trunc(1000*pgbms)
        writeRaster(x=pgbms, filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        if (models.keep==T) {ensemble.models$GBMSTEP <- results}
        ensemble <- ensemble + ws[3] * pgbms
        writeRaster(x=ensemble, filename=rasterfull, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        threshold <- eval1@t[which.max(eval1@TPR + eval1@TNR)] * 1000        
        pgbms <- pgbms > threshold
        enscount <- enscount + pgbms
        writeRaster(x=enscount, filename=rastercount, progress='window', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
        if (evaluation.strip == T) {
            ensemble.models$evaluation.strip$GBMSTEP <- as.numeric(predict(results, newdata=ResponseData, n.trees=results$n.trees, type="response"))
            ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + ws[3] * ensemble.models$evaluation.strip$GBMSTEP
        }    
#        detach(package:gbm)
    }
    if (RF > 0) {
#        require(randomForest, quietly=T)
        eval1 <- eval2 <- results <- NULL
        if(is.null(RF.OLD) == T) {        
        results <- randomForest(formula=RF.formula, ntree=RF.ntree, mtry=RF.mtry, data=TrainData)
        }else{ results <- RF.OLD}
        cat(paste("\n", "RF calibration","\n", "\n", sep = ""))
#       treat factors as factors during evaluations
        TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="response")
        TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="response")
        eval1 <- evaluate(p=TrainPres, a=TrainAbs)
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "RF evaluation","\n", "\n", sep = ""))
            TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="response")
            TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="response")
            eval2 <-  evaluate(p=TestPres, a=TestAbs)
            print(eval2)
        }
        fullname <- paste("models/", RASTER.file.name, "_RF", sep="")
        prf <- predict(xn, results, ext=ext,
             filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format)
        prf <- trunc(1000*prf)
        writeRaster(x=prf, filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        if (models.keep==T) {ensemble.models$RF <- results}
        ensemble <- ensemble + ws[4] * prf
        writeRaster(x=ensemble, filename=rasterfull, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        threshold <- eval1@t[which.max(eval1@TPR + eval1@TNR)] * 1000        
        prf <- prf > threshold
        enscount <- enscount + prf
        writeRaster(x=enscount, filename=rastercount, progress='window', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
        if (evaluation.strip == T) {
            ensemble.models$evaluation.strip$RF <- as.numeric(predict(results, newdata=ResponseData, type="response"))
            ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + ws[4] * ensemble.models$evaluation.strip$RF
        } 
#        detach(package:randomForest)
    } 
    if (GLM > 0) {
        eval1 <- eval2 <- results <- NULL
        if(is.null(GLM.OLD) == T) { 
        results <- glm(formula=GLM.formula, family=GLM.family, data=TrainData, weights=Yweights, maxit=100)
        }else{ results <- GLM.OLD}
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
        }
        fullname <- paste("models/", RASTER.file.name, "_GLM", sep="")
        pglm <- predict(xn, results, ext=ext, type="response", 
             filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format)
        pglm <- trunc(1000*pglm)
        writeRaster(x=pglm, filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        if (models.keep==T) {ensemble.models$GLM <- results}
        ensemble <- ensemble + ws[5] * pglm
        writeRaster(x=ensemble, filename=rasterfull, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        threshold <- eval1@t[which.max(eval1@TPR + eval1@TNR)] * 1000        
        pglm <- pglm > threshold
        enscount <- enscount + pglm
        writeRaster(x=enscount, filename=rastercount, progress='window', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
        if (evaluation.strip == T) {
            ensemble.models$evaluation.strip$GLM <- as.numeric(predict(results, newdata=ResponseData, type="response"))
            ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + ws[5] * ensemble.models$evaluation.strip$GLM
        }    
    }
    if (GLMSTEP > 0) {
#        require(MASS, quietly=T)
        eval1 <- eval2 <- results <- NULL
        if(is.null(GLMSTEP.OLD) == T) {
            results <- glm(formula=STEP.formula, family=GLM.family, data=TrainData, weights=Yweights, maxit=100)
            results <- stepAIC(results, scope=GLMSTEP.scope, direction="both", trace=F, steps=GLMSTEP.steps, k=GLMSTEP.k, maxit=100)
        }else{ results <- GLMSTEP.OLD}
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
        }
        fullname <- paste("models/", RASTER.file.name, "_GLMSTEP", sep="")
        pglms <- predict(xn, results, ext=ext, type="response", 
             filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format)
        pglms <- trunc(1000*pglms)
        writeRaster(x=pglms, filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        if (models.keep==T) {ensemble.models$GLMSTEP <- results}
        ensemble <- ensemble + ws[6] * pglms
        writeRaster(x=ensemble, filename=rasterfull, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        threshold <- eval1@t[which.max(eval1@TPR + eval1@TNR)] * 1000        
        pglms <- pglms > threshold
        enscount <- enscount + pglms
        writeRaster(x=enscount, filename=rastercount, progress='window', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
        if (evaluation.strip == T) {
            ensemble.models$evaluation.strip$GLMSTEP <- as.numeric(predict(results, newdata=ResponseData, type="response"))
            ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + ws[6] * ensemble.models$evaluation.strip$GLMSTEP
        }  
#        detach(package:MASS)
    }
    if (GAM > 0) {
        require(gam, quietly=T)
        eval1 <- eval2 <- results <- NULL
        if(is.null(GAM.OLD) == T) {
        results <- gam(formula=GAM.formula, family=GAM.family, data=TrainData, weights=Yweights, maxit=50, bf.maxit=50)
        }else{ results <- GAM.OLD}
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
        }
        fullname <- paste("models/", RASTER.file.name, "_GAM", sep="")
        pgam <- predict(xn, results, ext=ext, type="response", 
             filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format)
        pgam <- trunc(1000*pgam)
        writeRaster(x=pgam, filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        if (models.keep==T) {ensemble.models$GAM <- results}
        ensemble <- ensemble + ws[7] * pgam
        writeRaster(x=ensemble, filename=rasterfull, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        threshold <- eval1@t[which.max(eval1@TPR + eval1@TNR)] * 1000        
        pgam <- pgam > threshold
        enscount <- enscount + pgam
        writeRaster(x=enscount, filename=rastercount, progress='window', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
        if (evaluation.strip == T) {
            ensemble.models$evaluation.strip$GAM <- as.numeric(predict(results, newdata=ResponseData, type="response"))
            ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + ws[7] * ensemble.models$evaluation.strip$GAM
        }
        detach(package:gam)
    }
    if (GAMSTEP > 0) {
        cat(paste("\n", "NOTE: stepwise GAM seems to require assignments of family and data to pos 1", "\n", "\n", sep=""))
        GAMSTEP <- 0
    } 
    if (GAMSTEP > 0) {
        require(gam, quietly=T)
        eval1 <- eval2 <- results <- NULL
        if(is.null(GAMSTEP.OLD) == T) {
            results <- gam(formula=STEP.formula, family=GAM.family, data=TrainData, weights=Yweights, maxit=50, bf.maxit=50)
            results <- step.gam(results, scope=GAMSTEP.scope, direction="both", trace=F, steps=GAMSTEP.steps, maxit=50, bf.maxit=50)
        }else{ results <- GAMSTEP.OLD}
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
        }
        fullname <- paste("models/", RASTER.file.name, "_GAMSTEP", sep="")
        pgams <- predict(xn, results, ext=ext, type="response", 
             filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format)
        pgams <- trunc(1000*pgams)
        writeRaster(x=pgams, filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        if (models.keep==T) {ensemble.models$GAMSTEP <- results}
        ensemble <- ensemble + ws[8] * pgams
        writeRaster(x=ensemble, filename=rasterfull, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        threshold <- eval1@t[which.max(eval1@TPR + eval1@TNR)] * 1000        
        pgams <- pgams > threshold
        enscount <- enscount + pgams
        writeRaster(x=enscount, filename=rastercount, progress='window', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
        if (evaluation.strip == T) {
            ensemble.models$evaluation.strip$GAMSTEP <- as.numeric(predict(results, newdata=ResponseData, type="response"))
            ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + ws[8] * ensemble.models$evaluation.strip$GAMSTEP
        }
        detach(package:gam)
    }
    if (MGCV > 0) {
        eval1 <- eval2 <- results <- NULL
        require(mgcv, quietly=T)
        if(is.null(MGCV.OLD) == T) {
        results <- gam(formula=MGCV.formula, family=GAM.family, data=TrainData, weights=Yweights)
        }else{ results <- MGCV.OLD}
        cat(paste("\n", "GAM calibration (mgcv package)","\n", "\n", sep = ""))
        eval1 <- evaluate(p=p, a=a, model=results, x=x, type="response")
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "GAM evaluation (mgcv package)","\n", "\n", sep = ""))
            eval2 <- evaluate(p=pt, a=at, model=results, x=x, type="response")
            print(eval2)
        }
        fullname <- paste("models/", RASTER.file.name, "_MGCV", sep="")
        pmgcv <- predict(xn, results, ext=ext, type="response", 
             filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format)
        pmgcv <- trunc(1000*pmgcv)
        writeRaster(x=pmgcv, filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        if (models.keep==T) {ensemble.models$MGCV <- results}
        ensemble <- ensemble + ws[9] * pmgcv
        writeRaster(x=ensemble, filename=rasterfull, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        threshold <- eval1@t[which.max(eval1@TPR + eval1@TNR)] * 1000        
        pmgcv <- pmgcv > threshold
        enscount <- enscount + pmgcv
        writeRaster(x=enscount, filename=rastercount, progress='window', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
        if (evaluation.strip == T) {
            ensemble.models$evaluation.strip$MGCV <- as.numeric(predict(results, newdata=ResponseData, type="response"))
            ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + ws[9] * ensemble.models$evaluation.strip$MGCV
        }
        detach(package:mgcv)
    }
    if (!is.null(factors) && EARTH > 0) {
        cat(paste("\n", "MARS model (earth package) with factors may require explicit dummy variables", sep=""))
    } 
    if (EARTH > 0) {
#        require(earth, quietly=T)
        eval1 <- eval2 <- results <- NULL
        if(is.null(EARTH.OLD) == T) {
        results <- earth(formula=EARTH.formula, glm=EARTH.glm, data=TrainData, weights=Yweights, degree=2)
        }else{ results <- EARTH.OLD}
        cat(paste("\n", "MARS calibration (earth package)","\n", "\n", sep = ""))
        eval1 <- evaluate(p=p, a=a, model=results, x=x, type="response")
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "MARS evaluation (earth package)","\n", "\n", sep = ""))
            eval2 <- evaluate(p=pt, a=at, model=results, x=x, type="response")
            print(eval2)
        }
        fullname <- paste("models/", RASTER.file.name, "_EARTH", sep="")
        pearth <- predict(xn, results, ext=ext, type="response", 
             filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format)
        pearth <- trunc(1000*pearth)
        writeRaster(x=pearth, filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        if (models.keep==T) {ensemble.models$EARTH <- results}
        ensemble <- ensemble + ws[10] * pearth
        writeRaster(x=ensemble, filename=rasterfull, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        threshold <- eval1@t[which.max(eval1@TPR + eval1@TNR)] * 1000        
        pearth <- pearth > threshold
        enscount <- enscount + pearth
        writeRaster(x=enscount, filename=rastercount, progress='window', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
        if (evaluation.strip == T) {
            ensemble.models$evaluation.strip$EARTH <- as.numeric(predict(results, newdata=ResponseData, type="response"))
            ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + ws[10] * ensemble.models$evaluation.strip$EARTH
        }
#        detach(package:earth)
    }
    if (RPART > 0) {
#        require(rpart, quietly=T)
        eval1 <- eval2 <- results <- NULL
        if(is.null(RPART.OLD) == T) {
        results <- rpart(formula=RPART.formula, data=TrainData, weights=Yweights,
            control=rpart.control(xval=RPART.xval, minbucket=5, minsplit=5, cp=0.001, maxdepth=25))
        }else{ results <- RPART.OLD}
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
        }
        fullname <- paste("models/", RASTER.file.name, "_RPART", sep="")
        prpart <- predict(xn, results, ext=ext, type="prob", index=2,
             filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format)
        prpart <- trunc(1000*prpart)
        writeRaster(x=prpart, filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        if (models.keep==T) {ensemble.models$RPART <- results}
        ensemble <- ensemble + ws[11] * prpart
        writeRaster(x=ensemble, filename=rasterfull, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        threshold <- eval1@t[which.max(eval1@TPR + eval1@TNR)] * 1000        
        prpart <- prpart > threshold
        enscount <- enscount + prpart
        writeRaster(x=enscount, filename=rastercount, progress='window', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
        if (evaluation.strip == T) {
            ensemble.models$evaluation.strip$RPART <- as.numeric(predict(results, newdata=ResponseData, type="prob")[,2])
            ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + ws[11] * ensemble.models$evaluation.strip$RPART
        }
#        detach(package:rpart)
    }
    if (!is.null(factors) && NNET > 0) {
        cat(paste("\n", "ANN model with factors may require explicit dummy variables", sep=""))
    } 
    if (NNET > 0) {
#        require(nnet, quietly=T)
        eval1 <- eval2 <- results <- NULL
        if(is.null(NNET.OLD) == T) {
        results <- nnet(formula=NNET.formula, size=NNET.size, decay=NNET.decay, data=TrainData, weights=Yweights, 
            rang=0.1, maxit=200, trace=F)
        }else{ results <- NNET.OLD}
        cat(paste("\n", "ANN calibration (nnet package)","\n", "\n", sep = ""))
        eval1 <- evaluate(p=p, a=a, model=results, x=x, type="raw")
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "ANN evaluation (nnet package)","\n", "\n", sep = ""))
            eval2 <- evaluate(p=pt, a=at, model=results, x=x, type="raw")
            print(eval2)
        }
        fullname <- paste("models/", RASTER.file.name, "_NNET", sep="")
        pnnet <- predict(xn, results, ext=ext, type="raw", 
             filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format)
        pnnet <- trunc(1000*pnnet)
        writeRaster(x=pnnet, filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        if (models.keep==T) {ensemble.models$NNET <- results}
        ensemble <- ensemble + ws[12] * pnnet
        writeRaster(x=ensemble, filename=rasterfull, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        threshold <- eval1@t[which.max(eval1@TPR + eval1@TNR)] * 1000        
        pnnet <- pnnet > threshold
        enscount <- enscount + pnnet
        writeRaster(x=enscount, filename=rastercount, progress='window', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
        if (evaluation.strip == T) {
            ensemble.models$evaluation.strip$NNET <- as.numeric(predict(results, newdata=ResponseData, type="raw"))
            ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + ws[12] * ensemble.models$evaluation.strip$NNET
        }
#        detach(package:nnet)
    }
    if (FDA > 0) {
#        require(mda, quietly=T)
        eval1 <- eval2 <- results <- NULL
        if(is.null(FDA.OLD) == T) {
        results <- fda(formula=FDA.formula, method=mars, data=TrainData, weights=Yweights)
        }else{ results <- FDA.OLD}
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
        }
        fullname <- paste("models/", RASTER.file.name, "_FDA", sep="")
        pfda <- predict(xn, results, ext=ext, type="posterior", index=2,
             filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format)
        pfda <- trunc(1000*pfda)
        writeRaster(x=pfda, filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        if (models.keep==T) {ensemble.models$FDA <- results}
        ensemble <- ensemble + ws[13] * pfda
        writeRaster(x=ensemble, filename=rasterfull, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        threshold <- eval1@t[which.max(eval1@TPR + eval1@TNR)] * 1000        
        pfda <- pfda > threshold
        enscount <- enscount + pfda
        writeRaster(x=enscount, filename=rastercount, progress='window', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
        if (evaluation.strip == T) {
            ensemble.models$evaluation.strip$FDA <- as.numeric(predict(results, newdata=ResponseData, type="posterior")[,2])
            ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + ws[13] * ensemble.models$evaluation.strip$FDA
        }
#        detach(package:mda)
    }
    if (!is.null(factors) && SVM > 0) {
        cat(paste("\n", "SVM model with factors may require explicit dummy variables", sep=""))
    } 
    if (SVM > 0) {
#        require(kernlab, quietly=T)
        eval1 <- eval2 <- results <- NULL
        if(is.null(SVM.OLD) == T) {
        results <- ksvm(SVM.formula, data=TrainData, type="C-svc", prob.model=T)
        }else{ results <- SVM.OLD}
        cat(paste("\n", "SVM calibration (kernlab package)","\n", "\n", sep = "")) 
        TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="probabilities")[,2]
        TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="probabilities")[,2]
        eval1 <- evaluate(p=TrainPres, a=TrainAbs)
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "SVM evaluation (kernlab package)","\n", "\n", sep = ""))
            TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="probabilities")[,2]
            TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="probabilities")[,2]
            eval2 <-  evaluate(p=TestPres, a=TestAbs)
            print(eval2)
        }
        fullname <- paste("models/", RASTER.file.name, "_SVM", sep="")
        psvm <- predict(xn, results, ext=ext, type="probabilities", index=2,
             filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format)
        psvm <- trunc(1000*psvm)
        writeRaster(x=psvm, filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        if (models.keep==T) {ensemble.models$SVM <- results}
        ensemble <- ensemble + ws[14] * psvm
        writeRaster(x=ensemble, filename=rasterfull, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        threshold <- eval1@t[which.max(eval1@TPR + eval1@TNR)] * 1000        
        psvm <- psvm > threshold
        enscount <- enscount + psvm
        writeRaster(x=enscount, filename=rastercount, progress='window', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
        if (evaluation.strip == T) {
            ensemble.models$evaluation.strip$SVM <- as.numeric(predict(results, newdata=ResponseData, type="probabilities")[,2])
            ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + ws[14] * ensemble.models$evaluation.strip$SVM
        }
#        detach(package:kernlab)
    }
    if (BIOCLIM > 0) {  
        eval1 <- eval2 <- results <- NULL  
        xnum <- x
        xnnum <- xn
        if(is.null(factors)==F) {
            for (i in 1:length(factors)) {
                xnum <- dropLayer(xnum, which(names(xnum) == factors[i]))
                xnnum <- dropLayer(xnnum, which(names(xnnum) == factors[i]))                
            }
        }
        if(is.null(BIOCLIM.OLD) == T) {
        results <- bioclim(x=xnum, p=p)
        }else{ results <- BIOCLIM.OLD}
        cat(paste("\n", "BIOCLIM calibration","\n", "\n", sep = ""))
        eval1 <- evaluate(p=p, a=a, model=results, x=xnum)
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "BIOCLIM evaluation","\n", "\n", sep = ""))
            eval2 <- evaluate(p=pt, a=at, model=results, x=xnum)
            print(eval2)
        }
        fullname <- paste("models/", RASTER.file.name, "_BIOCLIM", sep="")
        pbio <- predict(xnnum, results, ext=ext, 
             filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format)
        pbio <- trunc(1000*pbio)
        writeRaster(x=pbio, filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        if (models.keep==T) {ensemble.models$BIOCLIM <- results}
        ensemble <- ensemble + ws[15] * pbio
        writeRaster(x=ensemble, filename=rasterfull, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        threshold <- eval1@t[which.max(eval1@TPR + eval1@TNR)] * 1000        
        pbio <- pbio > threshold
        enscount <- enscount + pbio
        writeRaster(x=enscount, filename=rastercount, progress='window', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
        if (evaluation.strip == T) {
            ensemble.models$evaluation.strip$BIOCLIM <- as.numeric(predict(results, ResponseData))
            ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + ws[15] * ensemble.models$evaluation.strip$BIOCLIM
        }
    }
    if (DOMAIN > 0) {
        eval1 <- eval2 <- results <- NULL
        xnum <- x
        xnnum <- xn
        if(is.null(factors)==F) {
            for (i in 1:length(factors)) {
                xnum <- dropLayer(xnum, which(names(xnum) == factors[i]))
                xnnum <- dropLayer(xnnum, which(names(xnnum) == factors[i]))                
            }
        }
        if(is.null(DOMAIN.OLD) == T) {
        results <- domain(x=xnum, p=p, factors=factors)
        }else{ results <- DOMAIN.OLD}
        cat(paste("\n", "DOMAIN calibration","\n", "\n", sep = ""))
        eval1 <- evaluate(p=p, a=a, model=results, x=xnum)
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "DOMAIN evaluation","\n", "\n", sep = ""))
            eval2 <- evaluate(p=pt, a=at, model=results, x=xnum)
            print(eval2)
        }
        fullname <- paste("models/", RASTER.file.name, "_DOMAIN", sep="")
        pdom <- predict(xnnum, results, ext=ext, 
             filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format)
        pdom <- trunc(1000*pdom)
        writeRaster(x=pdom, filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        if (models.keep==T) {ensemble.models$DOMAIN <- results}
        ensemble <- ensemble + ws[16] * pdom
        writeRaster(x=ensemble, filename=rasterfull, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        threshold <- eval1@t[which.max(eval1@TPR + eval1@TNR)] * 1000        
        pdom <- pdom > threshold
        enscount <- enscount + pdom
        writeRaster(x=enscount, filename=rastercount, progress='window', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
        if (evaluation.strip == T) {
            ensemble.models$evaluation.strip$DOMAIN <- as.numeric(predict(results, ResponseData))
            ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + ws[16] * ensemble.models$evaluation.strip$DOMAIN
        }
    }
    if (MAHAL > 0) {  
        eval1 <- eval2 <- results <- NULL  
        xnum <- x
        xnnum <- xn
        if(is.null(factors)==F) {
            for (i in 1:length(factors)) {
                xnum <- dropLayer(xnum, which(names(xnum) == factors[i]))
                xnnum <- dropLayer(xnnum, which(names(xnnum) == factors[i]))                
            }
        }
        if(is.null(MAHAL.OLD) == T) {
        results <- mahal(x=xnum, p=p)
        }else{ results <- MAHAL.OLD}
        cat(paste("\n", "Mahalanobis calibration","\n", "\n", sep = ""))
        eval1 <- evaluate(p=p, a=a, model=results, x=xnum)
        print(eval1)
        if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
            cat(paste("\n", "Mahalanobis evaluation","\n", "\n", sep = ""))
            eval2 <- evaluate(p=pt, a=at, model=results, x=xnum)
            print(eval2)
        }
        fullname <- paste("models/", RASTER.file.name, "_MAHAL", sep="")
# first calculate presence
        pmahal <- predict(xnnum, results, ext=ext, 
             filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format)
        threshold <- eval1@t[which.max(eval1@TPR + eval1@TNR)]        
        pmahal <- pmahal > threshold
        enscount <- enscount + pmahal
        writeRaster(x=enscount, filename=rastercount, progress='window', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
        pmahal <- predict(xnnum, results, ext=ext, 
             filename=fullname, progress='window', overwrite=TRUE)
# avoid large values
        setMinMax(pmahal)
        mahalmin <- minValue(pmahal)
        mahalrange <- (1 - mahalmin)
        pmahal <- pmahal - mahalmin
        pmahal <- pmahal / mahalrange
        pmahal <- trunc(1000*pmahal)
        writeRaster(x=pmahal, filename=fullname, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        if (models.keep==T) {ensemble.models$MAHAL <- results}
        ensemble <- ensemble + ws[17] * pmahal
        writeRaster(x=ensemble, filename=rasterfull, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        if (evaluation.strip == T) {
            ensemble.models$evaluation.strip$MAHAL <- as.numeric(predict(results, ResponseData))
            mahalmin <- min(ensemble.models$evaluation.strip$MAHAL)
            mahalrange <- (1 - mahalmin)
            ensemble.models$evaluation.strip$MAHAL <- ensemble.models$evaluation.strip$MAHAL - mahalmin
            ensemble.models$evaluation.strip$MAHAL <- ensemble.models$evaluation.strip$MAHAL / mahalrange
            ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + ws[17] * ensemble.models$evaluation.strip$MAHAL
        }
    }
    ensemble <- trunc(ensemble)
    writeRaster(x=ensemble, filename=rasterfull, progress='window', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
    cat(paste("\n", "Ensemble calibration","\n", "\n", sep = ""))
    pres_consensus <- extract(ensemble, p)
    abs_consensus <- extract(ensemble, a)
    eval1 <- evaluate(p=pres_consensus, a=abs_consensus)
    print(eval1)
    l1 <- minValue(ensemble)
    l5 <- maxValue(ensemble)
    l3 <- eval1@t[which.max(eval1@TPR + eval1@TNR)]
    r1 <- l3 - l1
    r2 <- l5 - l3 
    l2 <- l1 + r1/2
    l4 <- l3 + r2/2
    cat(paste("\n", "Suggested thresholds for absence", sep = ""))
    cat(paste("\n", "(includes minimum, midpoint and threshold)","\n", "\n", sep = ""))     
    cat(paste(l1, ", ", l2, ", ", l3, "\n", sep=""))
    cat(paste("\n", "Suggested thresholds for presence", sep = ""))
    cat(paste("\n", "(includes threshold, midpoint and maximum)","\n", "\n", sep = ""))      
    cat(paste(l3, ", ", l4, ", ", l5, "\n", "\n", sep=""))
    enspresence <- ensemble > l3
    writeRaster(x=enspresence, filename=rasterpresence, progress='window', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
        cat(paste("\n", "Ensemble evaluation","\n", "\n", sep = ""))
        pres_consensus <- extract(ensemble, pt)
        abs_consensus <- extract(ensemble, at)
        eval2 <- evaluate(p=pres_consensus, a=abs_consensus)
        print(eval2)
        l1 <- minValue(ensemble)
        l5 <- maxValue(ensemble)
        l3 <- eval2@t[which.max(eval2@TPR + eval2@TNR)]
        r1 <- l3 - l1
        r2 <- l5 - l3 
        l2 <- l1 + r1/2
        l4 <- l3 + r2/2
        cat(paste("\n", "Suggested thresholds for absence (evaluation)", sep = ""))
        cat(paste("\n", "(includes minimum, midpoint and threshold)","\n", "\n", sep = ""))       
        cat(paste(l1, ", ", l2, ", ", l3, "\n", sep=""))
        cat(paste("\n", "Suggested thresholds for presence (evaluation)", sep = ""))
        cat(paste("\n", "(includes threshold, midpoint and maximum)","\n", "\n", sep = ""))     
        cat(paste(l3, ", ", l4, ", ", l5, "\n", "\n", sep=""))
    }
    remove(Yweights, envir=.BiodiversityR)
    remove(TrainData, envir=.BiodiversityR)
    remove(TestData,envir=.BiodiversityR)
    if (evaluation.strip==T) {remove(ResponseData, envir=.BiodiversityR)}
    if (models.keep==T || evaluation.strip==T) {return(ensemble.models)}
}




