`ensemble.raster` <- function(
    x, p, a=NULL, an=1000, ext=NULL, k=0, pt=NULL, at=NULL, xn=x, 
    layer.drops=NULL,
    models.keep=FALSE, 
    RASTER.species.name="Species001", RASTER.stack.name=xn@title, RASTER.format="raster", RASTER.datatype="INT2S", RASTER.NAflag=-32767,
    RASTER.models.overwrite=TRUE,
    threshold.method="spec_sens",
    ENSEMBLE.decay=1, ENSEMBLE.multiply=TRUE, ENSEMBLE.best=0,
    input.weights=NULL,
    MAXENT=1, GBM=1, GBMSTEP=1, RF=1, GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=1, MGCV=1, MGCVFIX=0,
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, BIOCLIM=1, DOMAIN=1, MAHAL=1,  
    Yweights="BIOMOD", factors=NULL, dummy.vars=NULL,
    evaluation.strip=TRUE, 
    formulae.defaults=TRUE, maxit=100,
    MAXENT.OLD=NULL,
    GBM.formula=NULL, GBM.n.trees=2000, GBM.OLD=NULL,    
    GBMSTEP.gbm.x=2:(1+nlayers(x)), GBMSTEP.tree.complexity=5, GBMSTEP.learning.rate=0.005, 
    GBMSTEP.bag.fraction=0.5, GBMSTEP.step.size=100, GBMSTEP.OLD=NULL,
    RF.formula=NULL, RF.ntree=750, RF.mtry=log(nlayers(x)), RF.OLD=NULL,
    GLM.formula=NULL, GLM.family=binomial(link="logit"), GLM.OLD=NULL,
    GLMSTEP.steps=1000, STEP.formula=NULL, GLMSTEP.scope=NULL, GLMSTEP.k=2, GLMSTEP.OLD=NULL,
    GAM.formula=NULL, GAM.family=binomial(link="logit"), GAM.OLD=NULL, 
    GAMSTEP.steps=1000, GAMSTEP.scope=NULL, GAMSTEP.OLD=NULL,
    MGCV.formula=NULL, MGCV.select=FALSE, MGCV.OLD=NULL,
    MGCVFIX.formula=NULL, MGCVFIX.OLD=NULL,
    EARTH.formula=NULL, EARTH.glm=list(family=binomial(link="logit"), maxit=maxit), EARTH.OLD=NULL,
    RPART.formula=NULL, RPART.xval=50, RPART.OLD=NULL,
    NNET.formula=NULL, NNET.size=8, NNET.decay=0.01, NNET.OLD=NULL,
    FDA.formula=NULL, FDA.OLD=NULL, 
    SVM.formula=NULL, SVM.OLD=NULL, 
    SVME.formula=NULL, SVME.OLD=NULL,
    BIOCLIM.OLD=NULL, DOMAIN.OLD=NULL, MAHAL.OLD=NULL    
)
{
    .BiodiversityR <- new.env()
    if (! require(dismo)) {stop("Please install the dismo package")}
    if (is.null(xn) == T) {
        cat(paste("\n", "NOTE: new rasterStack assumed to be equal to the base rasterStack", sep = ""))
        xn <- x
    }
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
                cat(paste("\n", "WARNING: categorical variable '", factors[i], "' not among grid layers", "\n\n", sep = ""))
            }
        }
    }
    if (is.null(dummy.vars) == F) {
        vars <- names(x)
        dummy.vars <- as.character(dummy.vars)
        nf <- length(dummy.vars)
        for (i in 1:nf) {
            if (any(vars==dummy.vars[i])==FALSE) {
                cat(paste("\n", "WARNING: dummy variable '", dummy.vars[i], "' not among grid layers", "\n\n", sep = ""))
            }
        }
    }
    if (is.null(input.weights)==F) {
# use the last column in case output from the ensemble.test.splits function is used
        if (length(dim(input.weights)) == 2) {input.weights <- input.weights[,"MEAN"]}
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
    weights <- as.numeric(c(MAXENT, GBM, GBMSTEP, RF, GLM, GLMSTEP, GAM, GAMSTEP, MGCV, MGCVFIX, 
        EARTH, RPART, NNET, FDA, SVM, SVME, BIOCLIM, DOMAIN, MAHAL))
    ws <- ensemble.weights(weights, decay=ENSEMBLE.decay, multiply=ENSEMBLE.multiply,
        best=ENSEMBLE.best, scale=TRUE)
    ws <- round(ws, digits=4)
    names(ws) <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX", 
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL")
    cat(paste("\n", "Weights for ensemble forecasting", "\n", sep = ""))
    print(ws)
    prediction.failures <- FALSE
    if(is.null(factors)==F) {
        for (i in 1:length(factors)) {
            j <- which(names(x) == factors[i])
            x[[j]] <- raster::as.factor(x[[j]])
#            levels1 <- levels(as.factor(x[[j]]))[[1]]$ID
            j <- which(names(xn) == factors[i])
            xn[[j]] <- raster::as.factor(xn[[j]])
#            levels2 <- levels(as.factor(xn[[j]]))[[1]]$ID
        }
    }
    TrainData <- prepareData(x, p, b=a, factors=factors, xy=FALSE)
    if(any(is.na(TrainData[TrainData[,"pb"]==1,]))) {
        cat(paste("\n", "WARNING: presence locations with missing data removed from calibration data","\n\n", sep = ""))
    }
    TrainValid <- complete.cases(TrainData[TrainData[,"pb"]==1,])
    p <- p[TrainValid,]
    if(any(is.na(TrainData[TrainData[,"pb"]==0,]))) {
        cat(paste("\n", "WARNING: background locations with missing data removed from calibration data","\n\n", sep = ""))
    }
    TrainValid <- complete.cases(TrainData[TrainData[,"pb"]==0,])
    a <- a[TrainValid,]
    TrainData <- prepareData(x, p, b=a, factors=factors, xy=FALSE)   
    TestData <- prepareData(x, p=pt, b=at, factors=factors, xy=FALSE)
    if(any(is.na(TestData[TestData[,"pb"]==1,]))) {
        cat(paste("\n", "WARNING: presence locations with missing data removed from evaluation data","\n\n", sep = ""))
    }
    TestValid <- complete.cases(TestData[TestData[,"pb"]==1,])
    pt <- pt[TestValid,]
    if(any(is.na(TestData[TestData[,"pb"]==0,]))) {
        cat(paste("\n", "WARNING: background locations with missing data removed from evaluation data","\n\n", sep = ""))
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
    if(evaluation.strip==T) {
        ResponseData <- evaluation.strip.data(x, ext=ext, vars=names(x), factors=factors)
    }
    if(models.keep==T  && evaluation.strip==T) {
        ensemble.models <- list(evaluation.strip=ResponseData, p=p, a=a, pt=pt, at=at, 
            MAXENT=NULL, GBM=NULL, GBMSTEP=NULL, RF=NULL, GLM=NULL, GLMSTEP=NULL, GAM=NULL, GAMSTEP=NULL, MGCV=NULL, MGCVFIX=NULL, 
            EARTH=NULL, RPART=NULL, NNET=NULL, FDA=NULL, SVM=NULL, SVME=NULL, BIOCLIM=NULL, DOMAIN=NULL, MAHAL=NULL)
    }
    if(models.keep==T  && evaluation.strip==F) {
        ensemble.models <- list(p=p, a=a, pt=pt, at=at, MAXENT=NULL, GBM=NULL, GBMSTEP=NULL, RF=NULL, GLM=NULL, 
        GLMSTEP=NULL, GAM=NULL, GAMSTEP=NULL, MGCV=NULL, EARTH=NULL, RPART=NULL, 
        NNET=NULL, FDA=NULL, SVM=NULL, SVME=NULL, BIOCLIM=NULL, DOMAIN=NULL, MAHAL=NULL)
    }
    if(models.keep==F  && evaluation.strip==T) {
        ensemble.models <- list(evaluation.strip=ResponseData)
    }
    if(evaluation.strip==T) {
# create data set only with explanatory variables - avoid problems with SVME
        vars <- names(x)
        ResponseData <- ResponseData[,vars]
        assign("ResponseData", ResponseData, envir=.BiodiversityR)
    }
    dir.create("models", showWarnings = F)
    dir.create("ensembles", showWarnings = F)
    if(length(x@title) == 0) {x@title <- "base"}
    if(length(xn@title) == 0) {xn@title <- "newstack"}
    rasterfull <- paste("ensembles/", RASTER.species.name, "_ENSEMBLE_", RASTER.stack.name , sep="")
    rastercount <- paste("ensembles/", RASTER.species.name, "_COUNT_", RASTER.stack.name , sep="")
    rasterpresence <- paste("ensembles/", RASTER.species.name, "_PRESENCE_", RASTER.stack.name, sep="")
    RASTER.species.orig <- RASTER.species.name
    if (RASTER.models.overwrite==T) {
        RASTER.species.name <- "working"
    }else{
        RASTER.species.name <- paste(RASTER.species.name, "_", RASTER.stack.name, sep="")
    }
    ensemble <- raster(xn[[1]])
    ensemble[,] <- 0
    if (is.null(ext) == F) {ensemble <- crop(ensemble, y=ext, snap="in")}
#    writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag, NAflag=RASTER.NAflag)
    enscount <- ensemble
    enspresence <- ensemble
#    writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
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
    cat(paste("\n", "Start of modelling for organism: ", RASTER.species.orig, "\n\n", sep = ""))
    if (MAXENT > 0) {
        # Put the file 'maxent.jar' in the 'java' folder of dismo
        # the file 'maxent.jar' can be obtained from from http://www.cs.princeton.edu/~schapire/maxent/.
        jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
        eval1 <- eval2 <- results <- pmaxent <- NULL
        if(is.null(MAXENT.OLD) == T) {
            tryCatch(results <- maxent(x=x, p=p, a=a, factors=factors),
                error= function(err) {print(paste("MAXENT calibration failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["MAXENT"] <- -1
            }
        }else{ results <- MAXENT.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with rasterStack: ", x@title, " and calibration data for MAXENT", "\n\n", sep = ""))
            eval1 <- evaluate(p=p, a=a, model=results, x=x)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "MAXENT evaluation with test data","\n\n", sep = ""))
                tryCatch(eval2 <- evaluate(p=pt, a=at, model=results, x=x),
                    error= function(err) {print(paste("MAXENT evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                }
            }
            fullname <- paste("models/", RASTER.species.name, "_MAXENT", sep="")
            tryCatch(pmaxent <- raster::predict(object=results, x=xn, ext=ext, na.rm=TRUE, 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("MAXENT prediction failed"))},
                silent=T)
            if (is.null(pmaxent) == F) {
                pmaxent <- trunc(1000*pmaxent)
                writeRaster(x=pmaxent, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                if (models.keep==T) {ensemble.models$MAXENT <- results}
                wsa <- as.numeric(ws["MAXENT"])
                ensemble <- ensemble + wsa * pmaxent
                writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                thresholdvalue <- threshold(eval1, threshold.method) * 1000        
                pmaxent <- pmaxent > thresholdvalue
                enscount <- enscount + pmaxent
                writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
                if (evaluation.strip == T) {
                    ensemble.models$evaluation.strip$MAXENT <- as.numeric(dismo::predict(results, ResponseData))
                    ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + wsa * ensemble.models$evaluation.strip$MAXENT            
                }
            }else{
                cat(paste("\n", "WARNING: MAXENT prediction failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["MAXENT"] <- -1 
            }
        }
    }
    if (GBM > 0) {
#       require(gbm, quietly=T) 
        eval1 <- eval2 <- results <- pgbm <- NULL
        if(is.null(GBM.OLD) == T) {       
            tryCatch(results <- gbm(formula=GBM.formula, data=TrainData, weights=Yweights, distribution = "bernoulli", 
                interaction.depth = 7, shrinkage = 0.001, bag.fraction = 0.5, train.fraction = 1, 
                n.trees = GBM.n.trees, verbose = F, cv.folds = 5),
            error= function(err) {print(paste("GBM calibration failed"))},
            silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["GBM"] <- -1
            }
        }else{ results <- GBM.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with rasterStack: ", x@title, " and calibration data for GBM",  "\n\n", sep = ""))
            eval1 <- evaluate(p=p, a=a, model=results, x=x, n.trees=results$n.trees, type="response")
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "GBM evaluation with test data","\n\n", sep = ""))
                tryCatch(eval2 <- evaluate(p=pt, a=at, model=results, x=x, n.trees=results$n.trees, type="response"),
                    error= function(err) {print(paste("GBM evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                }
            }
            fullname <- paste("models/", RASTER.species.name, "_GBM", sep="")
            tryCatch(pgbm <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, n.trees=results$n.trees, type="response", 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("GBM prediction failed"))},
                silent=T)
            if (is.null(pgbm) == F) {
                pgbm <- trunc(1000*pgbm)
                writeRaster(x=pgbm, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                if (models.keep==T) {ensemble.models$GBM <- results}
                wsb <- as.numeric(ws["GBM"])
                ensemble <- ensemble + wsb * pgbm
                writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                thresholdvalue <- threshold(eval1, threshold.method) * 1000       
                pgbm <- pgbm > thresholdvalue
                enscount <- enscount + pgbm
                writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
                if (evaluation.strip == T) {
                    ensemble.models$evaluation.strip$GBM <- as.numeric(predict(results, newdata=ResponseData, n.trees=results$n.trees, type="response"))
                    ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + wsb * ensemble.models$evaluation.strip$GBM
                }
            }else{
                cat(paste("\n", "WARNING: GBM prediction failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["GBM"] <- -1 
            }
        }    
#       detach(package:gbm)
    }
    if (GBMSTEP > 0) {
#       require(gbm, quietly=T) 
        eval1 <- eval2 <- results <- pgbms <- NULL
        if(is.null(GBMSTEP.OLD) == T) {       
            tryCatch(results <- gbm.step(data=TrainData, gbm.y=1, gbm.x=GBMSTEP.gbm.x, family="bernoulli",
                site.weights=Yweights, tree.complexity=GBMSTEP.tree.complexity, learning.rate = GBMSTEP.learning.rate, 
                bag.fraction=GBMSTEP.bag.fraction, GBMSTEP.step.size=GBMSTEP.step.size, verbose=F, silent=T, plot.main=F), 
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
            cat(paste("\n", RASTER.species.orig, ": Evaluation with rasterStack: ", x@title, " and calibration data for stepwise GBM",  "\n\n", sep = ""))
            eval1 <- evaluate(p=p, a=a, model=results, x=x, n.trees=results$n.trees, type="response")
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "stepwise GBM evaluation with test data","\n\n", sep = ""))
                tryCatch(eval2 <- evaluate(p=pt, a=at, model=results, x=x, n.trees=results$n.trees, type="response"),
                    error= function(err) {print(paste("stepwise GBM evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                }
            }
            fullname <- paste("models/", RASTER.species.name, "_GBMSTEP", sep="")
            tryCatch(pgbms <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, n.trees=results$n.trees, type="response", 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("stepwise GBM prediction failed"))},
                silent=T)
            if (is.null(pgbms) == F) {
                pgbms <- trunc(1000*pgbms)
                writeRaster(x=pgbms, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                if (models.keep==T) {ensemble.models$GBMSTEP <- results}
                wsc <- as.numeric(ws["GBMSTEP"])
                ensemble <- ensemble + wsc * pgbms
                writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                thresholdvalue <- threshold(eval1, threshold.method) * 1000        
                pgbms <- pgbms > thresholdvalue
                enscount <- enscount + pgbms
                writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
                if (evaluation.strip == T) {
                    ensemble.models$evaluation.strip$GBMSTEP <- as.numeric(predict(results, newdata=ResponseData, n.trees=results$n.trees, type="response"))
                    ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + wsc * ensemble.models$evaluation.strip$GBMSTEP
                }
            }else{
                cat(paste("\n", "WARNING: stepwise GBM prediction failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["GBMSTEP"] <- -1 
            }
        }
#       detach(package:gbm)
    }
    if (RF > 0) {
#       require(randomForest, quietly=T)
        eval1 <- eval2 <- results <- prf <- TestPres <- TestAbs <- NULL
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
            cat(paste("\n", RASTER.species.orig, ": Evaluation with rasterStack: ", x@title, " and calibration data for RF",  "\n\n", sep = ""))
#   treat factors as factors during evaluations
            TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="response")
            TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="response")
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "RF evaluation with test data","\n\n", sep = ""))
                tryCatch(TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="response"),
                    error= function(err) {print(paste("RF evaluation (step 1) failed"))},
                    silent=T)
                tryCatch(TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="response"),
                    error= function(err) {print(paste("RF evaluation (step 2) failed"))},
                    silent=T)
                if (is.null(TestPres) == F && is.null(TestAbs) == F) {
                    eval2 <-  evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
                }
            }
            fullname <- paste("models/", RASTER.species.name, "_RF", sep="")
            tryCatch(prf <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("random forest prediction failed"))},
                silent=T)
            if (is.null(prf) == F) {
                prf <- trunc(1000*prf)
                writeRaster(x=prf, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                if (models.keep==T) {ensemble.models$RF <- results}
                wsd <- as.numeric(ws["RF"])
                ensemble <- ensemble + wsd * prf
                writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                thresholdvalue <- threshold(eval1, threshold.method) * 1000        
                prf <- prf > thresholdvalue
                enscount <- enscount + prf
                writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
                if (evaluation.strip == T) {
                    ensemble.models$evaluation.strip$RF <- as.numeric(predict(results, newdata=ResponseData, type="response"))
                    ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + wsd * ensemble.models$evaluation.strip$RF
                }
            }else{
                cat(paste("\n", "WARNING: random forest prediction failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["RF"] <- -1 
            }
        }
#       detach(package:randomForest)
    } 
    if (GLM > 0) {
        eval1 <- eval2 <- results <- pglm <- TestPres <- TestAbs <- NULL
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
            cat(paste("\n", RASTER.species.orig, ": Evaluation with rasterStack: ", x@title, " and calibration data for GLM",  "\n\n", sep = ""))
#           treat factors as factors during evaluations
            TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="response")
            TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="response")
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "GLM evaluation with test data","\n\n", sep = ""))
                tryCatch(TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="response"),
                    error= function(err) {print(paste("GLM evaluation (step 1) failed"))},
                    silent=T)
                tryCatch(TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="response"),
                    error= function(err) {print(paste("GLM evaluation (step 2) failed"))},
                    silent=T)
                if (is.null(TestPres) == F && is.null(TestAbs) == F) {
                    eval2 <-  evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
                }
            }
            fullname <- paste("models/", RASTER.species.name, "_GLM", sep="")
            tryCatch(pglm <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, type="response", 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("GLM prediction failed"))},
                silent=T)
            if (is.null(pglm) == F) {
                pglm <- trunc(1000*pglm)
                writeRaster(x=pglm, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                if (models.keep==T) {ensemble.models$GLM <- results}
                wse <- as.numeric(ws["GLM"])
                ensemble <- ensemble + wse * pglm
                writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                thresholdvalue <- threshold(eval1, threshold.method) * 1000        
                pglm <- pglm > thresholdvalue
                enscount <- enscount + pglm
                writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
                if (evaluation.strip == T) {
                    ensemble.models$evaluation.strip$GLM <- as.numeric(predict(results, newdata=ResponseData, type="response"))
                    ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + wse * ensemble.models$evaluation.strip$GLM
                }
            }else{
                cat(paste("\n", "WARNING: GLM prediction failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["GLM"] <- -1 
            }
        }
    }
    if (GLMSTEP > 0) {
#        require(MASS, quietly=T)
        eval1 <- eval2 <- results <- results2 <- pglms <- TestPres <- TestAbs <- NULL
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
            cat(paste("\n", RASTER.species.orig, ": Evaluation with rasterStack: ", x@title, " and calibration data for stepwise GLM", "\n\n", sep = ""))
            TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="response")
            TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="response")
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "stepwise GLM evaluation with test data","\n\n", sep = ""))
                tryCatch(TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="response"),
                    error= function(err) {print(paste("stepwise GLM evaluation (step 1) failed"))},
                    silent=T)
                tryCatch(TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="response"),
                    error= function(err) {print(paste("stepwise GLM evaluation (step 2) failed"))},
                    silent=T)
                if (is.null(TestPres) == F && is.null(TestAbs) == F) {
                    eval2 <-  evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
                }
            }
            fullname <- paste("models/", RASTER.species.name, "_GLMSTEP", sep="")
            tryCatch(pglms <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, type="response", 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("stepwise GLM prediction failed"))},
                silent=T)
            if (is.null(pglms) == F) {
                pglms <- trunc(1000*pglms)
                writeRaster(x=pglms, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                if (models.keep==T) {ensemble.models$GLMSTEP <- results}
                wsf <- as.numeric(ws["GLMSTEP"])
                ensemble <- ensemble + wsf * pglms
                writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                thresholdvalue <- eval1@t[which.max(eval1@TPR + eval1@TNR)] * 1000        
                pglms <- pglms > thresholdvalue
                enscount <- enscount + pglms
                writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
                if (evaluation.strip == T) {
                    ensemble.models$evaluation.strip$GLMSTEP <- as.numeric(predict(results, newdata=ResponseData, type="response"))
                    ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + wsf * ensemble.models$evaluation.strip$GLMSTEP
                }
            }else{
                cat(paste("\n", "WARNING: stepwise GLM prediction failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["GLMSTEP"] <- -1 
            } 
        }
#       detach(package:MASS)
    }
    if (GAM > 0) {
        require(gam, quietly=T)
        eval1 <- eval2 <- results <- pgam <- TestPres <- TestAbs <- NULL
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
            cat(paste("\n", RASTER.species.orig, ": Evaluation with rasterStack: ", x@title, " and calibration data for GAM (gam package)", "\n\n", sep = ""))
            TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="response")
            TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="response")
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "GAM evaluation with test data","\n\n", sep = ""))
                tryCatch(TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="response"),
                    error= function(err) {print(paste("GAM evaluation (step 1, gam package) failed"))},
                    silent=T)
                tryCatch(TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="response"),
                    error= function(err) {print(paste("GAM evaluation (step 2, gam package) failed"))},
                    silent=T)
                if (is.null(TestPres) == F && is.null(TestAbs) == F) {
                    eval2 <-  evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
                }
            }
            fullname <- paste("models/", RASTER.species.name, "_GAM", sep="")
            tryCatch(pgam <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, type="response", 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("GAM prediction (gam package) failed"))},
                silent=T)
            if (is.null(pgam) == F) {
                pgam <- trunc(1000*pgam)
                writeRaster(x=pgam, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                if (models.keep==T) {ensemble.models$GAM <- results}
                wsg <- as.numeric(ws["GAM"])
                ensemble <- ensemble + wsg * pgam
                writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                thresholdvalue <- threshold(eval1, threshold.method) * 1000       
                pgam <- pgam > thresholdvalue
                enscount <- enscount + pgam
                writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
                if (evaluation.strip == T) {
                    ensemble.models$evaluation.strip$GAM <- as.numeric(predict(results, newdata=ResponseData, type="response"))
                    ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + wsg * ensemble.models$evaluation.strip$GAM
               }
            }else{
                cat(paste("\n", "WARNING: GAM prediction (gam package) failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["GAM"] <- -1 
            } 
        }
        detach(package:gam)
    }
    if (GAMSTEP > 0) {
        cat(paste("\n", "NOTE: stepwise GAM seems to require assignments of family and data to pos 1", "\n\n", sep=""))
    } 
    if (GAMSTEP > 0) {
        require(gam, quietly=T)
        eval1 <- eval2 <- results <- results2 <- pgams <- TestPres <- TestAbs <- NULL
        if(is.null(GAMSTEP.OLD) == T) {
            tryCatch(results <- gam(formula=STEP.formula, family=GAM.family, data=TrainData, weights=Yweights, control=gam.control(maxit=maxit, bf.maxit=50)),
                error= function(err) {print(paste("first step of stepwise GAM calibration (gam package) failed"))},
                silent=T)
            tryCatch(results2 <- step.gam(results, scope=GAMSTEP.scope, direction="both", trace=F, steps=GAMSTEP.steps),
                error= function(err) {print(paste("stepwise GAM calibration (gam package) failed"))},
                silent=T)
            if(is.null(results2) == T) {
                prediction.failures <- TRUE
                ws["GAMSTEP"] <- -1
            }
        }else{ results2 <- GAMSTEP.OLD}
        if (is.null(results2) == F) {
            results <- results2
            cat(paste("\n", RASTER.species.orig, ": Evaluation with rasterStack: ", x@title, " and calibration data for stepwise GAM (gam package)", "\n\n", sep = ""))
            TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="response")
            TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="response")
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "stepwise GAM evaluation with test data","\n\n", "\n", sep = ""))
                tryCatch(TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="response"),
                    error= function(err) {print(paste("stepwise GAM evaluation (step 1, gam package) failed"))},
                    silent=T)
                tryCatch(TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="response"),
                    error= function(err) {print(paste("stepwise GAM evaluation (step 2, gam package) failed"))},
                    silent=T)
                if (is.null(TestPres) == F && is.null(TestAbs) == F) {
                    eval2 <-  evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
                }
            }
            fullname <- paste("models/", RASTER.species.name, "_GAMSTEP", sep="")
            tryCatch(pgams <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, type="response", 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("stepwise GAM prediction (gam package) failed"))},
                silent=T)
            if (is.null(pgams) == F) {
                pgams <- trunc(1000*pgams)
                writeRaster(x=pgams, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                if (models.keep==T) {ensemble.models$GAMSTEP <- results}
                wsh <- as.numeric(ws["GAMSTEP"])
                ensemble <- ensemble + wsh * pgams
                writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                thresholdvalue <- threshold(eval1, threshold.method) * 1000        
                pgams <- pgams > thresholdvalue
                enscount <- enscount + pgams
                writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
                if (evaluation.strip == T) {
                    ensemble.models$evaluation.strip$GAMSTEP <- as.numeric(predict(results, newdata=ResponseData, type="response"))
                    ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + wsh * ensemble.models$evaluation.strip$GAMSTEP
                }
            }else{
                cat(paste("\n", "WARNING: stepwise GAM prediction (gam package) failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["GAMSTEP"] <- -1 
            } 
        }
        detach(package:gam)
    }
    if (MGCV > 0) {
        eval1 <- eval2 <- results <- pmgcv <- NULL
        require(mgcv, quietly=T)
        if(is.null(MGCV.OLD) == T) {
            tryCatch(results <- gam(formula=MGCV.formula, family=GAM.family, data=TrainData, weights=Yweights, na.action=na.omit, 
                select=MGCV.select, control=gam.control(maxit=maxit)),
                error= function(err) {print(paste("GAM calibration (mgcv package) failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["MGCV"] <- -1
            }
        }else{ results <- MGCV.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with rasterStack: ", x@title, " and calibration data for GAM (mgcv package)", "\n\n", sep = ""))
            eval1 <- evaluate(p=p, a=a, model=results, x=x, type="response")
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "MGCV evaluation with test data","\n\n", sep = ""))
                tryCatch(eval2 <- evaluate(p=pt, a=at, model=results, x=x, type="response", na.action=na.omit),
                    error= function(err) {print(paste("GAM evaluation (mgcv package) failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                }
            }
            fullname <- paste("models/", RASTER.species.name, "_MGCV", sep="")
            tryCatch(pmgcv <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, type="response", 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("GAM prediction (mgcv package) failed"))},
                silent=T)
            if (is.null(pmgcv) == F) {
                pmgcv <- trunc(1000*pmgcv)
                writeRaster(x=pmgcv, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                if (models.keep==T) {ensemble.models$MGCV <- results}
                wsi <- as.numeric(ws["MGCV"])
                ensemble <- ensemble + wsi * pmgcv
                writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                thresholdvalue <- threshold(eval1, threshold.method) * 1000        
                pmgcv <- pmgcv > thresholdvalue
                enscount <- enscount + pmgcv
                writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
                if (evaluation.strip == T) {
                    ensemble.models$evaluation.strip$MGCV <- as.numeric(predict(results, newdata=ResponseData, type="response"))
                    ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + wsi * ensemble.models$evaluation.strip$MGCV
                }
            }else{
                cat(paste("\n", "WARNING: GAM prediction (mgcv package) failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["MGCV"] <- -1 
            } 
        }
        detach(package:mgcv)
    }
    if (MGCVFIX > 0) {
        eval1 <- eval2 <- results <- pmgcvF <- NULL
        require(mgcv, quietly=T)
        if(is.null(MGCVFIX.OLD) == T) {
            tryCatch(results <- gam(formula=MGCVFIX.formula, family=GAM.family, data=TrainData, weights=Yweights, na.action=na.omit, 
                select=FALSE, control=gam.control(maxit=maxit)),
                error= function(err) {print(paste("GAM calibration with fixed d.f. regression splines (mgcv package) failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["MGCVFIX"] <- -1
            }
        }else{ results <- MGCVFIX.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with rasterStack: ", x@title, " and calibration data for MGCVFIX (mgcv package)", "\n\n", sep = ""))
            eval1 <- evaluate(p=p, a=a, model=results, x=x, type="response")
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "MGCVFIX evaluation with test data","\n\n", sep = ""))
                tryCatch(eval2 <- evaluate(p=pt, a=at, model=results, x=x, type="response", na.action=na.omit),
                    error= function(err) {print(paste("MGCVFIX evaluation (mgcv package) failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                }
            }
            fullname <- paste("models/", RASTER.species.name, "_MGCVFIX", sep="")
            tryCatch(pmgcvf <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, type="response", 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("MGCVFIX prediction (mgcv package) failed"))},
                silent=T)
            if (is.null(pmgcvf) == F) {
                pmgcvf <- trunc(1000*pmgcvf)
                writeRaster(x=pmgcvf, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                if (models.keep==T) {ensemble.models$MGCVFIX <- results}
                wsj <- as.numeric(ws["MGCVFIX"])
                ensemble <- ensemble + wsj * pmgcvf
                writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                thresholdvalue <- threshold(eval1, threshold.method) * 1000        
                pmgcvf <- pmgcvf > thresholdvalue
                enscount <- enscount + pmgcvf
                writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
                if (evaluation.strip == T) {
                    ensemble.models$evaluation.strip$MGCVFIX <- as.numeric(predict(results, newdata=ResponseData, type="response"))
                    ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + wsj * ensemble.models$evaluation.strip$MGCVFIX
                }
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
    if (EARTH > 0) {
#        require(earth, quietly=T)
        eval1 <- eval2 <- results <- pearth <- NULL
        if(is.null(EARTH.OLD) == T) {
            tryCatch(results <- earth(formula=EARTH.formula, glm=EARTH.glm, data=TrainData, weights=Yweights, degree=2),
                error= function(err) {print(paste("MARS calibration (earth package) failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["EARTH"] <- -1
            }
        }else{ results <- EARTH.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with rasterStack: ", x@title, " and calibration data for MARS (earth package)", "\n\n", sep = ""))
            eval1 <- evaluate(p=p, a=a, model=results, x=x, type="response")
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "MARS evaluation with test data","\n\n", sep = ""))
                tryCatch(eval2 <- evaluate(p=pt, a=at, model=results, x=x, type="response"),
                    error= function(err) {print(paste("MARS evaluation (earth package) failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                }
            }
            fullname <- paste("models/", RASTER.species.name, "_EARTH", sep="")
            tryCatch(pearth <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, type="response", 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("MARS prediction (earth package) failed"))},
                silent=T)
            if (is.null(pearth) == F) {
                pearth <- trunc(1000*pearth)
                writeRaster(x=pearth, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                if (models.keep==T) {ensemble.models$EARTH <- results}
                wsk <- as.numeric(ws["EARTH"])
                ensemble <- ensemble + wsk * pearth
                writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                thresholdvalue <- threshold(eval1, threshold.method) * 1000        
                pearth <- pearth > thresholdvalue
                enscount <- enscount + pearth
                writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
                if (evaluation.strip == T) {
                    ensemble.models$evaluation.strip$EARTH <- as.numeric(predict(results, newdata=ResponseData, type="response"))
                    ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + wsk * ensemble.models$evaluation.strip$EARTH
                }
            }else{
                cat(paste("\n", "WARNING: MARS prediction (earth package) failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["EARTH"] <- -1 
            } 
        }
#       detach(package:earth)
    }
    if (RPART > 0) {
#        require(rpart, quietly=T)
        eval1 <- eval2 <- results <- prpart <- TestPres <- TestAbs <- NULL
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
            cat(paste("\n", RASTER.species.orig, ": Evaluation with rasterStack: ", x@title, " and calibration data for RPART", "\n\n", sep = "")) 
            TrainPres <- predict(results, newdata=TrainData[TrainData[,"pb"]==1,], type="prob")[,2]
            TrainAbs <- predict(results, newdata=TrainData[TrainData[,"pb"]==0,], type="prob")[,2]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "RPART evaluation with test data","\n\n", sep = ""))
                tryCatch(TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="prob")[,2],
                    error= function(err) {print(paste("RPART evaluation (step 1) failed"))},
                    silent=T)
                tryCatch(TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="prob")[,2],
                    error= function(err) {print(paste("RPART evaluation (step 2) failed"))},
                    silent=T)
                if (is.null(TestPres) == F && is.null(TestAbs) == F) {
                    eval2 <-  evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
                }
            }
            fullname <- paste("models/", RASTER.species.name, "_RPART", sep="")
            tryCatch(prpart <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, type="prob", index=2,
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("RPART prediction failed"))},
                silent=T)
            if (is.null(prpart) == F) {
                prpart <- trunc(1000*prpart)
                writeRaster(x=prpart, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                if (models.keep==T) {ensemble.models$RPART <- results}
                wsl <- as.numeric(ws["RPART"])
                ensemble <- ensemble + wsl * prpart
                writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                thresholdvalue <- threshold(eval1, threshold.method) * 1000        
                prpart <- prpart > thresholdvalue
                enscount <- enscount + prpart
                writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
                if (evaluation.strip == T) {
                    ensemble.models$evaluation.strip$RPART <- as.numeric(predict(results, newdata=ResponseData, type="prob")[,2])
                    ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + wsl * ensemble.models$evaluation.strip$RPART
                }
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
    if (NNET > 0) {
#        require(nnet, quietly=T)
        eval1 <- eval2 <- results <- pnnet <- NULL
        if(is.null(NNET.OLD) == T) {
            tryCatch(results <- nnet(formula=NNET.formula, size=NNET.size, decay=NNET.decay, data=TrainData, weights=Yweights, 
                rang=0.1, maxit=maxit, trace=F, na.action=na.omit),
                error= function(err) {print(paste("ANN calibration (nnet package) failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["NNET"] <- -1
            }
        }else{ results <- NNET.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with rasterStack: ", x@title, " and calibration data for ANN (nnet package)", "\n\n", sep = ""))
            eval1 <- evaluate(p=p, a=a, model=results, x=x, type="raw")
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "ANN evaluation with test data","\n\n", sep = ""))
                tryCatch(eval2 <- evaluate(p=pt, a=at, model=results, x=x, type="raw"),
                    error= function(err) {print(paste("ANN evaluation (nnet package) failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                }
            }
            fullname <- paste("models/", RASTER.species.name, "_NNET", sep="")
            tryCatch(pnnet <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, type="raw", 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("ANN prediction (nnet package) failed"))},
                silent=T)
            if (is.null(pnnet) == F) {
                pnnet <- trunc(1000*pnnet)
                writeRaster(x=pnnet, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                if (models.keep==T) {ensemble.models$NNET <- results}
                wsm <- as.numeric(ws["NNET"])
                ensemble <- ensemble + wsm * pnnet
                writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                thresholdvalue <- threshold(eval1, threshold.method) * 1000        
                pnnet <- pnnet > thresholdvalue
                enscount <- enscount + pnnet
                writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
                if (evaluation.strip == T) {
                    ensemble.models$evaluation.strip$NNET <- as.numeric(predict(results, newdata=ResponseData, type="raw"))
                    ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + wsm * ensemble.models$evaluation.strip$NNET
                }
            }else{
                cat(paste("\n", "WARNING: ANN prediction (nnet package) failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["NNET"] <- -1 
            } 
        }
#       detach(package:nnet)
    }
    if (FDA > 0) {
#        require(mda, quietly=T)
        eval1 <- eval2 <- results <- pfda <- TestPres <- TestAbs <- NULL
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
            cat(paste("\n", RASTER.species.orig, ": Evaluation with rasterStack: ", x@title, " and calibration data for FDA", "\n\n", sep = "")) 
            TrainPres <- predict(results,newdata=TrainData[TrainData[,"pb"]==1,], type="posterior")[,2]
            TrainAbs <- predict(results,newdata=TrainData[TrainData[,"pb"]==0,], type="posterior")[,2]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "FDA evaluation with test data","\n\n", sep = ""))
                tryCatch(TestPres <- predict(results,newdata=TestData[TestData[,"pb"]==1,], type="posterior")[,2],
                    error= function(err) {print(paste("FDA evaluation (step 1) failed"))},
                    silent=T)
                tryCatch(TestAbs <- predict(results,newdata=TestData[TestData[,"pb"]==0,], type="posterior")[,2],
                    error= function(err) {print(paste("FDA evaluation (step 2) failed"))},
                    silent=T)
                if (is.null(TestPres) == F && is.null(TestAbs) == F) {
                    eval2 <-  evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
                }
            }
            fullname <- paste("models/", RASTER.species.name, "_FDA", sep="")
            tryCatch(pfda <- raster::predict(object=xn, model=results, ext=ext, na.rm=TRUE, type="posterior", index=2,
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("FDA prediction failed"))},
                silent=T)
            if (is.null(pfda) == F) {
                pfda <- trunc(1000*pfda)
                writeRaster(x=pfda, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                if (models.keep==T) {ensemble.models$FDA <- results}
                wsn <- as.numeric(ws["FDA"])
                ensemble <- ensemble + wsn * pfda
                writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                thresholdvalue <- threshold(eval1, threshold.method) * 1000        
                pfda <- pfda > thresholdvalue
                enscount <- enscount + pfda
                writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
                if (evaluation.strip == T) {
                    ensemble.models$evaluation.strip$FDA <- as.numeric(predict(results, newdata=ResponseData, type="posterior")[,2])
                    ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + wsn * ensemble.models$evaluation.strip$FDA
                }
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
    if (SVM > 0) {
#        require(kernlab, quietly=T)
        eval1 <- eval2 <- results <- psvm <- TestPres <- TestAbs <- NULL
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
            cat(paste("\n", RASTER.species.orig, ": Evaluation with rasterStack: ", x@title, " and calibration data for SVM (kernlab package)", "\n\n", sep = "")) 
            TrainPres <- kernlab::predict(results, newdata=TrainData[TrainData[,"pb"]==1,], type="probabilities")[,2]
            TrainAbs <- kernlab::predict(results, newdata=TrainData[TrainData[,"pb"]==0,], type="probabilities")[,2]
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "SVM evaluation with test data","\n\n", sep = ""))
                tryCatch(TestPres <- kernlab::predict(results, newdata=TestData[TestData[,"pb"]==1,], type="probabilities")[,2],
                    error= function(err) {print(paste("SVM evaluation (step 1, kernlab package) failed"))},
                    silent=T)
                tryCatch(TestAbs <- kernlab::predict(results, newdata=TestData[TestData[,"pb"]==0,], type="probabilities")[,2],
                    error= function(err) {print(paste("SVM evaluation (step 2, kernlab package) failed"))},
                    silent=T)
                if (is.null(TestPres) == F && is.null(TestAbs) == F) {
                    eval2 <-  evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
                }
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
                if (models.keep==T) {ensemble.models$SVM <- results}
                wso <- as.numeric(ws["SVM"])
                ensemble <- ensemble + wso * psvm
                writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                thresholdvalue <- threshold(eval1, threshold.method) * 1000        
                psvm <- psvm > thresholdvalue
                enscount <- enscount + psvm
                writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
                if (evaluation.strip == T) {
                    ensemble.models$evaluation.strip$SVM <- as.numeric(kernlab::predict(results, newdata=ResponseData, type="probabilities")[,2])
                    ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + wso * ensemble.models$evaluation.strip$SVM
                }
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
    if (SVME > 0) {
#        require(e1071, quietly=T)
        eval1 <- eval2 <- results <- psvme <- TestPres <- TestAbs <- NULL
        if(is.null(SVME.OLD) == T) {
            tryCatch(results <- svm(SVME.formula, data=TrainData, type="C-classification", kernel="polynomial", degree=3, probability=T),
                error= function(err) {print(paste("SVM calibration (e1071 package) failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["SVME"] <- -1
            }
        }else{ results <- SVME.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with rasterStack: ", x@title, " and calibration data for SVM (e1071 package)", "\n\n", sep = ""))
            TrainPres <- predict.svme(results, newdata=TrainData[TrainData[,"pb"]==1,])
            TrainAbs <- predict.svme(results, newdata=TrainData[TrainData[,"pb"]==0,])
            eval1 <- evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "SVME evaluation with test data","\n\n", sep = ""))
                tryCatch(TestPres <- predict.svme(results, newdata=TestData[TestData[,"pb"]==1,]),
                    error= function(err) {print(paste("SVM evaluation (step 1, e1071 package) failed"))},
                    silent=T)
                tryCatch(TestAbs <- predict.svme(results, newdata=TestData[TestData[,"pb"]==0,]),
                    error= function(err) {print(paste("SVM evaluation (step 2, e1071 package) failed"))},
                    silent=T)
                if (is.null(TestPres) == F && is.null(TestAbs) == F) {
                    eval2 <-  evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
                }
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
                if (models.keep==T) {ensemble.models$SVME <- results}
                wsp <- as.numeric(ws["SVME"])
                ensemble <- ensemble + wsp * psvme
                writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                thresholdvalue <- threshold(eval1, threshold.method) * 1000        
                psvme <- psvme > thresholdvalue
                enscount <- enscount + psvme
                writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
                if (evaluation.strip == T) {
                    ensemble.models$evaluation.strip$SVME <- as.numeric(predict.svme(results, newdata=ResponseData))
                    ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + wsp * ensemble.models$evaluation.strip$SVME
                }
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
                x <- dropLayer(x, which(names(x) == factors[i]))
                xn <- dropLayer(xn, which(names(xn) == factors[i]))                
            }
        }
    }
    if (BIOCLIM > 0) {  
        eval1 <- eval2 <- results <- pbio <- NULL  
        if(is.null(BIOCLIM.OLD) == T) {
            tryCatch(results <- bioclim(x=x, p=p),
                error= function(err) {print(paste("BIOCLIM calibration failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["BIOCLIM"] <- -1
            }
            }else{ results <- BIOCLIM.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with rasterStack: ", x@title, " and calibration data for BIOCLIM", "\n\n", sep = ""))
            eval1 <- evaluate(p=p, a=a, model=results, x=x)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "BIOCLIM evaluation with test data","\n\n", sep = ""))
                tryCatch(eval2 <- evaluate(p=pt, a=at, model=results, x=x),
                    error= function(err) {print(paste("BIOCLIM evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                }
            }
            fullname <- paste("models/", RASTER.species.name, "_BIOCLIM", sep="")
            tryCatch(pbio <- dismo::predict(object=results, x=xn, ext=ext, na.rm=TRUE, 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("BIOCLIM prediction failed"))},
                silent=T)
            if (is.null(pbio) == F) {
                pbio <- trunc(1000*pbio)
                writeRaster(x=pbio, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                if (models.keep==T) {ensemble.models$BIOCLIM <- results}
                wsq <- as.numeric(ws["BIOCLIM"])
                ensemble <- ensemble + wsq * pbio
                writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                thresholdvalue <- threshold(eval1, threshold.method) * 1000        
                pbio <- pbio > thresholdvalue
                enscount <- enscount + pbio
                writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
                if (evaluation.strip == T) {
                    ensemble.models$evaluation.strip$BIOCLIM <- as.numeric(dismo::predict(results, ResponseData))
                    ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + wsq * ensemble.models$evaluation.strip$BIOCLIM
                }
            }else{
                cat(paste("\n", "WARNING: BIOCLIM prediction failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["BIOCLIM"] <- -1 
            }
        }
    }
    if (DOMAIN > 0) {
        eval1 <- eval2 <- results <- pdom <- NULL
        if(is.null(DOMAIN.OLD) == T) {
            tryCatch(results <- domain(x=x, p=p, factors=factors),
                error= function(err) {print(paste("DOMAIN calibration failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["DOMAIN"] <- -1
            }
        }else{ results <- DOMAIN.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with rasterStack: ", x@title, " and calibration data for DOMAIN", "\n\n", sep = ""))
            eval1 <- evaluate(p=p, a=a, model=results, x=x)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "DOMAIN evaluation with test data","\n\n", sep = ""))
                tryCatch(eval2 <- evaluate(p=pt, a=at, model=results, x=x),
                    error= function(err) {print(paste("DOMAIN evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                }
            }
            fullname <- paste("models/", RASTER.species.name, "_DOMAIN", sep="")
            tryCatch(pdom <- dismo::predict(object=results, x=xn, ext=ext, na.rm=TRUE, 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("DOMAIN prediction failed"))},
                silent=T)
            if (is.null(pdom) == F) {
                pdom <- trunc(1000*pdom)
                writeRaster(x=pdom, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                if (models.keep==T) {ensemble.models$DOMAIN <- results}
                wsr <- as.numeric(ws["DOMAIN"])
                ensemble <- ensemble + wsr * pdom
                writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                thresholdvalue <- threshold(eval1, threshold.method) * 1000        
                pdom <- pdom > thresholdvalue
                enscount <- enscount + pdom
                writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
                if (evaluation.strip == T) {
                    ensemble.models$evaluation.strip$DOMAIN <- as.numeric(dismo::predict(results, ResponseData))
                    ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + wsr * ensemble.models$evaluation.strip$DOMAIN
                }
            }else{
                cat(paste("\n", "WARNING: DOMAIN prediction failed","\n\n", sep = ""))
                prediction.failures <- TRUE
                ws["DOMAIN"] <- -1 
            }
        }
    }
    if (MAHAL > 0) {  
        eval1 <- eval2 <- results <- pmahal <- NULL  
        if(is.null(MAHAL.OLD) == T) {
            tryCatch(results <- mahal(x=x, p=p),
                error= function(err) {print(paste("Mahalanobis calibration failed"))},
                silent=T)
            if(is.null(results) == T) {
                prediction.failures <- TRUE
                ws["MAHAL"] <- -1
            }
        }else{ results <- MAHAL.OLD}
        if (is.null(results) == F) {
            cat(paste("\n", RASTER.species.orig, ": Evaluation with rasterStack: ", x@title, " and calibration data for Mahalanobis", "\n\n", sep = ""))
            eval1 <- evaluate(p=p, a=a, model=results, x=x)
            print(eval1)
            if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
                cat(paste("\n", "Mahalanobis evaluation with test data","\n\n", sep = ""))
                tryCatch(eval2 <- evaluate(p=pt, a=at, model=results, x=x),
                    error= function(err) {print(paste("Mahalanobis evaluation failed"))},
                    silent=T)
                if (is.null(eval2) == F) {
                    print(eval2)
                }
            }
            fullname <- paste("models/", RASTER.species.name, "_MAHAL", sep="")
# first calculate presence
            tryCatch(pmahal <- dismo::predict(object=results, x=xn, ext=ext, na.rm=TRUE, 
                    filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format),
                error= function(err) {print(paste("Mahalanobis prediction failed"))},
                silent=T)
            if (is.null(pmahal) == F) {
                thresholdvalue <- threshold(eval1, threshold.method)        
                pmahal2 <- pmahal > thresholdvalue
                enscount <- enscount + pmahal2
                writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
#                setMinMax(pmahal)
#                mahalmin <- minValue(pmahal)
#                mahalrange <- (1 - mahalmin)
#                pmahal <- pmahal - mahalmin
#                pmahal <- pmahal / mahalrange
#                pmahal <- trunc(1000*pmahal)
                pmahal <- pmahal - 2
                pmahal <- abs(pmahal)
                pmahal <- 1 / pmahal
                pmahal <- trunc(1000*pmahal)
                writeRaster(x=pmahal, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                if (models.keep==T) {ensemble.models$MAHAL <- results}
                wss <- as.numeric(ws["MAHAL"])
                ensemble <- ensemble + wss * pmahal
                writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
                if (evaluation.strip == T) {
                    ensemble.models$evaluation.strip$MAHAL <- as.numeric(dismo::predict(results, ResponseData))
                    mahalmin <- min(ensemble.models$evaluation.strip$MAHAL)
                    mahalrange <- (1 - mahalmin)
                    ensemble.models$evaluation.strip$MAHAL <- ensemble.models$evaluation.strip$MAHAL - mahalmin
                    ensemble.models$evaluation.strip$MAHAL <- ensemble.models$evaluation.strip$MAHAL / mahalrange
                    ensemble.models$evaluation.strip$ENSEMBLE <- ensemble.models$evaluation.strip$ENSEMBLE + wss * ensemble.models$evaluation.strip$MAHAL
                }
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
    }
    ensemble <- trunc(ensemble)
    writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
    cat(paste("\n", RASTER.species.orig, ": Ensemble evaluation with calibration data", "\n\n", sep = ""))
    pres_consensus <- extract(ensemble, p)
    abs_consensus <- extract(ensemble, a)
    eval1 <- evaluate(p=pres_consensus, a=abs_consensus)
    print(eval1)
    l1 <- minValue(ensemble)
    l5 <- maxValue(ensemble)
    l3 <- threshold(eval1, threshold.method)
    r1 <- l3 - l1
    r2 <- l5 - l3 
    l2 <- l1 + r1/2
    l4 <- l3 + r2/2
    cat(paste("\n", "Suggested thresholds for absence", sep = ""))
    cat(paste("\n", "(includes minimum, midpoint and threshold)","\n\n", sep = ""))     
    cat(paste(l1, ", ", l2, ", ", l3, "\n", sep=""))
    cat(paste("\n", "Suggested thresholds for presence", sep = ""))
    cat(paste("\n", "(includes threshold, midpoint and maximum)","\n\n", sep = ""))      
    cat(paste(l3, ", ", l4, ", ", l5, "\n\n", sep=""))
    enspresence <- ensemble > l3
    writeRaster(x=enspresence, filename=rasterpresence, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    if ((identical(pt, p) == F) || (identical(at, a) == F) == T) {
        cat(paste("\n", RASTER.species.orig, ": Ensemble evaluation with testing data", "\n\n", sep = ""))
        pres_consensus <- extract(ensemble, pt)
        abs_consensus <- extract(ensemble, at)
        eval2 <- evaluate(p=pres_consensus, a=abs_consensus)
        print(eval2)
    }
    cat(paste("\n", "End of modelling for organism: ", RASTER.species.orig, "\n\n", sep = ""))      
    remove(Yweights, envir=.BiodiversityR)
    remove(TrainData, envir=.BiodiversityR)
    remove(TestData,envir=.BiodiversityR)
    if (evaluation.strip==T) {remove(ResponseData, envir=.BiodiversityR)}
    if (models.keep==T || evaluation.strip==T) {return(ensemble.models)}
}




