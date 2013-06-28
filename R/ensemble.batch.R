`ensemble.batch` <- function(
    x=NULL, species.presence=NULL, species.absence=NULL, an=1000, excludep=FALSE, ext=NULL, 
    k.raster=0, k.splits=5, n.ensembles=1, xn=c(x), 
    TrainData=NULL, TestData=NULL,
    layer.drops=NULL,
    presence.min=20,
    RASTER.format="raster", RASTER.datatype="INT2S", RASTER.NAflag=-32767,
    threshold.method="spec_sens",
    ENSEMBLE.decay=1, ENSEMBLE.best=0, ENSEMBLE.min=0.7,
    input.weights=NULL,
    MAXENT=1, GBM=1, GBMSTEP=1, RF=1, GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=1, MGCV=1, MGCVFIX=0,
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, BIOCLIM=1, DOMAIN=1, MAHAL=1, 
    GEODIST=0,
    Yweights="BIOMOD", factors=NULL, dummy.vars=NULL, 
    formulae.defaults=TRUE, maxit=100,
    MAXENT.a=NULL, MAXENT.an=10000, MAXENT.BackData=NULL, MAXENT.path=paste(getwd(), "/models/maxent", sep=""),
    GBM.formula=NULL, GBM.n.trees=2001,
    GBMSTEP.gbm.x=2:(1+nlayers(x)), GBMSTEP.tree.complexity=5, GBMSTEP.learning.rate=0.005, 
    GBMSTEP.bag.fraction=0.5, GBMSTEP.step.size=100,
    RF.formula=NULL, RF.ntree=751, RF.mtry=floor(sqrt(nlayers(x))), 
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
    SVM.formula=NULL, SVME.formula=NULL,
    MAHAL.shape=1  
)
{
    .BiodiversityR <- new.env()
    k.raster <- as.integer(k.raster)
    k.splits <- as.integer(k.splits)
    if (k.splits < 1) {
        cat(paste("\n", "NOTE: parameter k.splits was set to be smaller than 1", sep = ""))
        cat(paste("\n", "default value of 5 therefore set for parameter k.splits", sep = ""))
        k.splits <- 5
    }
    n.ensembles <- as.integer(n.ensembles)
    if (n.ensembles < 1) {n.ensembles <- 1}
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
                for (j in 1:length(xn)) {
                    xn[[j]] <- dropLayer(xn[[j]], which(names(xn[[j]]) == layer.drops[i]))
                }
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
# set minimum and maximum values
    for (i in 1:nlayers(x)) {
        x[[i]] <- setMinMax(x[[i]])
    }
    species.presence <- data.frame(species.presence)
    if (ncol(species.presence) < 2) {stop("species.presence expected to be 3-column data.frame with species, x (e.g., lon) and y (e.g., lat) columns")}
    if (ncol(species.presence) == 2) {
        cat(paste("\n", "species.presence was expected to be 3-column data.frame with columns representing species, x (e.g., lon) and y (e.g., lat)", sep = ""))        
        cat(paste("\n", "only two columns were provided, it is therefore assumed that these reflect x and y coordinates for a single species", "\n\n", sep = ""))
        species.name <- rep("Species001", nrow(species.presence))
        species.presence <- cbind(species.name, species.presence)
        species.presence <- data.frame(species.presence)
    }
    if (ncol(species.presence) > 3) {
        cat(paste("\n", "species.presence was expected to be 3-column data.frame with species, x (e.g., lon) and y (e.g., lat) columns", sep = ""))        
        cat(paste("\n", "only first three columns used", "\n\n", sep = ""))
        species.presence <- species.presence[,c(1:3)] 
    }
    if (is.null(species.absence)==F && ncol(species.absence) < 2) {stop("species.absence expected to be a 2-column data.frame with x (e.g., lon) and y (e.g., lat),  or 3-column data.frame with species, x (e.g., lon) and y (e.g., lat) columns")}
    if (is.null(species.absence)==F && ncol(species.absence)> 3) {
        cat(paste("\n", "species.absence was expected to be 3-column data.frame with species, x (e.g., lon) and y (e.g., lat) columns", sep = ""))        
        cat(paste("\n", "only first three columns used", "\n\n", sep = ""))
        species.absence <- species.absence[,c(1:3)] 
    }
    if (is.null(species.absence)==F && ncol(species.absence) == 2) {
        cat(paste("\n", "species.absence was expected to be 3-column data.frame with species, x (e.g., lon) and y (e.g., lat) columns", sep = ""))        
        cat(paste("\n", "only two columns were provided, it is therefore assumed that these reflect x and y coordinates for absence locations to be used for each species run", "\n\n", sep = ""))
        as <- data.frame(species.absence)
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
#
# create background data for MAXENT
# same background data used for all species
    if (MAXENT > 0) {
        if (is.null(MAXENT.BackData) == T) {
            if (is.null(x) == T) {
                cat(paste("\n", "WARNING: not possible to create MAXENT.BackData as RasterStack x is missing", sep = "")) 
                cat(paste("\n", "MAXENT model will not be calibrated", "\n", sep = "")) 
                MAXENT <- 0
            }else{
# default option of MAXENT is to exclude presence locations
# not possible now as different species are modelled
                if (is.null(MAXENT.a) == T) {
                    MAXENT.a <- randomPoints(x[[1]], n=MAXENT.an, p=NULL, ext=ext, excludep=F)
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
        assign("MAXENT.BackData", MAXENT.BackData, envir=.BiodiversityR)
    }
# 
# process species by species
    species.names <- levels(factor(species.presence[,1,drop=T]))

# output only for raster evaluations
    output <- data.frame(array(dim=c(length(species.names)*20, n.ensembles*2+2), NA))
    output[,1] <- as.factor(output[,1])
    output[,2] <- as.factor(output[,2])
    colnames(output)[1] <- "species"
    colnames(output)[2] <- "model"
    if (n.ensembles > 1) {
       for (j in 1:n.ensembles) {
           colnames(output)[2*(j-1)+3] <- paste("calibration_ENSEMBLE_", j, sep="")
           colnames(output)[2*(j-1)+4] <- paste("test_ENSEMBLE_", j, sep="")
        }
    }else{
           colnames(output)[3] <- "calibration"
           colnames(output)[4] <- "test"
    }
    output[,1] <- rep(species.names, each=length(species.names))
    modelnames <- c("ENSEMBLE", "MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL")
    output[,2] <- rep(modelnames, length(species.names))    

    for (s in 1:length(species.names)) {
        focal.species <- species.names[s]

# check if species has required minimum number of presence points
        n.pres <- nrow(species.presence[species.presence[,1]==focal.species,])
        if (n.pres < presence.min) { 
            cat(paste("\n", "Species: ", focal.species, " only has ", n.pres, " presence locations", sep = ""))
            cat(paste("\n", "This species therefore not included in batch processing", "\n\n", sep = ""))

        }else{

        cat(paste("\n", "Evaluations for species: ", focal.species, "\n\n", sep = ""))
        ps <- species.presence[species.presence[,1]==focal.species, c(2:3)]
        if (is.null(species.absence)==F && ncol(species.absence) == 3) {
            as <- species.absence[species.absence[,1]==focal.species, c(2:3)]
        }
        if (is.null(species.absence)==T) {
            if (excludep == T) {
                as <- randomPoints(x[[1]], n=an, p=ps, ext=ext, excludep=T)
            }else{
                as <- randomPoints(x[[1]], n=an, p=NULL, ext=ext, excludep=F)
            }
        }


# repeat the whole process for n.ensembles

    RASTER.species.name1 <- focal.species
    for (runs in 1:n.ensembles) {
        if (n.ensembles > 1) { 
            cat(paste("\n", focal.species, ": ENSEMBLE ", runs, "\n\n", sep = ""))
            RASTER.species.name1 <- paste(focal.species, "_ENSEMBLE_", runs, sep="")
        }

#1. first ensemble tests
    splits <- ensemble.test.splits(x=x, p=ps, a=as, ext=ext, k=k.splits, 
        VIF=T,
        PLOTS=F,
        ENSEMBLE.decay=ENSEMBLE.decay, ENSEMBLE.best=ENSEMBLE.best, ENSEMBLE.min=ENSEMBLE.min, 
        input.weights=input.weights,
        MAXENT=MAXENT, GBM=GBM, GBMSTEP=GBMSTEP, RF=RF, GLM=GLM, GLMSTEP=GLMSTEP, 
        GAM=GAM, GAMSTEP=GAMSTEP, MGCV=MGCV, EARTH=EARTH, RPART=RPART, 
        NNET=NNET, FDA=FDA, SVM=SVM, SVME=SVME, BIOCLIM=BIOCLIM, DOMAIN=DOMAIN, MAHAL=MAHAL,
        GEODIST=GEODIST, GEODIST.file.name=RASTER.species.name1, RASTER.format=RASTER.format,  
        Yweights=Yweights, factors=factors, dummy.vars=dummy.vars,
        maxit=maxit,
        MAXENT.BackData=MAXENT.BackData, MAXENT.path=MAXENT.path,
        GBM.formula=GBM.formula, GBM.n.trees=GBM.n.trees,
        GBMSTEP.gbm.x=GBMSTEP.gbm.x, GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
        GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction,
        GBMSTEP.step.size=GBMSTEP.step.size,
        RF.formula=RF.formula, RF.ntree=RF.ntree, RF.mtry=RF.mtry, 
        GLM.formula=GLM.formula, GLM.family=GLM.family, 
        GLMSTEP.k=GLMSTEP.k, GLMSTEP.steps=GLMSTEP.steps, STEP.formula=STEP.formula, GLMSTEP.scope=GLMSTEP.scope, 
        GAM.formula=GAM.formula, GAM.family=GAM.family, 
        GAMSTEP.steps=GAMSTEP.steps, GAMSTEP.scope=GAMSTEP.scope, GAMSTEP.pos=GAMSTEP.pos,
        MGCV.formula=MGCV.formula, MGCV.select=MGCV.select,
        MGCVFIX.formula=MGCVFIX.formula, 
        EARTH.formula=EARTH.formula, EARTH.glm=EARTH.glm,
        RPART.formula=RPART.formula, RPART.xval=RPART.xval, 
        NNET.formula=NNET.formula, NNET.size=NNET.size, NNET.decay=NNET.decay,
        FDA.formula=FDA.formula, SVM.formula=SVM.formula, SVME.formula=SVME.formula,
        MAHAL.shape=MAHAL.shape)

#2. calibrate models and predict first raster
#    xn <- as.character(xn)    
    xn.f <- xn[[1]]
    if(length(xn.f@title) == 0) {xn.f@title <- "stack1"}
#    xn.f <- eval(as.name(xn.focal))
    cat(paste("\n", "Predictions for species: ", RASTER.species.name1, " for rasterStack: ", xn.f@title,  "\n\n", sep = ""))

    rasters <- ensemble.raster(
        x=x, p=ps, a=as, an=1000, ext=ext, k=k.raster, pt=NULL, at=NULL, xn=xn.f, 
        models.keep=TRUE, 
        RASTER.species.name=RASTER.species.name1, RASTER.stack.name=xn.f@title, RASTER.format=RASTER.format, RASTER.datatype=RASTER.datatype, RASTER.NAflag=RASTER.NAflag,
        RASTER.models.overwrite=TRUE,
        threshold.method=threshold.method,
        ENSEMBLE.decay=1, ENSEMBLE.best=0, ENSEMBLE.min=0, 
        input.weights=splits$output.weights, models.list=NULL,
        MAXENT=0, GBM=0, GBMSTEP=0, RF=0, GLM=0, GLMSTEP=0, GAM=0, GAMSTEP=0, MGCV=0, MGCVFIX=0,
        EARTH=0, RPART=0, NNET=0, FDA=0, SVM=0, SVME=0, BIOCLIM=0, DOMAIN=0, MAHAL=0,  
        Yweights=Yweights, factors=factors, dummy.vars=dummy.vars,
        evaluation.strip=FALSE, 
        maxit=maxit,
        MAXENT.BackData=MAXENT.BackData, MAXENT.OLD=NULL, MAXENT.path=MAXENT.path,
        GBM.formula=GBM.formula, GBM.n.trees=GBM.n.trees, GBM.OLD=NULL,    
        GBMSTEP.gbm.x=GBMSTEP.gbm.x, GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
        GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction, GBMSTEP.OLD=NULL,
        RF.formula=RF.formula, RF.ntree=RF.ntree, RF.mtry=RF.mtry, RF.OLD=NULL,
        GLM.formula=GLM.formula, GLM.family=GLM.family, GLM.OLD=NULL,
        GLMSTEP.k=GLMSTEP.k, GLMSTEP.steps=GLMSTEP.steps, STEP.formula=STEP.formula, GLMSTEP.scope=GLMSTEP.scope, GLMSTEP.OLD=NULL,
        GAM.formula=GAM.formula, GAM.family=GAM.family, GAM.OLD=NULL, 
        GAMSTEP.steps=GAMSTEP.steps, GAMSTEP.scope=GAMSTEP.scope, GAMSTEP.OLD=NULL, GAMSTEP.pos=GAMSTEP.pos,
        MGCV.formula=MGCV.formula, MGCV.select=MGCV.select, MGCV.OLD=NULL,
        MGCVFIX.formula=MGCVFIX.formula, MGCVFIX.OLD=NULL,
        EARTH.formula=EARTH.formula, EARTH.glm=EARTH.glm, EARTH.OLD=NULL,
        RPART.formula=RPART.formula, RPART.xval=RPART.xval, RPART.OLD=NULL,
        NNET.formula=NNET.formula, NNET.size=NNET.size, NNET.decay=NNET.decay, NNET.OLD=NULL,
        FDA.formula=FDA.formula,  FDA.OLD=NULL, 
        SVM.formula=SVM.formula, SVM.OLD=NULL, 
        SVME.formula=SVME.formula, SVME.OLD=NULL,
        BIOCLIM.OLD=NULL, DOMAIN.OLD=NULL, MAHAL.OLD=NULL,
        MAHAL.shape=MAHAL.shape) 

    output[(s-1)*20+1, (runs-1)*2+3] <- rasters$evaluations$evaluations["ENSEMBLE", 1]
    output[(s-1)*20+1, (runs-1)*2+4] <- rasters$evaluations$evaluations["ENSEMBLE", 2]
    output[(s-1)*20+2, (runs-1)*2+3] <- rasters$evaluations$evaluations["MAXENT", 1]
    output[(s-1)*20+2, (runs-1)*2+4] <- rasters$evaluations$evaluations["MAXENT", 2]
    output[(s-1)*20+3, (runs-1)*2+3] <- rasters$evaluations$evaluations["GBM", 1]
    output[(s-1)*20+3, (runs-1)*2+4] <- rasters$evaluations$evaluations["GBM", 2]
    output[(s-1)*20+4, (runs-1)*2+3] <- rasters$evaluations$evaluations["GBMSTEP", 1]
    output[(s-1)*20+4, (runs-1)*2+4] <- rasters$evaluations$evaluations["GBMSTEP", 2]
    output[(s-1)*20+5, (runs-1)*2+3] <- rasters$evaluations$evaluations["RF", 1]
    output[(s-1)*20+5, (runs-1)*2+4] <- rasters$evaluations$evaluations["RF", 2]
    output[(s-1)*20+6, (runs-1)*2+3] <- rasters$evaluations$evaluations["GLM", 1]
    output[(s-1)*20+6, (runs-1)*2+4] <- rasters$evaluations$evaluations["GLM", 2]
    output[(s-1)*20+7, (runs-1)*2+3] <- rasters$evaluations$evaluations["GLMSTEP", 1]
    output[(s-1)*20+7, (runs-1)*2+4] <- rasters$evaluations$evaluations["GLMSTEP", 2]
    output[(s-1)*20+8, (runs-1)*2+3] <- rasters$evaluations$evaluations["GAM", 1]
    output[(s-1)*20+8, (runs-1)*2+4] <- rasters$evaluations$evaluations["GAM", 2]
    output[(s-1)*20+9, (runs-1)*2+3] <- rasters$evaluations$evaluations["GAMSTEP", 1]
    output[(s-1)*20+9, (runs-1)*2+4] <- rasters$evaluations$evaluations["GAMSTEP", 2]
    output[(s-1)*20+10, (runs-1)*2+3] <- rasters$evaluations$evaluations["MGCV", 1]
    output[(s-1)*20+10, (runs-1)*2+4] <- rasters$evaluations$evaluations["MGCV", 2]
    output[(s-1)*20+11, (runs-1)*2+3] <- rasters$evaluations$evaluations["MGCVFIX", 1]
    output[(s-1)*20+11, (runs-1)*2+4] <- rasters$evaluations$evaluations["MGCVFIX", 2]
    output[(s-1)*20+12, (runs-1)*2+3] <- rasters$evaluations$evaluations["EARTH", 1]
    output[(s-1)*20+12, (runs-1)*2+4] <- rasters$evaluations$evaluations["EARTH", 2]
    output[(s-1)*20+13, (runs-1)*2+3] <- rasters$evaluations$evaluations["RPART", 1]
    output[(s-1)*20+13, (runs-1)*2+4] <- rasters$evaluations$evaluations["RPART", 2]
    output[(s-1)*20+14, (runs-1)*2+3] <- rasters$evaluations$evaluations["NNET", 1]
    output[(s-1)*20+14, (runs-1)*2+4] <- rasters$evaluations$evaluations["NNET", 2]
    output[(s-1)*20+15, (runs-1)*2+3] <- rasters$evaluations$evaluations["FDA", 1]
    output[(s-1)*20+15, (runs-1)*2+4] <- rasters$evaluations$evaluations["FDA", 2]
    output[(s-1)*20+16, (runs-1)*2+3] <- rasters$evaluations$evaluations["SVM", 1]
    output[(s-1)*20+16, (runs-1)*2+4] <- rasters$evaluations$evaluations["SVM", 2]
    output[(s-1)*20+17, (runs-1)*2+3] <- rasters$evaluations$evaluations["SVME", 1]
    output[(s-1)*20+17, (runs-1)*2+4] <- rasters$evaluations$evaluations["SVME", 2]
    output[(s-1)*20+18, (runs-1)*2+3] <- rasters$evaluations$evaluations["BIOCLIM", 1]
    output[(s-1)*20+18, (runs-1)*2+4] <- rasters$evaluations$evaluations["BIOCLIM", 2]
    output[(s-1)*20+19, (runs-1)*2+3] <- rasters$evaluations$evaluations["DOMAIN", 1]
    output[(s-1)*20+19, (runs-1)*2+4] <- rasters$evaluations$evaluations["DOMAIN", 2]
    output[(s-1)*20+20, (runs-1)*2+3] <- rasters$evaluations$evaluations["MAHAL", 1]
    output[(s-1)*20+20, (runs-1)*2+4] <- rasters$evaluations$evaluations["MAHAL", 2]

#3. if several rasters, continue to predict the other rasters

    if (length(xn) > 1) {
        for (n in 2:length(xn)) {
            xn.f <- xn[[n]]
            if(length(xn.f@title) == 0) {xn.f@title <- paste("stack", n, sep="")}
#            xn.f <- eval(as.name(xn.focal))
            cat(paste("\n", "Predictions for species: ", RASTER.species.name1, " for rasterStack: ", xn.f@title, sep = ""))
            cat(paste("\n", "NOTE: models are not re-calibrated (except if they failed earlier)", "\n\n", sep = ""))
            rasters2 <- ensemble.raster(
                x=x, p=ps, a=as, an=1000, ext=ext, k=k.raster, pt=NULL, at=NULL, xn=xn.f, models.keep=FALSE, 
                RASTER.species.name=RASTER.species.name1, RASTER.stack.name=xn.f@title, RASTER.format=RASTER.format, RASTER.datatype=RASTER.datatype, RASTER.NAflag=RASTER.NAflag,
                RASTER.models.overwrite=TRUE,
                threshold.method=threshold.method,
                ENSEMBLE.decay=1, ENSEMBLE.best=0, ENSEMBLE.min=0, 
                input.weights=splits$output.weights, models.list=rasters$models,
                MAXENT=0, GBM=0, GBMSTEP=0, RF=0, GLM=0, GLMSTEP=0, GAM=0, GAMSTEP=0, MGCV=0, MGCVFIX=0,
                EARTH=0, RPART=0, NNET=0, FDA=0, SVM=0, SVME=0, BIOCLIM=0, DOMAIN=0, MAHAL=0,  
                Yweights=Yweights, factors=factors, dummy.vars=dummy.vars,
                evaluation.strip=FALSE, 
                MAXENT.path=MAXENT.path,
                GBM.formula=GBM.formula, GBM.n.trees=GBM.n.trees, 
                GBMSTEP.gbm.x=GBMSTEP.gbm.x, GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
                GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction, 
                RF.formula=RF.formula, RF.ntree=RF.ntree, RF.mtry=RF.mtry, 
                GLM.formula=GLM.formula, GLM.family=GLM.family,
                GLMSTEP.k=GLMSTEP.k, GLMSTEP.steps=GLMSTEP.steps, STEP.formula=STEP.formula, GLMSTEP.scope=GLMSTEP.scope,
                GAM.formula=GAM.formula, GAM.family=GAM.family,
                GAMSTEP.steps=GAMSTEP.steps, GAMSTEP.scope=GAMSTEP.scope, GAMSTEP.OLD=rasters$GAMSTEP, GAMSTEP.pos=GAMSTEP.pos,
                MGCV.formula=MGCV.formula, MGCV.select=MGCV.select,
                MGCVFIX.formula=MGCVFIX.formula, 
                EARTH.formula=EARTH.formula, EARTH.glm=EARTH.glm, 
                RPART.formula=RPART.formula, RPART.xval=RPART.xval, 
                NNET.formula=NNET.formula, NNET.size=NNET.size, NNET.decay=NNET.decay, 
                FDA.formula=FDA.formula,  
                SVM.formula=SVM.formula, 
                SVME.formula=SVME.formula, 
                MAHAL.shape=MAHAL.shape)  
            }
        }

# end for all the species and all runs

    }
    }
    }

    cat(paste("\n", "summary of results", "\n\n", sep = ""))
    print(output)
    cat(paste("\n\n"))
    return(list(table=output, call=match.call() ))

}



