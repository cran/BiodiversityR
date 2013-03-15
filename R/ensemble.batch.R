`ensemble.batch` <- function(
    x, species.presence, species.absence=NULL, an=1000, ext=NULL, k=0, k.splits=5, xn=c(x), 
    layer.drops=NULL,
    digits=2,
    RASTER.format="raster", RASTER.datatype="INT2S", RASTER.NAflag=-32767,
    threshold.method="spec_sens",
    ENSEMBLE.decay=1, ENSEMBLE.multiply=TRUE, ENSEMBLE.best=0,
    input.weights=NULL,
    MAXENT=1, GBM=1, GBMSTEP=1, RF=1, GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=1, MGCV=1, MGCVFIX=0,
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, BIOCLIM=1, DOMAIN=1, MAHAL=1, 
    GEODIST=0,
    Yweights="BIOMOD", factors=NULL, dummy.vars=NULL, 
    formulae.defaults=TRUE, maxit=100,
    GBM.formula=NULL, GBM.n.trees=2000,
    GBMSTEP.gbm.x=2:(1+nlayers(x)), GBMSTEP.tree.complexity=5, GBMSTEP.learning.rate=0.005, 
    GBMSTEP.bag.fraction=0.5, GBMSTEP.step.size=100,
    RF.formula=NULL, RF.ntree=750, RF.mtry=log(nlayers(x)), 
    GLM.formula=NULL, GLM.family=binomial(link="logit"), 
    GLMSTEP.steps=1000, STEP.formula=NULL, GLMSTEP.scope=NULL, GLMSTEP.k=2,
    GAM.formula=NULL, GAM.family=binomial(link="logit"), 
    GAMSTEP.steps=1000, GAMSTEP.scope=NULL,
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
    if (k < 1) {stop("Parameter k smaller than 1")}
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
    species.presence <- data.frame(species.presence)
    if (ncol(species.presence) < 3) {stop("species.presence expected to be 3-column data.frame with species, x (lon) and y (lat) columns")}
    if (ncol(species.presence) > 3) {
        cat(paste("\n", "species.presence was expected to be 3-column data.frame with species, x (lon) and y (lat) columns", "\n\n", sep = ""))        
        cat(paste("\n", "only first three columns used", "\n\n", sep = ""))
        species.presence <- species.presence[,c(1:3)] 
    }
    if (is.null(species.absence)==F && ncol(species.absence) < 2) {stop("species.absence expected to be a 2-column data.frame with x (lon) and y (lat),  or 3-column data.frame with species, x (lon) and y (lat) columns")}
    if (is.null(species.absence)==F && ncol(species.absence)> 3) {
        cat(paste("\n", "species.absence was expected to be 3-column data.frame with species, x (lon) and y (lat) columns", "\n\n", sep = ""))        
        cat(paste("\n", "only first three columns used", "\n\n", sep = ""))
        species.absence <- species.absence[,c(1:3)] 
    }
    if (is.null(species.absence)==F && ncol(species.absence) == 2) {as <- data.frame(species.absence)}
    if (is.null(species.absence)==T) {as <- randomPoints(x, n=an, ext=ext)}
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
# process species by species
    species.names <- levels(as.factor(species.presence[,1]))
    output <- NULL
    for (s in 1:length(species.names)) {
        focal.species <- species.names[s]
        cat(paste("\n", "Evaluations for species: ", focal.species, "\n\n", sep = ""))
        ps <- species.presence[species.presence[,1]==focal.species, c(2:3)]
        if (is.null(species.absence)==F && ncol(species.absence) == 3) {
            as <- species.absence[species.absence[,1]==focal.species, c(2:3)]
        }

#1. first ensemble tests
    splits <- ensemble.test.splits(x=x, p=ps, a=as, ext=ext, k=k.splits, 
        VIF=T,
        digits=digits, PLOTS=F, 
        input.weights=input.weights,
        MAXENT=MAXENT, GBM=GBM, GBMSTEP=GBMSTEP, RF=RF, GLM=GLM, GLMSTEP=GLMSTEP, 
        GAM=GAM, GAMSTEP=GAMSTEP, MGCV=MGCV, EARTH=EARTH, RPART=RPART, 
        NNET=NNET, FDA=FDA, SVM=SVM, SVME=SVME, BIOCLIM=BIOCLIM, DOMAIN=DOMAIN, MAHAL=MAHAL,
        GEODIST=GEODIST, GEODIST.file.name=focal.species, RASTER.format=RASTER.format,  
        Yweights=Yweights, factors=factors, dummy.vars=dummy.vars,
        maxit=maxit,
        GBM.formula=GBM.formula, GBM.n.trees=GBM.n.trees,
        GBMSTEP.gbm.x=GBMSTEP.gbm.x, GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
        GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction,
        GBMSTEP.step.size=GBMSTEP.step.size,
        RF.formula=RF.formula, RF.ntree=RF.ntree, RF.mtry=RF.mtry, 
        GLM.formula=GLM.formula, GLM.family=GLM.family, 
        GLMSTEP.k=GLMSTEP.k, GLMSTEP.steps=GLMSTEP.steps, STEP.formula=STEP.formula, GLMSTEP.scope=GLMSTEP.scope, 
        GAM.formula=GAM.formula, GAM.family=GAM.family, 
        GAMSTEP.steps=GAMSTEP.steps, GAMSTEP.scope=GAMSTEP.scope,
        MGCV.formula=MGCV.formula, MGCV.select=MGCV.select,
        MGCVFIX.formula=MGCVFIX.formula, 
        EARTH.formula=EARTH.formula, EARTH.glm=EARTH.glm,
        RPART.formula=RPART.formula, RPART.xval=RPART.xval, 
        NNET.formula=NNET.formula, NNET.size=NNET.size, NNET.decay=NNET.decay,
        FDA.formula=FDA.formula, SVM.formula=SVM.formula, SVME.formula=SVME.formula)

#2. calibrate models and predict first raster
#    xn <- as.character(xn)    
    xn.f <- xn[[1]]
    if(length(xn.f@title) == 0) {xn.f@title <- "stack1"}
#    xn.f <- eval(as.name(xn.focal))
    cat(paste("\n", "Predictions for species: ", focal.species, " for rasterStack: ", xn.f@title,  "\n\n", sep = ""))
    rasters <- ensemble.raster(
        x=x, p=ps, a=as, an=1000, ext=ext, k=k, pt=NULL, at=NULL, xn=xn.f, models.keep=TRUE, 
        RASTER.species.name=focal.species, RASTER.stack.name=xn.f@title, RASTER.format=RASTER.format, RASTER.datatype=RASTER.datatype, RASTER.NAflag=RASTER.NAflag,
        RASTER.models.overwrite=TRUE,
        threshold.method=threshold.method,
        ENSEMBLE.decay=ENSEMBLE.decay, ENSEMBLE.multiply=ENSEMBLE.multiply, ENSEMBLE.best=ENSEMBLE.best,
        input.weights=splits[,"MEAN"],
        MAXENT=0, GBM=0, GBMSTEP=0, RF=0, GLM=0, GLMSTEP=0, GAM=0, GAMSTEP=0, MGCV=0, MGCVFIX=0,
        EARTH=0, RPART=0, NNET=0, FDA=0, SVM=0, SVME=0, BIOCLIM=0, DOMAIN=0, MAHAL=0,  
        Yweights=Yweights, factors=factors, dummy.vars=dummy.vars,
        evaluation.strip=FALSE, 
        maxit=maxit,
        MAXENT.OLD=NULL,
        GBM.formula=GBM.formula, GBM.n.trees=GBM.n.trees, GBM.OLD=NULL,    
        GBMSTEP.gbm.x=GBMSTEP.gbm.x, GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
        GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction, GBMSTEP.OLD=NULL,
        RF.formula=RF.formula, RF.ntree=RF.ntree, RF.mtry=RF.mtry, RF.OLD=NULL,
        GLM.formula=GLM.formula, GLM.family=GLM.family, GLM.OLD=NULL,
        GLMSTEP.k=GLMSTEP.k, GLMSTEP.steps=GLMSTEP.steps, STEP.formula=STEP.formula, GLMSTEP.scope=GLMSTEP.scope, GLMSTEP.OLD=NULL,
        GAM.formula=GAM.formula, GAM.family=GAM.family, GAM.OLD=NULL, 
        GAMSTEP.steps=GAMSTEP.steps, GAMSTEP.scope=GAMSTEP.scope, GAMSTEP.OLD=NULL,
        MGCV.formula=MGCV.formula, MGCV.select=MGCV.select, MGCV.OLD=NULL,
        MGCVFIX.formula=MGCVFIX.formula, MGCVFIX.OLD=NULL,
        EARTH.formula=EARTH.formula, EARTH.glm=EARTH.glm, EARTH.OLD=NULL,
        RPART.formula=RPART.formula, RPART.xval=RPART.xval, RPART.OLD=NULL,
        NNET.formula=NNET.formula, NNET.size=NNET.size, NNET.decay=NNET.decay, NNET.OLD=NULL,
        FDA.formula=FDA.formula,  FDA.OLD=NULL, 
        SVM.formula=SVM.formula, SVM.OLD=NULL, 
        SVME.formula=SVME.formula, SVME.OLD=NULL,
        BIOCLIM.OLD=NULL, DOMAIN.OLD=NULL, MAHAL.OLD=NULL)   

#3. if several rasters, continue to predict the other rasters
# (In this case, x is also set to be the new layer to generate evaluations for the new layer)

    if (length(xn) > 1) {
        for (n in 2:length(xn)) {
            xn.f <- xn[[n]]
            if(length(xn.f@title) == 0) {xn.f@title <- paste("stack", n, sep="")}
#            xn.f <- eval(as.name(xn.focal))
            cat(paste("\n", "Predictions for species: ", focal.species, " for rasterStack: ", xn.f@title,  "\n\n", sep = ""))
            rasters2 <- ensemble.raster(
                x=xn.f, p=ps, a=as, an=1000, ext=ext, k=k, pt=NULL, at=NULL, xn=xn.f, models.keep=FALSE, 
                RASTER.species.name=focal.species, RASTER.stack.name=xn.f@title, RASTER.format=RASTER.format, RASTER.datatype=RASTER.datatype, RASTER.NAflag=RASTER.NAflag,
                RASTER.models.overwrite=TRUE,
                threshold.method=threshold.method,
                ENSEMBLE.decay=ENSEMBLE.decay, ENSEMBLE.multiply=ENSEMBLE.multiply, ENSEMBLE.best=ENSEMBLE.best,
                input.weights=splits[,"MEAN"],
                MAXENT=0, GBM=0, GBMSTEP=0, RF=0, GLM=0, GLMSTEP=0, GAM=0, GAMSTEP=0, MGCV=0, MGCVFIX=0,
                EARTH=0, RPART=0, NNET=0, FDA=0, SVM=0, SVME=0, BIOCLIM=0, DOMAIN=0, MAHAL=0,  
                Yweights=Yweights, factors=factors, dummy.vars=dummy.vars,
                evaluation.strip=FALSE, 
                MAXENT.OLD=rasters$MAXENT,
                GBM.formula=GBM.formula, GBM.n.trees=GBM.n.trees, GBM.OLD=rasters$GBM,    
                GBMSTEP.gbm.x=GBMSTEP.gbm.x, GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
                GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction, GBMSTEP.OLD=rasters$GBMSTEP,
                RF.formula=RF.formula, RF.ntree=RF.ntree, RF.mtry=RF.mtry, RF.OLD=rasters$RF,
                GLM.formula=GLM.formula, GLM.family=GLM.family, GLM.OLD=rasters$GLM,
                GLMSTEP.k=GLMSTEP.k, GLMSTEP.steps=GLMSTEP.steps, STEP.formula=STEP.formula, GLMSTEP.scope=GLMSTEP.scope, GLMSTEP.OLD=rasters$GLMSTEP,
                GAM.formula=GAM.formula, GAM.family=GAM.family, GAM.OLD=rasters$GAM, 
                GAMSTEP.steps=GAMSTEP.steps, GAMSTEP.scope=GAMSTEP.scope, GAMSTEP.OLD=rasters$GAMSTEP,
                MGCV.formula=MGCV.formula, MGCV.select=MGCV.select, MGCV.OLD=rasters$MGCV,
                MGCVFIX.formula=MGCVFIX.formula, MGCVFIX.OLD=rasters$MGCVFIX,
                EARTH.formula=EARTH.formula, EARTH.glm=EARTH.glm, EARTH.OLD=rasters$EARTH,
                RPART.formula=RPART.formula, RPART.xval=RPART.xval, RPART.OLD=rasters$RPART,
                NNET.formula=NNET.formula, NNET.size=NNET.size, NNET.decay=NNET.decay, NNET.OLD=rasters$NNET,
                FDA.formula=FDA.formula,  FDA.OLD=rasters$FDA, 
                SVM.formula=SVM.formula, SVM.OLD=rasters$SVM, 
                SVME.formula=SVME.formula, SVME.OLD=rasters$SVME,
                BIOCLIM.OLD=rasters$BIOCLIM, DOMAIN.OLD=rasters$DOMAIN, MAHAL.OLD=rasters$MAHAL)  
            }
        }

# end for all the species
    }

    cat(paste("\n", "end of batch processing", "\n\n", sep = ""))
}




