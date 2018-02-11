`ensemble.batch` <- function(
    x=NULL, xn=c(x), 
    species.presence=NULL, species.absence=NULL, 
    presence.min=20, thin.km=0.1,
    an=1000, excludep=FALSE, 
    get.block=FALSE,
    SSB.reduce=FALSE, CIRCLES.d=250000,
    k.splits=4, k.test=0, 
    n.ensembles=1, 
    VIF.max=10, VIF.keep=NULL,
    SINK=FALSE, CATCH.OFF=FALSE,
    RASTER.format="raster", RASTER.datatype="INT2S", RASTER.NAflag=-32767,
    KML.out=FALSE, KML.maxpixels=100000, KML.blur=10,
    models.save=FALSE,
    threshold.method="spec_sens", threshold.sensitivity=0.9, threshold.PresenceAbsence=FALSE,
    ENSEMBLE.best=0, ENSEMBLE.min=0.7, ENSEMBLE.exponent=1, ENSEMBLE.weight.min=0.05,
    input.weights=NULL,
    MAXENT=1, MAXLIKE=1, GBM=1, GBMSTEP=0, RF=1, GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=1, MGCV=1, 
    MGCVFIX=0, EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, GLMNET=1,
    BIOCLIM.O=0, BIOCLIM=1, DOMAIN=1, 
    MAHAL=1, MAHAL01=1, 
    PROBIT=FALSE,
    Yweights="BIOMOD", 
    layer.drops=NULL, factors=NULL, dummy.vars=NULL, 
    formulae.defaults=TRUE, maxit=100,
    MAXENT.a=NULL, MAXENT.an=10000, MAXENT.path=paste(getwd(), "/models/maxent", sep=""),
    MAXLIKE.formula=NULL, MAXLIKE.method="BFGS",
    GBM.formula=NULL, GBM.n.trees=2001,
    GBMSTEP.gbm.x=2:(1+raster::nlayers(x)), GBMSTEP.tree.complexity=5, GBMSTEP.learning.rate=0.005, 
    GBMSTEP.bag.fraction=0.5, GBMSTEP.step.size=100,
    RF.formula=NULL, RF.ntree=751, RF.mtry=floor(sqrt(raster::nlayers(x))), 
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
    GLMNET.nlambda=100, GLMNET.class=FALSE,
    BIOCLIM.O.fraction=0.9,
    MAHAL.shape=1  
)
{
    .BiodiversityR <- new.env()
#
    k.test <- as.integer(k.test)
    k.splits <- as.integer(k.splits)
    if (k.splits < 1) {
        cat(paste("\n", "NOTE: parameter k.splits was set to be smaller than 1", sep = ""))
        cat(paste("\n", "default value of 4 therefore set for parameter k.splits", sep = ""))
        k.splits <- 4
    }
    n.ensembles <- as.integer(n.ensembles)
    if (n.ensembles < 1) {n.ensembles <- 1}
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if (is.null(xn) == T) {
        cat(paste("\n", "NOTE: new rasterStack assumed to be equal to the base rasterStack", sep = ""))
        xn <- x
    }
    xn <- c(xn)
# need to recalculate threshold for mean of ensembles
# therefore put x as first of new stacks
    if (n.ensembles > 1) {
        xn <- c(x, xn)
        i <- 1
        while (i < length(xn)) {
            i <- i+1
            if(identical(x, xn[[i]])) {xn[[i]] <- NULL}
        }
    }
    species.presence <- data.frame(species.presence)
    if (ncol(species.presence) < 2) {stop("species.presence expected to be 3-column data.frame with species, x (e.g., lon) and y (e.g., lat) columns")}
    if (ncol(species.presence) == 2) {
        cat(paste("\n", "species.presence was expected to be 3-column data.frame with columns representing species, x (e.g., lon) and y (e.g., lat)", sep = ""))        
        cat(paste("\n", "only two columns were provided, it is therefore assumed that these reflect x and y coordinates for a single species", "\n\n", sep = ""))
        species.name <- rep("Species001", nrow(species.presence))
        species.presence <- cbind(species.name, species.presence)
        species.presence <- data.frame(species.presence)
        species.presence[,2] <- as.numeric(species.presence[,2])
        species.presence[,3] <- as.numeric(species.presence[,3])
        names(species.presence) <- c("species", "x", "y")
    }
    if (ncol(species.presence) > 3) {
        cat(paste("\n", "species.presence was expected to be 3-column data.frame with species, x (e.g., lon) and y (e.g., lat) columns", sep = ""))        
        cat(paste("\n", "only first three columns used", "\n\n", sep = ""))
        species.presence <- species.presence[,c(1:3)]
        species.presence[,2] <- as.numeric(species.presence[,2])
        species.presence[,3] <- as.numeric(species.presence[,3])
        names(species.presence) <- c("species", "x", "y")
    }
    if (is.null(species.absence)==F) {species.absence <- data.frame(species.absence)}
    if (is.null(species.absence)==F && ncol(species.absence) < 2) {stop("species.absence expected to be a 2-column data.frame with x (e.g., lon) and y (e.g., lat),  or 3-column data.frame with species, x (e.g., lon) and y (e.g., lat) columns")}
    if (is.null(species.absence)==F && ncol(species.absence)> 3) {
        cat(paste("\n", "species.absence was expected to be 3-column data.frame with species, x (e.g., lon) and y (e.g., lat) columns", sep = ""))        
        cat(paste("\n", "only first three columns used", "\n\n", sep = ""))
        species.absence <- species.absence[,c(1:3)]
        species.absence[,2] <- as.numeric(species.absence[,2])
        species.absence[,3] <- as.numeric(species.absence[,3])
        names(species.absence) <- c("species", "x", "y")
    }
    if (is.null(species.absence)==F && ncol(species.absence) == 2) {
        cat(paste("\n", "species.absence was expected to be 3-column data.frame with species, x (e.g., lon) and y (e.g., lat) columns", sep = ""))        
        cat(paste("\n", "only two columns were provided, it is therefore assumed that these reflect x and y coordinates for absence locations to be used for each species run", "\n\n", sep = ""))
        species.absence[,1] <- as.numeric(species.absence[,1])
        species.absence[,2] <- as.numeric(species.absence[,2])
        names(species.absence) <- c("x", "y")
        as <- species.absence
    }

#
# check for variables below the maximum VIF
# note that the ensemble.VIF function does not remove factor variables
#
    VIF.result <- ensemble.VIF(x=x, VIF.max=VIF.max, keep=VIF.keep,
        layer.drops=layer.drops, factors=factors, dummy.vars=dummy.vars)

    layer.drops <- VIF.result$var.drops
    factors <- VIF.result$factors
    dummy.vars <- VIF.result$dummy.vars
# 
# process species by species
    species.names <- levels(droplevels(factor(species.presence[,1])))

    AUC.table.out <- AUC.ensemble.out <- output.weights.out <- ensemble.highest <- NULL

    for (s in 1:length(species.names)) {
        focal.species <- species.names[s]

        ps <- species.presence[species.presence[,1]==focal.species, c(2:3)]
        n.pres <- nrow(ps)

# check after spatial thinning if species has required minimum number of presence points
# already calculate all the spatially thinned data sets for each run

	if (thin.km > 0) {

        cat(paste("\n", "Generation of spatially thinned presence data sets for each ensemble", "\n\n", sep = ""))    

            ps.thins <- vector("list", n.ensembles)
            for (i in 1:n.ensembles) {
                ps.thins[[i]] <- ensemble.spatialThin(ps, thin.km=thin.km)
                if (nrow(ps.thins[[i]]) < n.pres) {n.pres <- nrow(ps.thins[[i]])}
            }
        }

        if (n.pres < presence.min) {
            if (thin.km > 0) {
                cat(paste("\n", "Species: ", focal.species, " only has ", n.pres, " presence locations in one of the spatially thinned data sets", sep = ""))
            }else{
                cat(paste("\n", "Species: ", focal.species, " only has ", n.pres, " presence locations", sep = ""))
            }
            cat(paste("\n", "This species therefore not included in batch processing", "\n\n", sep = ""))

        }else{

# create output file
    if (s==1) {dir.create("outputs", showWarnings = F)}
    paste.file <- paste(getwd(), "/outputs/", focal.species, "_output.txt", sep="")
    OLD.SINK <- TRUE
    if (sink.number(type="output") == 0) {OLD.SINK <- F}
    if (SINK==T && OLD.SINK==F) {
        if (file.exists(paste.file) == F) {
            cat(paste("\n", "NOTE: results captured in file: ", paste.file, "\n", sep = ""))
        }else{
            cat(paste("\n", "NOTE: results appended in file: ", paste.file, "\n", sep = ""))
        }
        cat(paste("\n\n", "RESULTS (ensemble.batch function)", "\n\n", sep=""), file=paste.file, append=T)
        sink(file=paste.file, append=T)
        cat(paste(date(), "\n", sep=""))
        print(match.call())
    }
#
    cat(paste("\n", "Evaluations for species: ", focal.species, "\n", sep = ""))
    ps <- species.presence[species.presence[,1]==focal.species, c(2:3)]

# repeat the whole process for n.ensembles

    RASTER.species.name1 <- focal.species
    for (runs in 1:n.ensembles) {
        if (n.ensembles > 1) { 
            cat(paste("\n", focal.species, ": ENSEMBLE ", runs, "\n\n", sep = ""))
            RASTER.species.name1 <- paste(focal.species, "_ENSEMBLE_", runs, sep="")
        }

    if (thin.km > 0) {
        ps <- ps.thins[[runs]]
    }

    if (is.null(species.absence)==F && ncol(species.absence) == 3) {
        as <- species.absence[species.absence[,1]==focal.species, c(2:3)]
    }

# random selection of background locations for each run
    if (is.null(species.absence)==T) {
        if (excludep == T) {
            as <- dismo::randomPoints(x[[1]], n=an, p=ps, excludep=T)
        }else{
            as <- dismo::randomPoints(x[[1]], n=an, p=NULL, excludep=F)
        }
    }

    assign("ps", ps, envir=.BiodiversityR)
    assign("as", as, envir=.BiodiversityR)

#1. first ensemble tests
    calibration1 <- ensemble.calibrate.weights(x=x, p=ps, a=as, k=k.splits, get.block=get.block,
        SSB.reduce=SSB.reduce, CIRCLES.d=CIRCLES.d,
        CATCH.OFF=CATCH.OFF,
        ENSEMBLE.tune=T,
        ENSEMBLE.best=ENSEMBLE.best, ENSEMBLE.min=ENSEMBLE.min, 
        ENSEMBLE.exponent=ENSEMBLE.exponent, ENSEMBLE.weight.min=ENSEMBLE.weight.min,
        species.name = RASTER.species.name1,
        threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence,
        input.weights=input.weights,
        MAXENT=MAXENT, MAXLIKE=MAXLIKE, GBM=GBM, GBMSTEP=GBMSTEP, RF=RF, GLM=GLM, GLMSTEP=GLMSTEP, 
        GAM=GAM, GAMSTEP=GAMSTEP, MGCV=MGCV, MGCVFIX=MGCVFIX, EARTH=EARTH, RPART=RPART, 
        NNET=NNET, FDA=FDA, SVM=SVM, SVME=SVME, GLMNET=GLMNET,
        BIOCLIM.O=BIOCLIM.O, BIOCLIM=BIOCLIM, DOMAIN=DOMAIN, MAHAL=MAHAL, MAHAL01=MAHAL01,
        PROBIT=PROBIT, VIF=T,
        Yweights=Yweights, 
        layer.drops=layer.drops, factors=factors, dummy.vars=dummy.vars,
        maxit=maxit,
        MAXENT.a=MAXENT.a, MAXENT.an=MAXENT.an, 
        MAXENT.path=MAXENT.path,
        MAXLIKE.formula=MAXLIKE.formula, MAXLIKE.method=MAXLIKE.method,
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
        GLMNET.nlambda=GLMNET.nlambda, GLMNET.class=GLMNET.class,
        BIOCLIM.O.fraction=BIOCLIM.O.fraction,
        MAHAL.shape=MAHAL.shape)

    x.batch <- calibration1$x
    p.batch <- calibration1$p
    a.batch <- calibration1$a
    MAXENT.a.batch <- calibration1$MAXENT.a
    var.names.batch <- calibration1$var.names
    factors.batch <- calibration1$factors
    dummy.vars.batch <- calibration1$dummy.vars
    dummy.vars.noDOMAIN.batch <- calibration1$dummy.vars.noDOMAIN

    AUC.table <- calibration1$AUC.table
    rownames(AUC.table) <- paste(rownames(AUC.table), "_", runs, sep="")
    AUC.table <- cbind(rep(runs, nrow(AUC.table)), AUC.table, rep(NA, nrow(AUC.table)))
    colnames(AUC.table)[1] <- "ensemble"
    colnames(AUC.table)[ncol(AUC.table)] <- "final.calibration"

    if (runs == 1) {
        AUC.table.out <- AUC.table
    }else{
        AUC.table.out <- rbind(AUC.table.out, AUC.table)
    }

    AUC.ensemble <- calibration1$AUC.with.suggested.weights
    AUC.ensemble <- c(runs, AUC.ensemble)
    names(AUC.ensemble)[1] <- "ensemble"

    if (runs == 1) {
        AUC.ensemble.out <- AUC.ensemble
    }else{
        AUC.ensemble.out <- rbind(AUC.ensemble.out, AUC.ensemble)
    }    

#2. calibrate final model

#    xn.f <- eval(as.name(xn.focal))
    cat(paste("\n", "Final model calibrations for species: ", RASTER.species.name1,  "\n", sep = ""))

    cat(paste("\n\n", "Input weights for ensemble.calibrate.models are average weights determined by ensemble.calibrate.weights function", "\n", sep=""))        
    output.weights <- calibration1$output.weights
    print(output.weights)

    output.weights <- c(runs, output.weights)
    names(output.weights)[1] <- "ensemble"

    if (runs == 1) {
        output.weights.out <- output.weights
    }else{
        output.weights.out <- rbind(output.weights.out, output.weights)
        rownames(output.weights.out) <- c(1:nrow(output.weights.out))
    }  

    if (sum(output.weights) > 0) {

    calibration2 <- ensemble.calibrate.models(
        x=x.batch, p=p.batch, a=a.batch, k=k.test, pt=NULL, at=NULL,
        models.keep=TRUE, evaluations.keep=TRUE,
        PLOTS=FALSE, CATCH.OFF=CATCH.OFF, 
        models.save=models.save, species.name=RASTER.species.name1,
        ENSEMBLE.tune=F,
        input.weights=output.weights,
        threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence,
        PROBIT=PROBIT, VIF=T,
        Yweights=Yweights, 
        factors=factors.batch, dummy.vars=dummy.vars.batch,
        maxit=maxit,
        MAXENT.a=MAXENT.a.batch, MAXENT.path=MAXENT.path,
        MAXLIKE.formula=MAXLIKE.formula, MAXLIKE.method=MAXLIKE.method,
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
        GLMNET.nlambda=GLMNET.nlambda, GLMNET.class=GLMNET.class,
        BIOCLIM.O.fraction=BIOCLIM.O.fraction,
        MAHAL.shape=MAHAL.shape) 

        AUC.table.out[which(rownames(AUC.table.out) == paste("MAXENT", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["MAXENT"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("MAXLIKE", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["MAXLIKE"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("GBM", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["GBM"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("GBMSTEP", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["GBMSTEP"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("RF", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["RF"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("GLM", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["GLM"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("GLMSTEP", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["GLMSTEP"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("GAM", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["GAM"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("GAMSTEP", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["GAMSTEP"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("MGCV", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["MGCV"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("MGCVFIX", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["MGCVFIX"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("EARTH", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["EARTH"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("RPART", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["RPART"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("NNET", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["NNET"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("FDA", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["FDA"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("SVM", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["SVM"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("SVME", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["SVME"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("GLMNET", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["GLMNET"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("BIOCLIM.O", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["BIOCLIM.O"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("BIOCLIM", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["BIOCLIM"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("DOMAIN", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["DOMAIN"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("MAHAL", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["MAHAL"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("MAHAL01", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["MAHAL01"]
        AUC.table.out[which(rownames(AUC.table.out) == paste("ENSEMBLE", "_", runs, sep="")), ncol(AUC.table.out)] <- calibration2$AUC.calibration["ENSEMBLE"]

#3. predict for all the rasters

    for (n in 1:length(xn)) {
        xn.f <- raster::stack(xn[[n]])
        xn.f <- raster::subset(xn.f, subset=var.names.batch)
        xn.f <- raster::stack(xn.f)
        if(length(xn.f@title) == 0) {xn.f@title <- paste("stack_", n, sep="")}
        if (gsub(".", "_", xn.f@title, fixed=T) != xn.f@title) {cat(paste("\n", "WARNING: title of stack (", xn.f@title, ") contains '.'", "\n\n", sep = ""))}
        cat(paste("\n", "Predictions for species: ", RASTER.species.name1, " for rasterStack: ", xn.f@title, "\n\n", sep = ""))

        if (n == 1) {
            rasters2 <- ensemble.raster(xn=xn.f, 
                models.list=calibration2$models,            
                RASTER.species.name=RASTER.species.name1, 
                evaluate=T, p=p.batch, a=a.batch,
                RASTER.format=RASTER.format, RASTER.datatype=RASTER.datatype, RASTER.NAflag=RASTER.NAflag,
                KML.out=KML.out, KML.maxpixels=KML.maxpixels, KML.blur=KML.blur)
        }else{
            rasters2 <- ensemble.raster(xn=xn.f, 
                models.list=calibration2$models,            
                RASTER.species.name=RASTER.species.name1, 
                RASTER.format=RASTER.format, RASTER.datatype=RASTER.datatype, RASTER.NAflag=RASTER.NAflag,
                KML.out=KML.out, KML.maxpixels=KML.maxpixels, KML.blur=KML.blur)
        }

        if(runs==n.ensembles && n.ensembles>1 && RASTER.format=="raster") {

# recalculate threshold for mean of predictions with calibration stack (xn[[1]])
# use threshold to calculate mean ensemble, ensemble count, ensemble presence and ensemble sd
            if (n == 1) {
                calibrate.mean <- NULL
                calibrate.mean <- ensemble.mean(RASTER.species.name=focal.species, RASTER.stack.name=xn.f@title,
                    positive.filters = c("grd", "_ENSEMBLE_"), negative.filters = c("xml"), 
                    RASTER.format=RASTER.format, RASTER.datatype=RASTER.datatype, RASTER.NAflag=RASTER.NAflag,
                    KML.out=KML.out, KML.maxpixels=KML.maxpixels, KML.blur=KML.blur,
                    p=p.batch, a=a.batch,
                    pt = NULL, at = NULL,
                    threshold = -1,
                    threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence)
                cat(paste("\n", "threshold for mean suitability: ", calibrate.mean$threshold, "\n", sep = ""))
            }else{
                ensemble.mean(RASTER.species.name=focal.species, RASTER.stack.name=xn.f@title,
                    positive.filters = c("grd", "_ENSEMBLE_"), negative.filters = c("xml"), 
                    RASTER.format=RASTER.format, RASTER.datatype=RASTER.datatype, RASTER.NAflag=RASTER.NAflag,
                    KML.out=KML.out, KML.maxpixels=KML.maxpixels, KML.blur=KML.blur,
                    p=NULL, a=NULL,
                    pt = NULL, at = NULL,
                    threshold = calibrate.mean$threshold,
                    threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence)
            }
        }
    }


# sum output weights > 0 loop
    }

# n ensembles loop
    }

    cat(paste("\n\n", "All AUC results for species: ", RASTER.species.name1, " for rasterStack: ", xn.f@title, "\n\n", sep=""))
    print(AUC.table.out)

    if (n.ensembles > 1) {
        ensemble.highest <- AUC.ensemble.out[which.max(AUC.ensemble.out[, "MEAN.T"]), "ensemble"]
        cat(paste("\n", "ensemble with highest average AUC is ensemble: ", ensemble.highest, "\n", sep = ""))
    }else{
        ensemble.highest <- 1
    }

# if (sufficient presence locations) loop
    if (SINK==T && OLD.SINK==F) {sink(file=NULL, append=T)}
    }

# s (species) loop
    }

    result <- list(species=species.names, AUC.table=AUC.table.out, AUC.ensemble.selected.weights=AUC.ensemble.out, output.weights=output.weights.out, ensemble.highest.AUC=ensemble.highest, call=match.call())

    cat(paste("\n\n", "(all calibrations and projections finalized by function ensemble.batch)", "\n\n", sep=""))

    return(result)

}
