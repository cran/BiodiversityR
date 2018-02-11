`ensemble.calibrate.weights` <- function(
    x=NULL, p=NULL,
    a=NULL, an=1000, 
    get.block=FALSE, SSB.reduce=FALSE, CIRCLES.d=100000,
    excludep=FALSE, k=4, 
    VIF=FALSE, COR=FALSE,
    SINK=FALSE, PLOTS=FALSE, CATCH.OFF=FALSE,
    data.keep=FALSE,
    species.name = "Species001",
    threshold.method="spec_sens", threshold.sensitivity=0.9, threshold.PresenceAbsence=FALSE,
    ENSEMBLE.tune=FALSE, 
    ENSEMBLE.best=0, ENSEMBLE.min=0.7, ENSEMBLE.exponent=1, ENSEMBLE.weight.min=0.05,
    input.weights=NULL,
    MAXENT=1, MAXLIKE=1, GBM=1, GBMSTEP=1, RF=1, GLM=1, GLMSTEP=1, 
    GAM=1, GAMSTEP=1, MGCV=1, MGCVFIX=0,
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, GLMNET=1,
    BIOCLIM.O=0, BIOCLIM=1, DOMAIN=1, MAHAL=1, MAHAL01=1,
    PROBIT=FALSE,
    Yweights="BIOMOD", 
    layer.drops=NULL, factors=NULL, dummy.vars=NULL,
    formulae.defaults=TRUE, maxit=100,
    MAXENT.a=NULL, MAXENT.an=10000, 
    MAXENT.path=paste(getwd(), "/models/maxent_", species.name,  sep=""), 
    MAXLIKE.formula=NULL, MAXLIKE.method="BFGS",
    GBM.formula=NULL, GBM.n.trees=2001, 
    GBMSTEP.gbm.x=2:(length(var.names)+1), GBMSTEP.tree.complexity=5, GBMSTEP.learning.rate=0.005, 
    GBMSTEP.bag.fraction=0.5, GBMSTEP.step.size=100, 
    RF.formula=NULL, RF.ntree=751, RF.mtry=floor(sqrt(length(var.names))), 
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
    GLMNET.nlambda=100, GLMNET.class=FALSE,
    BIOCLIM.O.fraction=0.9,
    MAHAL.shape=1 
)
{
    .BiodiversityR <- new.env()

#    if (! require(dismo)) {stop("Please install the dismo package")}

    if(is.null(x) == T) {stop("value for parameter x is missing (RasterStack object)")}
    if(inherits(x,"RasterStack") == F) {stop("x is not a RasterStack object")}
    if(is.null(p) == T) {stop("presence locations are missing (parameter p)")}
    if(is.null(p) == F) {names(p) <- c("x", "y")}
    if(is.null(a) == F) {names(a) <- c("x", "y")}
    if(is.null(MAXENT.a) == F) {names(MAXENT.a) <- c("x", "y")}
#
    k <- as.integer(k)
    if (k < 2) {
        cat(paste("\n", "NOTE: parameter k was set to be smaller than 2", sep = ""))
        cat(paste("\n", "default value of 4 therefore set for parameter k", "\n", sep = ""))
        k <- 4
    }
#
    if (is.null(layer.drops) == F) {
        layer.drops <- as.character(layer.drops)
        if (is.null(x)==F) {x <- raster::dropLayer(x, which(names(x) %in% layer.drops))}
        x <- raster::stack(x)
        factors <- as.character(factors)
        dummy.vars <- as.character(dummy.vars)
        nd <- length(layer.drops)
        for (i in 1:nd) {
            if (is.null(factors) == F) {
                factors <- factors[factors != layer.drops[i]]
                if(length(factors) == 0) {factors <- NULL}
            }
            if (is.null(dummy.vars) == F) {
                dummy.vars <- dummy.vars[dummy.vars != layer.drops[i]]
                if(length(dummy.vars) == 0) {dummy.vars <- NULL}
            }
        }
        if(length(layer.drops) == 0) {layer.drops <- NULL}
    }
#
    output.rownames <- c("MAXENT", "MAXLIKE", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "GLMNET", 
        "BIOCLIM.O", "BIOCLIM", "DOMAIN", "MAHAL", "MAHAL01", "ENSEMBLE")
#
    if(length(ENSEMBLE.exponent) > 1 || length(ENSEMBLE.best) > 1 || length(ENSEMBLE.min) > 1) {ENSEMBLE.tune <- TRUE}
    if(ENSEMBLE.tune == F) {
        output <- array(0, dim=c(length(output.rownames), k+1))
        rownames(output) <- output.rownames
        colnames(output) <- c(paste("T_", c(1:k), sep=""),"MEAN")
    }else{
        output <- array(0, dim=c(length(output.rownames), 2*k+2))
        rownames(output) <- output.rownames
        colnames(output) <- c(paste("T_", c(1:k), sep=""), "MEAN.T", paste("S_", c(1:k), sep=""), "MEAN")
    }

# keep data for final output.weights checks with suggested weights
    TestData.all <- vector("list", k)

# create output file
    dir.create("outputs", showWarnings = F)
    paste.file <- paste(getwd(), "/outputs/", species.name, "_output.txt", sep="")
    OLD.SINK <- TRUE
    if (sink.number(type="output") == 0) {OLD.SINK <- F}
    if (SINK==T && OLD.SINK==F) {
        if (file.exists(paste.file) == F) {
            cat(paste("\n", "NOTE: results captured in file: ", paste.file, "\n", sep = ""))
        }else{
            cat(paste("\n", "NOTE: results appended in file: ", paste.file, "\n", sep = ""))
        }
        cat(paste("\n\n", "RESULTS (ensemble.calibrate.weights function)", "\n\n", sep=""), file=paste.file, append=T)
        sink(file=paste.file, append=T)
        cat(paste(date(), "\n", sep=""))
        print(match.call())
    }

# 
# run ensemble.calibrate.models first to obtain MAXENT.a and var.names
    tests <- ensemble.calibrate.models(x=x, 
        p=p, a=a, an=an, pt=NULL, at=NULL, excludep=excludep, k=0, 
        TrainData=NULL, 
        VIF=F, COR=F,
        PLOTS=PLOTS, evaluations.keep=T, models.keep=F,
        ENSEMBLE.tune=F,
        ENSEMBLE.exponent=1, ENSEMBLE.best=1, ENSEMBLE.min=0.7,
        MAXENT=0, MAXLIKE=0, GBM=0, GBMSTEP=0, RF=0, GLM=0, GLMSTEP=0, 
        GAM=0, GAMSTEP=0, MGCV=0, MGCVFIX=0, EARTH=0, RPART=0, 
        NNET=0, FDA=0, SVM=0, SVME=0, GLMNET=0,
        BIOCLIM.O=0, BIOCLIM=0, DOMAIN=0, MAHAL=0, MAHAL01=0,
        MAXENT.a=MAXENT.a, MAXENT.an=MAXENT.an, 
        factors=factors)

    var.names <- tests$evaluations$var.names
    MAXENT.a <- tests$evaluations$MAXENT.a

    factors2 <- NULL
    if (is.null(factors) == F) {
        factors2 <- factors[which(factors %in% var.names)]
        if (length(factors2) == 0) {factors2 <- NULL}
    }
    dummy.vars2 <- NULL
    if (is.null(dummy.vars) == F) {
        dummy.vars2 <- dummy.vars[which(dummy.vars %in% var.names)]
        if (length(dummy.vars2) == 0) {dummy.vars2 <- NULL}
    }

    p.all <- tests$evaluations$p
    a.all <- tests$evaluations$a

    if (get.block == F) {
        groupp <- dismo::kfold(p.all, k=k)
        groupa <- dismo::kfold(a.all, k=k)
    }else{
        blocks <- ENMeval::get.block(occ=p.all, bg.coords=a.all)
        groupp <- blocks$occ.grp
        groupa <- blocks$bg.grp
        k <- 4
    }

# Start cross-validations
    
    for (i in 1:k){
        cat(paste(species.name, " ", k, "-FOLD CROSS-VALIDATION RUN: ", i, "\n", sep = ""))

            p1 <- p.all[groupp != i,]
            p2 <- p.all[groupp == i,]
            a1 <- a.all[groupa != i,]
            a2 <- a.all[groupa == i,]

            tests <- ensemble.calibrate.models(x=x, 
                TrainData=NULL, TestData=NULL,
                p=p1, a=a1, pt=p2, at=a2, SSB.reduce=SSB.reduce, CIRCLES.d=CIRCLES.d,
                VIF=VIF, COR=COR,
                threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence,
                PLOTS=PLOTS, CATCH.OFF=CATCH.OFF,
                evaluations.keep=T, models.keep=F,
                ENSEMBLE.tune=ENSEMBLE.tune,
                ENSEMBLE.best=ENSEMBLE.best, ENSEMBLE.min=ENSEMBLE.min, ENSEMBLE.exponent=ENSEMBLE.exponent, ENSEMBLE.weight.min=ENSEMBLE.weight.min,
                MAXENT=MAXENT, MAXLIKE=MAXLIKE, GBM=GBM, GBMSTEP=GBMSTEP, RF=RF, GLM=GLM, GLMSTEP=GLMSTEP, 
                GAM=GAM, GAMSTEP=GAMSTEP, MGCV=MGCV, MGCVFIX=MGCVFIX, EARTH=EARTH, RPART=RPART, 
                NNET=NNET, FDA=FDA, SVM=SVM, SVME=SVME, GLMNET=GLMNET,
                BIOCLIM.O=BIOCLIM.O, BIOCLIM=BIOCLIM, DOMAIN=DOMAIN, MAHAL=MAHAL, MAHAL01=MAHAL01,
                PROBIT=PROBIT,  
                Yweights=Yweights, 
                factors=factors2, dummy.vars=dummy.vars2,
                maxit=maxit,
                MAXENT.a=MAXENT.a, MAXENT.path=MAXENT.path,
                MAXLIKE.formula=MAXLIKE.formula, MAXLIKE.method=MAXLIKE.method,
                GBM.formula=GBM.formula, GBM.n.trees=GBM.n.trees, 
                GBMSTEP.gbm.x=GBMSTEP.gbm.x, GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
                GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction,
                GBMSTEP.step.size=GBMSTEP.step.size, 
                RF.formula=RF.formula, RF.ntree=RF.ntree, RF.mtry=RF.mtry, 
                GLM.formula=GLM.formula, GLM.family=GLM.family, 
                GLMSTEP.k=GLMSTEP.k, GLMSTEP.steps=GLMSTEP.steps, STEP.formula=STEP.formula, 
                GLMSTEP.scope=GLMSTEP.scope, 
                GAM.formula=GAM.formula, GAM.family=GAM.family, 
                GAMSTEP.steps=GAMSTEP.steps, GAMSTEP.scope=GAMSTEP.scope, GAMSTEP.pos=GAMSTEP.pos, 
                MGCV.formula=MGCV.formula, MGCV.select=MGCV.select, 
                MGCVFIX.formula=MGCVFIX.formula, 
                EARTH.formula=EARTH.formula, EARTH.glm=EARTH.glm, 
                RPART.formula=RPART.formula, RPART.xval=RPART.xval, 
                NNET.formula=NNET.formula, NNET.size=NNET.size, NNET.decay=NNET.decay,
                FDA.formula=FDA.formula, 
                SVM.formula=SVM.formula, 
                SVME.formula=SVME.formula,
                GLMNET.nlambda=GLMNET.nlambda, GLMNET.class=GLMNET.class,
                BIOCLIM.O.fraction=BIOCLIM.O.fraction,
                MAHAL.shape=MAHAL.shape)

        dummy.vars.noDOMAIN <- tests$evaluations$dummy.vars.noDOMAIN

        if(is.null(tests$evaluations$MAXENT.T)==F) {output["MAXENT",i] <- tests$evaluations$MAXENT.T@auc}
        if(is.null(tests$evaluations$MAXLIKE.T)==F) {output["MAXLIKE",i] <- tests$evaluations$MAXLIKE.T@auc}
        if(is.null(tests$evaluations$GBM.T)==F) {output["GBM",i] <- tests$evaluations$GBM.T@auc} 
        if(is.null(tests$evaluations$GBMSTEP.T)==F) {output["GBMSTEP",i] <- tests$evaluations$GBMSTEP.T@auc} 
        if(is.null(tests$evaluations$RF.T)==F) {output["RF",i] <- tests$evaluations$RF.T@auc}
        if(is.null(tests$evaluations$GLM.T)==F) {output["GLM",i] <- tests$evaluations$GLM.T@auc} 
        if(is.null(tests$evaluations$GLMS.T)==F) {output["GLMSTEP",i] <- tests$evaluations$GLMS.T@auc}
        if(is.null(tests$evaluations$GAM.T)==F) {output["GAM",i] <- tests$evaluations$GAM.T@auc} 
        if(is.null(tests$evaluations$GAMS.T)==F) {output["GAMSTEP",i] <- tests$evaluations$GAMS.T@auc}
        if(is.null(tests$evaluations$MGCV.T)==F) {output["MGCV",i] <- tests$evaluations$MGCV.T@auc} 
        if(is.null(tests$evaluations$MGCVF.T)==F) {output["MGCVFIX",i] <- tests$evaluations$MGCVF.T@auc} 
        if(is.null(tests$evaluations$EARTH.T)==F) {output["EARTH",i] <- tests$evaluations$EARTH.T@auc} 
        if(is.null(tests$evaluations$RPART.T)==F) {output["RPART",i] <- tests$evaluations$RPART.T@auc}
        if(is.null(tests$evaluations$NNET.T)==F) {output["NNET",i] <- tests$evaluations$NNET.T@auc} 
        if(is.null(tests$evaluations$FDA.T)==F) {output["FDA",i] <- tests$evaluations$FDA.T@auc}
        if(is.null(tests$evaluations$SVM.T)==F) {output["SVM",i] <- tests$evaluations$SVM.T@auc}
        if(is.null(tests$evaluations$SVME.T)==F) {output["SVME",i] <- tests$evaluations$SVME.T@auc}
        if(is.null(tests$evaluations$GLMNET.T)==F) {output["GLMNET",i] <- tests$evaluations$GLMNET.T@auc}
        if(is.null(tests$evaluations$BIOCLIM.O.T)==F) {output["BIOCLIM.O",i] <- tests$evaluations$BIOCLIM.O.T@auc}
        if(is.null(tests$evaluations$BIOCLIM.T)==F) {output["BIOCLIM",i] <- tests$evaluations$BIOCLIM.T@auc}
        if(is.null(tests$evaluations$DOMAIN.T)==F) {output["DOMAIN",i] <- tests$evaluations$DOMAIN.T@auc}
        if(is.null(tests$evaluations$MAHAL.T)==F) {output["MAHAL",i] <- tests$evaluations$MAHAL.T@auc}
        if(is.null(tests$evaluations$MAHAL01.T)==F) {output["MAHAL01",i] <- tests$evaluations$MAHAL01.T@auc}
        if(is.null(tests$evaluations$ENSEMBLE.T)==F) {output["ENSEMBLE",i] <- tests$evaluations$ENSEMBLE.T@auc}

        if(ENSEMBLE.tune == T) {
            output["MAXENT",k+1+i] <- tests$evaluations$STRATEGY.weights["MAXENT"]
            output["MAXLIKE",k+1+i] <- tests$evaluations$STRATEGY.weights["MAXLIKE"]
            output["GBM",k+1+i] <- tests$evaluations$STRATEGY.weights["GBM"]
            output["GBMSTEP",k+1+i] <- tests$evaluations$STRATEGY.weights["GBMSTEP"]
            output["RF",k+1+i] <- tests$evaluations$STRATEGY.weights["RF"]
            output["GLM",k+1+i] <- tests$evaluations$STRATEGY.weights["GLM"]
            output["GLMSTEP",k+1+i] <- tests$evaluations$STRATEGY.weights["GLMSTEP"]
            output["GAM",k+1+i] <- tests$evaluations$STRATEGY.weights["GAM"]
            output["GAMSTEP",k+1+i] <- tests$evaluations$STRATEGY.weights["GAMSTEP"]
            output["MGCV",k+1+i] <- tests$evaluations$STRATEGY.weights["MGCV"]
            output["MGCVFIX",k+1+i] <- tests$evaluations$STRATEGY.weights["MGCVFIX"]
            output["EARTH",k+1+i] <- tests$evaluations$STRATEGY.weights["EARTH"]
            output["RPART",k+1+i] <- tests$evaluations$STRATEGY.weights["RPART"]
            output["NNET",k+1+i] <- tests$evaluations$STRATEGY.weights["NNET"]
            output["FDA",k+1+i] <- tests$evaluations$STRATEGY.weights["FDA"]
            output["SVM",k+1+i] <- tests$evaluations$STRATEGY.weights["SVM"]
            output["SVME",k+1+i] <- tests$evaluations$STRATEGY.weights["SVME"]
            output["GLMNET",k+1+i] <- tests$evaluations$STRATEGY.weights["GLMNET"]
            output["BIOCLIM.O",k+1+i] <- tests$evaluations$STRATEGY.weights["BIOCLIM.O"]
            output["BIOCLIM",k+1+i] <- tests$evaluations$STRATEGY.weights["BIOCLIM"]
            output["DOMAIN",k+1+i] <- tests$evaluations$STRATEGY.weights["DOMAIN"]
            output["MAHAL",k+1+i] <- tests$evaluations$STRATEGY.weights["MAHAL"]
            output["MAHAL01",k+1+i] <- tests$evaluations$STRATEGY.weights["MAHAL01"]
        }

        TestData.all[[i]] <- tests$evaluations$TestData
    }

    output[,k+1] <- rowMeans(output[,c(1:k)], na.rm=T)
    output[is.na(output[,k+1]),(k+1)] <- 0
#
# Try to use exponent, min and best to calculate final weights
#
# in case there were several exponents, then do not change the weights
        ENSEMBLE.exponent1 <- ENSEMBLE.exponent
        if (length(ENSEMBLE.exponent) > 1) {ENSEMBLE.exponent1 <- 1}
        ENSEMBLE.min1 <- ENSEMBLE.min
        if (length(ENSEMBLE.min) > 1) {ENSEMBLE.min1 <- 0.5}
        ENSEMBLE.best1 <- ENSEMBLE.best
        if (length(ENSEMBLE.best) > 1) {ENSEMBLE.best1 <- 0}
#
    if(ENSEMBLE.tune == T) {
        output[,2*k+2] <- rowMeans(output[,c((k+2):(2*k+1))], na.rm=T)
        output[is.na(output[,2*k+2]),(2*k+2)] <- 0
    }
#
    output.weights <- output[, "MEAN"]
    output.weights <- output.weights[names(output.weights) != "ENSEMBLE"]
    output <- output[order(output[,k+1], decreasing=T),]
    cat(paste("Results of ensemble.calibrate.weights sorted by average AUC for tests T_1 to T_", k, "\n", sep = ""))
    cat(paste("\n", "columns T_1 to T_", k, " show the AUC for each ", k, "-fold cross-validation run", "\n", sep = ""))
    if(ENSEMBLE.tune == T) {
        cat(paste("column MEAN.T shows the mean of the AUC", "\n\n", sep = ""))
        cat(paste("columns S_1 to S_", k, " show the weights for the ensemble model with best AUC", "\n", sep = ""))
        cat(paste("column MEAN shows the mean of these weights", "\n\n", sep = ""))
    }else{
        cat(paste("column MEAN shows the mean of the AUC", "\n\n", sep = ""))
    }
    print(output)

    cat(paste("\n", "input weights for ensemble modelling based on MEAN column", "\n",  sep = ""))
    print(output.weights)

    if(ENSEMBLE.tune == F) {
# also downweight the average AUC with the same parameters
        cat(paste("\n\n", "parameters for next weighting: ENSEMBLE.min=", ENSEMBLE.min1, ", ENSEMBLE.best=", ENSEMBLE.best1, " and ENSEMBLE.exponent=", ENSEMBLE.exponent1, "\n\n", sep = ""))
        output.weights <- ensemble.weights(weights=output.weights, exponent=ENSEMBLE.exponent1, best=ENSEMBLE.best1, min.weight=ENSEMBLE.min1)
        print(output.weights)
    }else{
# if possible, select best models (models with highest weights)
        if (ENSEMBLE.best1 > 0) {
            cat(paste("\n\n", "parameters for next weighting: ENSEMBLE.min=0, ENSEMBLE.best=", ENSEMBLE.best1, " and ENSEMBLE.exponent=1", "\n\n", sep = ""))
            output.weights <- ensemble.weights(weights=output.weights, exponent=1, best=ENSEMBLE.best1, min.weight=0)
            print(output.weights)
        }
    }

# remove models with low input weights (mainly to reduce the number of models for final calibrations and mapping)
    cat(paste("\n", "Minimum input weight is ", ENSEMBLE.weight.min, "\n", sep=""))
    output.weights2 <- output.weights
    while(min(output.weights2) < ENSEMBLE.weight.min) {
        output.weights2 <- output.weights2[-which.min(output.weights2)]
        output.weights2 <- ensemble.weights(weights=output.weights2, exponent=1, best=0, min.weight=0)
    }
    output.weights[] <- 0
    for (i in 1:length(output.weights2)) {output.weights[which(names(output.weights) == names(output.weights2)[i])] <- output.weights2[i]}
    cat(paste("\n", "final suggested weights for ensemble forecasting", "\n", sep = ""))
    print(output.weights)

# test with suggested final weights
    output2 <- numeric(length=k+1)
    names(output2)[1:k] <- paste("T_", c(1:k), sep="")
    names(output2)[k+1] <- c("MEAN.T")
    for (i in 1:k) {
        TestData <- TestData.all[[i]]
        TestData[,"ENSEMBLE"] <- output.weights["MAXENT"]*TestData[,"MAXENT"] + output.weights["MAXLIKE"]*TestData[,"MAXLIKE"] + output.weights["GBM"]*TestData[,"GBM"] +
            output.weights["GBMSTEP"]*TestData[,"GBMSTEP"] + output.weights["RF"]*TestData[,"RF"] + output.weights["GLM"]*TestData[,"GLM"] +
            output.weights["GLMSTEP"]*TestData[,"GLMSTEP"] + output.weights["GAM"]*TestData[,"GAM"] + output.weights["GAMSTEP"]*TestData[,"GAMSTEP"] +
            output.weights["MGCV"]*TestData[,"MGCV"] + output.weights["MGCVFIX"]*TestData[,"MGCVFIX"] + output.weights["EARTH"]*TestData[,"EARTH"] +
            output.weights["RPART"]*TestData[,"RPART"] + output.weights["NNET"]*TestData[,"NNET"] + output.weights["FDA"]*TestData[,"FDA"] +
            output.weights["SVM"]*TestData[,"SVM"] + output.weights["SVME"]*TestData[,"SVME"] + output.weights["GLMNET"]*TestData[,"GLMNET"] + 
            output.weights["BIOCLIM.O"]*TestData[,"BIOCLIM.O"] + output.weights["BIOCLIM"]*TestData[,"BIOCLIM"] +
            output.weights["DOMAIN"]*TestData[,"DOMAIN"] + output.weights["MAHAL"]*TestData[,"MAHAL"]+ output.weights["MAHAL01"]*TestData[,"MAHAL01"]
        eval1 <- eval2 <- NULL
        TestPres <- as.numeric(TestData[TestData[,"pb"]==1, "ENSEMBLE"])
        TestAbs <- as.numeric(TestData[TestData[,"pb"]==0, "ENSEMBLE"])
        eval1 <- dismo::evaluate(p=TestPres, a=TestAbs)
        output2[i] <- eval1@auc
    }
    output2[k+1] <- mean(output2[1:k])
    cat(paste("\n", "AUC for ensemble models based on suggested input weights (using presence and background data sets generated for ", k, "-fold cross-validations)",  "\n", sep = ""))
    print(output2)

    if (SINK==T && OLD.SINK==F) {sink(file=NULL, append=T)}
    if (data.keep == F) {
        cat(paste("\n\n"))
        return(list(AUC.table=output, table=output, output.weights=output.weights, AUC.with.suggested.weights=output2, 
            data=TestData.all, 
            x=x, p=p.all, a=a.all, MAXENT.a=MAXENT.a,
            var.names=var.names, factors=factors2, dummy.vars=dummy.vars2, dummy.vars.noDOMAIN=dummy.vars.noDOMAIN,
            species.name=species.name, 
            call=match.call()))
    }else{
        cat(paste("\n\n"))
        return(list(data=TestData.all, 
            AUC.table=output, table=output, output.weights=output.weights, AUC.with.suggested.weights=output2, 
            data=TestData.all, 
            x=x, p=p.all, a=a.all, MAXENT.a=MAXENT.a,
            var.names=var.names, factors=factors2, dummy.vars=dummy.vars2, dummy.vars.noDOMAIN=dummy.vars.noDOMAIN,
            species.name=species.name, 
            call=match.call()))
    }
}

