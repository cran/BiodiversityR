`ensemble.drop1` <- function(
    x=NULL, p=NULL, 
    a=NULL, an=1000, excludep=FALSE, target.groups=FALSE,
    k=0, pt=NULL, at=NULL, SSB.reduce=FALSE, CIRCLES.d=100000,
    TrainData=NULL, TestData=NULL,
    VIF=FALSE, COR=FALSE,
    SINK=FALSE, species.name="Species001",
    difference=FALSE, variables.alone=FALSE,
    ENSEMBLE.tune=FALSE,
    ENSEMBLE.best=0, ENSEMBLE.min=0.7, ENSEMBLE.exponent=1,
    input.weights=NULL,
    MAXENT=1, MAXNET=1, MAXLIKE=1, GBM=1, GBMSTEP=0, RF=1, CF=1,
    GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=1, MGCV=1, 
    MGCVFIX=0, EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, GLMNET=1,
    BIOCLIM.O=0, BIOCLIM=1, DOMAIN=1, MAHAL=1, MAHAL01=1,
    PROBIT=FALSE,
    Yweights="BIOMOD", 
    layer.drops=NULL, factors=NULL, dummy.vars=NULL,
    maxit=100,
    MAXENT.a=NULL, MAXENT.an=10000, 
    MAXENT.path=paste(getwd(), "/models/maxent_", species.name,  sep=""), 
    MAXNET.classes="default", MAXNET.clamp=FALSE, MAXNET.type="cloglog",
    MAXLIKE.method="BFGS",
    GBM.n.trees=2001, 
    GBMSTEP.tree.complexity=5, GBMSTEP.learning.rate=0.005, 
    GBMSTEP.bag.fraction=0.5, GBMSTEP.step.size=100, 
    RF.ntree=751, 
    CF.ntree=751,
    GLM.family=binomial(link="logit"), 
    GLMSTEP.steps=1000, GLMSTEP.scope=NULL, GLMSTEP.k=2, 
    GAM.family=binomial(link="logit"), 
    GAMSTEP.steps=1000, GAMSTEP.scope=NULL, GAMSTEP.pos=1, 
    MGCV.select=FALSE,  
    EARTH.glm=list(family=binomial(link="logit"), maxit=maxit), 
    RPART.xval=50, 
    NNET.size=8, NNET.decay=0.01, 
    GLMNET.nlambda=100, GLMNET.class=FALSE,
    BIOCLIM.O.fraction=0.9,
    MAHAL.shape=1
)
{

    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if (MAXLIKE > 0) {
        cat(paste("\n", "WARNING: MAXLIKE algorithm will not be implemented as MAXLIKE does not accept (new) data.frames as input", "\n", sep = ""))
        MAXLIKE <- 0
    }
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
        cat(paste("\n\n", "RESULTS (ensemble.drop1 function)", "\n\n", sep=""), file=paste.file, append=T)
        sink(file=paste.file, append=T)
        cat(paste(date(), "\n", sep=""))
        print(match.call())
    }

# estimate deviance
    loglik.calculation <- function(obs=NULL, preds=NULL) {
        preds[preds < 0.0000000001] <- 0.0000000001
        preds[preds > 0.9999999999] <- 0.9999999999
        out <- dismo::calc.deviance(obs=obs, pred=preds, calc.mean=F)
        return(out)
    }

#
# first fit with all variables

    if (raster::nlayers(x) == 0) {
# MAXENT needs x to make MAXENT.TrainData, MAXLIKE can only be calibrated with x
        MAXENT <- 0
        MAXLIKE <- 0
    }

    cat(paste("\n\n", "RESULTS WITH ALL VARIABLES", "\n\n", sep=""))

    tests <- ensemble.calibrate.models(x=x, 
        p=p, a=a, an=an, excludep=excludep, target.groups=target.groups,
        k=k, pt=pt, at=at, SSB.reduce=SSB.reduce, CIRCLES.d=CIRCLES.d,
        TrainData=NULL, TestData=NULL,
        PLOTS=FALSE, evaluations.keep=T, models.keep=F,
        VIF=VIF, COR=COR,
        formulae.defaults=T, maxit=maxit,
        ENSEMBLE.tune=ENSEMBLE.tune,
        ENSEMBLE.best=ENSEMBLE.best, ENSEMBLE.min=ENSEMBLE.min,
        ENSEMBLE.exponent=ENSEMBLE.exponent,
        input.weights=input.weights,
        MAXENT=MAXENT, MAXNET=MAXNET, MAXLIKE=MAXLIKE, GBM=GBM, GBMSTEP=GBMSTEP, RF=RF, CF=CF,
        GLM=GLM, GLMSTEP=GLMSTEP, 
        GAM=GAM, GAMSTEP=GAMSTEP, MGCV=MGCV, MGCVFIX=MGCVFIX, EARTH=EARTH, RPART=RPART, 
        NNET=NNET, FDA=FDA, SVM=SVM, SVME=SVME, GLMNET=GLMNET,
        BIOCLIM.O=BIOCLIM.O, BIOCLIM=BIOCLIM, DOMAIN=DOMAIN, 
        MAHAL=MAHAL, MAHAL01=MAHAL01,
        PROBIT=PROBIT,
        Yweights=Yweights, 
        layer.drops=layer.drops, factors=factors, dummy.vars=dummy.vars,
        MAXENT.a=MAXENT.a, MAXENT.an=MAXENT.an, MAXENT.path=MAXENT.path,
        MAXNET.classes=MAXNET.classes, MAXNET.clamp=MAXNET.clamp, MAXNET.type=MAXNET.type,
        MAXLIKE.formula=NULL, MAXLIKE.method=MAXLIKE.method,
        GBM.formula=NULL, GBM.n.trees=GBM.n.trees,
        GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
        GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction,
        GBMSTEP.step.size=GBMSTEP.step.size,
        RF.formula=NULL, RF.ntree=RF.ntree, 
        CF.formula=NULL, CF.ntree=CF.ntree, 
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
        GLMNET.nlambda=GLMNET.nlambda, GLMNET.class=GLMNET.class,
        BIOCLIM.O.fraction=BIOCLIM.O.fraction,
        MAHAL.shape=MAHAL.shape)

# use output to get names of the variables 
    var.names <- tests$evaluations$var.names
    nv <- length(var.names)

# get locations for MAXLIKE
    p1 <- tests$evaluations$p
    p2 <- tests$evaluations$pt
    a1 <- tests$evaluations$a
    a2 <- tests$evaluations$at

    model.names <- c("MAXENT", "MAXNET", "MAXLIKE", "GBM", "GBMSTEP", "RF", "CF",
        "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", 
        "MGCVFIX", "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "GLMNET",
        "BIOCLIM.O", "BIOCLIM", "DOMAIN", "MAHAL", "MAHAL01", "ENSEMBLE")

    output.C <- array(NA, dim=c(length(model.names), nv+1))
    rownames(output.C) <- model.names
    colnames(output.C) <- c("all_vars", paste("without_", var.names, sep=""))

    output.T <- array(NA, dim=c(length(model.names), nv+1))
    rownames(output.T) <- model.names
    colnames(output.T) <- c("all_vars", paste("without_", var.names, sep=""))

    output.LLC <- array(NA, dim=c(length(model.names), nv+1))
    rownames(output.LLC) <- model.names
    colnames(output.LLC) <- c("all_vars", paste("without_", var.names, sep=""))

    output.LLT <- array(NA, dim=c(length(model.names), nv+1))
    rownames(output.LLT) <- model.names
    colnames(output.LLT) <- c("all_vars", paste("without_", var.names, sep=""))

    if(is.null(tests$evaluations$MAXENT.C)==F) {output.C["MAXENT",1] <- tests$evaluations$MAXENT.C@auc}
    if(is.null(tests$evaluations$MAXENT.T)==F) {output.T["MAXENT",1] <- tests$evaluations$MAXENT.T@auc}
    if(sum(tests$evaluations$TrainData$MAXENT) > 0) {output.LLC["MAXENT",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MAXENT)}
    if(sum(tests$evaluations$TestData$MAXENT) > 0) {output.LLT["MAXENT",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MAXENT)}

    if(is.null(tests$evaluations$MAXNET.C)==F) {output.C["MAXNET",1] <- tests$evaluations$MAXNET.C@auc}
    if(is.null(tests$evaluations$MAXNET.T)==F) {output.T["MAXNET",1] <- tests$evaluations$MAXNET.T@auc}
    if(sum(tests$evaluations$TrainData$MAXNET) > 0) {output.LLC["MAXNET",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MAXNET)}
    if(sum(tests$evaluations$TestData$MAXNET) > 0) {output.LLT["MAXNET",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MAXNET)}

    if(is.null(tests$evaluations$MAXLIKE.C)==F) {output.C["MAXLIKE",1] <- tests$evaluations$MAXLIKE.C@auc}
    if(is.null(tests$evaluations$MAXLIKE.T)==F) {output.T["MAXLIKE",1] <- tests$evaluations$MAXLIKE.T@auc}
    if(sum(tests$evaluations$TrainData$MAXLIKE) > 0) {output.LLC["MAXLIKE",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MAXLIKE)}
    if(sum(tests$evaluations$TestData$MAXLIKE) > 0) {output.LLT["MAXLIKE",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MAXLIKE)}

    if(is.null(tests$evaluations$GBM.C)==F) {output.C["GBM",1] <- tests$evaluations$GBM.C@auc} 
    if(is.null(tests$evaluations$GBM.T)==F) {output.T["GBM",1] <- tests$evaluations$GBM.T@auc} 
    if(sum(tests$evaluations$TrainData$GBM) > 0) {output.LLC["GBM",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GBM)}
    if(sum(tests$evaluations$TestData$GBM) > 0) {output.LLT["GBM",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GBM)}

    if(is.null(tests$evaluations$GBMSTEP.C)==F) {output.C["GBMSTEP",1] <- tests$evaluations$GBMSTEP.C@auc} 
    if(is.null(tests$evaluations$GBMSTEP.T)==F) {output.T["GBMSTEP",1] <- tests$evaluations$GBMSTEP.T@auc} 
    if(sum(tests$evaluations$TrainData$GBMSTEP) > 0) {output.LLC["GBMSTEP",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GBMSTEP)}
    if(sum(tests$evaluations$TestData$GBMSTEP) > 0) {output.LLT["GBMSTEP",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GBMSTEP)}

    if(is.null(tests$evaluations$RF.C)==F) {output.C["RF",1] <- tests$evaluations$RF.C@auc}
    if(is.null(tests$evaluations$RF.T)==F) {output.T["RF",1] <- tests$evaluations$RF.T@auc}
    if(sum(tests$evaluations$TrainData$RF) > 0) {output.LLC["RF",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$RF)}
    if(sum(tests$evaluations$TestData$RF) > 0) {output.LLT["RF",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$RF)}

    if(is.null(tests$evaluations$CF.C)==F) {output.C["CF",1] <- tests$evaluations$CF.C@auc}
    if(is.null(tests$evaluations$CF.T)==F) {output.T["CF",1] <- tests$evaluations$CF.T@auc}
    if(sum(tests$evaluations$TrainData$CF) > 0) {output.LLC["CF",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$CF)}
    if(sum(tests$evaluations$TestData$CF) > 0) {output.LLT["CF",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$CF)}

    if(is.null(tests$evaluations$GLM.C)==F) {output.C["GLM",1] <- tests$evaluations$GLM.C@auc} 
    if(is.null(tests$evaluations$GLM.T)==F) {output.T["GLM",1] <- tests$evaluations$GLM.T@auc} 
    if(sum(tests$evaluations$TrainData$GLM) > 0) {output.LLC["GLM",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GLM)}
    if(sum(tests$evaluations$TestData$GLM) > 0) {output.LLT["GLM",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GLM)}

    if(is.null(tests$evaluations$GLMS.C)==F) {output.C["GLMSTEP",1] <- tests$evaluations$GLMS.C@auc}
    if(is.null(tests$evaluations$GLMS.T)==F) {output.T["GLMSTEP",1] <- tests$evaluations$GLMS.T@auc}
    if(sum(tests$evaluations$TrainData$GLMSTEP) > 0) {output.LLC["GLMSTEP",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GLMSTEP)}
    if(sum(tests$evaluations$TestData$GLMSTEP) > 0) {output.LLT["GLMSTEP",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GLMSTEP)}

    if(is.null(tests$evaluations$GAM.C)==F) {output.C["GAM",1] <- tests$evaluations$GAM.C@auc} 
    if(is.null(tests$evaluations$GAM.T)==F) {output.T["GAM",1] <- tests$evaluations$GAM.T@auc} 
    if(sum(tests$evaluations$TrainData$GAM) > 0) {output.LLC["GAM",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GAM)}
    if(sum(tests$evaluations$TestData$GAM) > 0) {output.LLT["GAM",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GAM)}

    if(is.null(tests$evaluations$GAMS.C)==F) {output.C["GAMSTEP",1] <- tests$evaluations$GAMS.C@auc}
    if(is.null(tests$evaluations$GAMS.T)==F) {output.T["GAMSTEP",1] <- tests$evaluations$GAMS.T@auc}
    if(sum(tests$evaluations$TrainData$GAMSTEP) > 0) {output.LLC["GAMSTEP",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GAMSTEP)}
    if(sum(tests$evaluations$TestData$GAMSTEP) > 0) {output.LLT["GAMSTEP",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GAMSTEP)}

    if(is.null(tests$evaluations$MGCV.C)==F) {output.C["MGCV",1] <- tests$evaluations$MGCV.C@auc} 
    if(is.null(tests$evaluations$MGCV.T)==F) {output.T["MGCV",1] <- tests$evaluations$MGCV.T@auc} 
    if(sum(tests$evaluations$TrainData$MGCV) > 0) {output.LLC["MGCV",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MGCV)}
    if(sum(tests$evaluations$TestData$MGCV) > 0) {output.LLT["MGCV",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MGCV)}

    if(is.null(tests$evaluations$MGCVF.C)==F) {output.C["MGCVFIX",1] <- tests$evaluations$MGCVF.C@auc} 
    if(is.null(tests$evaluations$MGCVF.T)==F) {output.T["MGCVFIX",1] <- tests$evaluations$MGCVF.T@auc} 
    if(sum(tests$evaluations$TrainData$MGCVFIX) > 0) {output.LLC["MGCVFIX",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MGCVFIX)}
    if(sum(tests$evaluations$TestData$MGCVFIX) > 0) {output.LLT["MGCVFIX",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MGCVFIX)}

    if(is.null(tests$evaluations$EARTH.C)==F) {output.C["EARTH",1] <- tests$evaluations$EARTH.C@auc} 
    if(is.null(tests$evaluations$EARTH.T)==F) {output.T["EARTH",1] <- tests$evaluations$EARTH.T@auc} 
    if(sum(tests$evaluations$TrainData$EARTH) > 0) {output.LLC["EARTH",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$EARTH)}
    if(sum(tests$evaluations$TestData$EARTH) > 0) {output.LLT["EARTH",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$EARTH)}

    if(is.null(tests$evaluations$RPART.C)==F) {output.C["RPART",1] <- tests$evaluations$RPART.C@auc}
    if(is.null(tests$evaluations$RPART.T)==F) {output.T["RPART",1] <- tests$evaluations$RPART.T@auc}
    if(sum(tests$evaluations$TrainData$RPART) > 0) {output.LLC["RPART",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$RPART)}
    if(sum(tests$evaluations$TestData$RPART) > 0) {output.LLT["RPART",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$RPART)}

    if(is.null(tests$evaluations$NNET.C)==F) {output.C["NNET",1] <- tests$evaluations$NNET.C@auc} 
    if(is.null(tests$evaluations$NNET.T)==F) {output.T["NNET",1] <- tests$evaluations$NNET.T@auc} 
    if(sum(tests$evaluations$TrainData$NNET) > 0) {output.LLC["NNET",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$NNET)}
    if(sum(tests$evaluations$TestData$NNET) > 0) {output.LLT["NNET",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$NNET)}

    if(is.null(tests$evaluations$FDA.C)==F) {output.C["FDA",1] <- tests$evaluations$FDA.C@auc}
    if(is.null(tests$evaluations$FDA.T)==F) {output.T["FDA",1] <- tests$evaluations$FDA.T@auc}
    if(sum(tests$evaluations$TrainData$FDA) > 0) {output.LLC["FDA",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$FDA)}
    if(sum(tests$evaluations$TestData$FDA) > 0) {output.LLT["FDA",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$FDA)}

    if(is.null(tests$evaluations$SVM.C)==F) {output.C["SVM",1] <- tests$evaluations$SVM.C@auc}
    if(is.null(tests$evaluations$SVM.T)==F) {output.T["SVM",1] <- tests$evaluations$SVM.T@auc}
    if(sum(tests$evaluations$TrainData$SVM) > 0) {output.LLC["SVM",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$SVM)}
    if(sum(tests$evaluations$TestData$SVM) > 0) {output.LLT["SVM",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$SVM)}

    if(is.null(tests$evaluations$SVME.C)==F) {output.C["SVME",1] <- tests$evaluations$SVME.C@auc}
    if(is.null(tests$evaluations$SVME.T)==F) {output.T["SVME",1] <- tests$evaluations$SVME.T@auc}
    if(sum(tests$evaluations$TrainData$SVME) > 0) {output.LLC["SVME",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$SVME)}
    if(sum(tests$evaluations$TestData$SVME) > 0) {output.LLT["SVME",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$SVME)}

    if(is.null(tests$evaluations$GLMNET.C)==F) {output.C["GLMNET",1] <- tests$evaluations$GLMNET.C@auc}
    if(is.null(tests$evaluations$GLMNET.T)==F) {output.T["GLMNET",1] <- tests$evaluations$GLMNET.T@auc}
    if(sum(tests$evaluations$TrainData$GLMNET) > 0) {output.LLC["GLMNET",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GLMNET)}
    if(sum(tests$evaluations$TestData$GLMNET) > 0) {output.LLT["GLMNET",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GLMNET)}

    if(is.null(tests$evaluations$BIOCLIM.O.C)==F) {output.C["BIOCLIM.O",1] <- tests$evaluations$BIOCLIM.O.C@auc}
    if(is.null(tests$evaluations$BIOCLIM.O.T)==F) {output.T["BIOCLIM.O",1] <- tests$evaluations$BIOCLIM.O.T@auc}
    if(sum(tests$evaluations$TrainData$BIOCLIM.O) > 0) {output.LLC["BIOCLIM.O",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$BIOCLIM.O)}
    if(sum(tests$evaluations$TestData$BIOCLIM.O) > 0) {output.LLT["BIOCLIM.O",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$BIOCLIM.O)}

    if(is.null(tests$evaluations$BIOCLIM.C)==F) {output.C["BIOCLIM",1] <- tests$evaluations$BIOCLIM.C@auc}
    if(is.null(tests$evaluations$BIOCLIM.T)==F) {output.T["BIOCLIM",1] <- tests$evaluations$BIOCLIM.T@auc}
    if(sum(tests$evaluations$TrainData$BIOCLIM) > 0) {output.LLC["BIOCLIM",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$BIOCLIM)}
    if(sum(tests$evaluations$TestData$BIOCLIM) > 0) {output.LLT["BIOCLIM",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$BIOCLIM)}

    if(is.null(tests$evaluations$DOMAIN.C)==F) {output.C["DOMAIN",1] <- tests$evaluations$DOMAIN.C@auc}
    if(is.null(tests$evaluations$DOMAIN.T)==F) {output.T["DOMAIN",1] <- tests$evaluations$DOMAIN.T@auc}
    if(sum(tests$evaluations$TrainData$DOMAIN) > 0) {output.LLC["DOMAIN",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$DOMAIN)}
    if(sum(tests$evaluations$TestData$DOMAIN) > 0) {output.LLT["DOMAIN",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$DOMAIN)}

    if(is.null(tests$evaluations$MAHAL.C)==F) {output.C["MAHAL",1] <- tests$evaluations$MAHAL.C@auc}
    if(is.null(tests$evaluations$MAHAL.T)==F) {output.T["MAHAL",1] <- tests$evaluations$MAHAL.T@auc}
    if(sum(tests$evaluations$TrainData$MAHAL) > 0) {output.LLC["MAHAL",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MAHAL)}
    if(sum(tests$evaluations$TestData$MAHAL) > 0) {output.LLT["MAHAL",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MAHAL)}

    if(is.null(tests$evaluations$MAHAL01.C)==F) {output.C["MAHAL01",1] <- tests$evaluations$MAHAL01.C@auc}
    if(is.null(tests$evaluations$MAHAL01.T)==F) {output.T["MAHAL01",1] <- tests$evaluations$MAHAL01.T@auc}
    if(sum(tests$evaluations$TrainData$MAHAL01) > 0) {output.LLC["MAHAL01",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MAHAL01)}
    if(sum(tests$evaluations$TestData$MAHAL01) > 0) {output.LLT["MAHAL01",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MAHAL01)}

    if(is.null(tests$evaluations$ENSEMBLE.C)==F) {output.C["ENSEMBLE",1] <- tests$evaluations$ENSEMBLE.C@auc}
    if(is.null(tests$evaluations$ENSEMBLE.T)==F) {output.T["ENSEMBLE",1] <- tests$evaluations$ENSEMBLE.T@auc}
    if(sum(tests$evaluations$TrainData$ENSEMBLE) > 0) {output.LLC["ENSEMBLE",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$ENSEMBLE)}
    if(sum(tests$evaluations$TestData$ENSEMBLE) > 0) {output.LLT["ENSEMBLE",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$ENSEMBLE)}

# sequentially leave out the focal variable, then fit again
# the data sets are used from the "full" model

    var.names2 <- c("pb", var.names)
    TrainData1 <- tests$evaluations$TrainData
    TrainData1 <- TrainData1[, which(names(TrainData1) %in% var.names2), drop=F]

    glm1 <- glm(as.formula("pb~1"), data=TrainData1, family="binomial")
    preval.dev.cal <- glm1$deviance

# calculate worst possible model (predict absence where present and vice versa)
    numpres <- sum(TrainData1[, "pb"])
    numabs <- nrow(TrainData1) - numpres
    obs1 <- TrainData1[, "pb"]
    pred1 <- rep(0.000000001, numpres)
    pred2 <- rep(0.999999999, numabs)
    pred3 <- c(pred1, pred2)
    null.dev.cal <- loglik.calculation(obs=obs1, preds=pred3)

    TestData1 <- tests$evaluations$TestData
    TestData1 <- TestData1[, which(names(TestData1) %in% var.names2), drop=F]

    glm2 <- glm(as.formula("pb~1"), data=TestData1, family="binomial")
    preval.dev.test <- glm2$deviance

# calculate worst possible model (predict absence where present and vice versa)
    numpres <- sum(TestData1[, "pb"])
    numabs <- nrow(TestData1) - numpres
    obs1 <- TestData1[, "pb"]
    pred1 <- rep(0.000000001, numpres)
    pred2 <- rep(0.999999999, numabs)
    pred3 <- c(pred1, pred2)
    null.dev.test <- loglik.calculation(obs=obs1, preds=pred3)

    MAXENT.a <- tests$evaluations$MAXENT.a

    for (i in 1:nv) {
        var.f <- var.names[i]
        cat(paste("\n", "2.",  i, ". Leaving out variable: ", var.f, "\n\n", sep = ""))
        TrainData2 <- TrainData1[, which(names(TrainData1) != var.f), drop=F]
        TestData2 <- TestData1[, which(names(TestData1) != var.f), drop=F]
        factors2 <- NULL
        if (is.null(factors) == F) {
            factors2 <- factors[which(factors != var.f)]
            if (length(factors2) == 0) {factors2 <- NULL}
        }
        dummy.vars2 <- NULL
        if (is.null(dummy.vars) == F) {
            dummy.vars2 <- dummy.vars[which(dummy.vars != var.f)]
            if (length(dummy.vars2) == 0) {dummy.vars2 <- NULL}
        }

        if (is.null(layer.drops) == T) {
            layer.drops2 <- var.f
        }else{
            layer.drops2 <- c(layer.drops, var.f)
        }

        x2 <- raster::dropLayer(x, which(names(x) %in% layer.drops2))
        x2 <- raster::stack(x2)

        if (raster::nlayers(x2) == 0) {
            MAXENT <- 0
            MAXLIKE <- 0
        }

    tests <- ensemble.calibrate.models(x=x2, 
        p=p1, a=a1, an=an, excludep=excludep, 
        k=k, pt=p2, at=a2, SSB.reduce=FALSE, CIRCLES.d=CIRCLES.d,
        TrainData=TrainData2, TestData=TestData2, 
        PLOTS=FALSE, evaluations.keep=T,  
        VIF=VIF, COR=COR,
        formulae.defaults=T, maxit=maxit,
        ENSEMBLE.tune=ENSEMBLE.tune,
        ENSEMBLE.best=ENSEMBLE.best, ENSEMBLE.min=ENSEMBLE.min,
        ENSEMBLE.exponent=ENSEMBLE.exponent,
        input.weights=input.weights,
        MAXENT=MAXENT, MAXNET=MAXNET, MAXLIKE=MAXLIKE, GBM=GBM, GBMSTEP=GBMSTEP, RF=RF, CF=CF,
        GLM=GLM, GLMSTEP=GLMSTEP, GAM=GAM, GAMSTEP=GAMSTEP, MGCV=MGCV, MGCVFIX=MGCVFIX, 
        EARTH=EARTH, RPART=RPART, NNET=NNET, FDA=FDA, SVM=SVM, SVME=SVME, GLMNET=GLMNET,
        BIOCLIM.O=BIOCLIM.O, BIOCLIM=BIOCLIM, DOMAIN=DOMAIN, 
        MAHAL=MAHAL, MAHAL01=MAHAL01,
        PROBIT=PROBIT,
        Yweights=Yweights, 
        layer.drops=layer.drops2, factors=factors2, dummy.vars=dummy.vars2,
        MAXENT.a=MAXENT.a, MAXENT.path=MAXENT.path,
        MAXNET.classes=MAXNET.classes, MAXNET.clamp=MAXNET.clamp, MAXNET.type=MAXNET.type,
        MAXLIKE.formula=NULL, MAXLIKE.method=MAXLIKE.method,
        GBM.formula=NULL, GBM.n.trees=GBM.n.trees,
        GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
        GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction,
        GBMSTEP.step.size=GBMSTEP.step.size,
        RF.formula=NULL, RF.ntree=RF.ntree, 
        CF.formula=NULL, CF.ntree=CF.ntree, 
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
        GLMNET.nlambda=GLMNET.nlambda,
        BIOCLIM.O.fraction=BIOCLIM.O.fraction,
        MAHAL.shape=MAHAL.shape)

    if(is.null(tests$evaluations$MAXENT.C)==F) {output.C["MAXENT",i+1] <- tests$evaluations$MAXENT.C@auc}
    if(is.null(tests$evaluations$MAXENT.T)==F) {output.T["MAXENT",i+1] <- tests$evaluations$MAXENT.T@auc}
    if(sum(tests$evaluations$TrainData$MAXENT) > 0) {output.LLC["MAXENT",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MAXENT)}
    if(sum(tests$evaluations$TestData$MAXENT) > 0) {output.LLT["MAXENT",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MAXENT)}

    if(is.null(tests$evaluations$MAXNET.C)==F) {output.C["MAXNET",i+1] <- tests$evaluations$MAXNET.C@auc}
    if(is.null(tests$evaluations$MAXNET.T)==F) {output.T["MAXNET",i+1] <- tests$evaluations$MAXNET.T@auc}
    if(sum(tests$evaluations$TrainData$MAXNET) > 0) {output.LLC["MAXNET",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MAXNET)}
    if(sum(tests$evaluations$TestData$MAXNET) > 0) {output.LLT["MAXNET",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MAXNET)}

    if(is.null(tests$evaluations$MAXLIKE.C)==F) {output.C["MAXLIKE",i+1] <- tests$evaluations$MAXLIKE.C@auc}
    if(is.null(tests$evaluations$MAXLIKE.T)==F) {output.T["MAXLIKE",i+1] <- tests$evaluations$MAXLIKE.T@auc}
    if(sum(tests$evaluations$TrainData$MAXLIKE) > 0) {output.LLC["MAXLIKE",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MAXLIKE)}
    if(sum(tests$evaluations$TestData$MAXLIKE) > 0) {output.LLT["MAXLIKE",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MAXLIKE)}

    if(is.null(tests$evaluations$GBM.C)==F) {output.C["GBM",i+1] <- tests$evaluations$GBM.C@auc} 
    if(is.null(tests$evaluations$GBM.T)==F) {output.T["GBM",i+1] <- tests$evaluations$GBM.T@auc} 
    if(sum(tests$evaluations$TrainData$GBM) > 0) {output.LLC["GBM",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GBM)}
    if(sum(tests$evaluations$TestData$GBM) > 0) {output.LLT["GBM",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GBM)}

    if(is.null(tests$evaluations$GBMSTEP.C)==F) {output.C["GBMSTEP",i+1] <- tests$evaluations$GBMSTEP.C@auc} 
    if(is.null(tests$evaluations$GBMSTEP.T)==F) {output.T["GBMSTEP",i+1] <- tests$evaluations$GBMSTEP.T@auc} 
    if(sum(tests$evaluations$TrainData$GBMSTEP) > 0) {output.LLC["GBMSTEP",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GBMSTEP)}
    if(sum(tests$evaluations$TestData$GBMSTEP) > 0) {output.LLT["GBMSTEP",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GBMSTEP)}

    if(is.null(tests$evaluations$RF.C)==F) {output.C["RF",i+1] <- tests$evaluations$RF.C@auc}
    if(is.null(tests$evaluations$RF.T)==F) {output.T["RF",i+1] <- tests$evaluations$RF.T@auc}
    if(sum(tests$evaluations$TrainData$RF) > 0) {output.LLC["RF",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$RF)}
    if(sum(tests$evaluations$TestData$RF) > 0) {output.LLT["RF",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$RF)}

    if(is.null(tests$evaluations$CF.C)==F) {output.C["CF",i+1] <- tests$evaluations$CF.C@auc}
    if(is.null(tests$evaluations$CF.T)==F) {output.T["CF",i+1] <- tests$evaluations$CF.T@auc}
    if(sum(tests$evaluations$TrainData$CF) > 0) {output.LLC["CF",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$CF)}
    if(sum(tests$evaluations$TestData$CF) > 0) {output.LLT["CF",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$CF)}

    if(is.null(tests$evaluations$GLM.C)==F) {output.C["GLM",i+1] <- tests$evaluations$GLM.C@auc} 
    if(is.null(tests$evaluations$GLM.T)==F) {output.T["GLM",i+1] <- tests$evaluations$GLM.T@auc} 
    if(sum(tests$evaluations$TrainData$GLM) > 0) {output.LLC["GLM",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GLM)}
    if(sum(tests$evaluations$TestData$GLM) > 0) {output.LLT["GLM",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GLM)}

    if(is.null(tests$evaluations$GLMS.C)==F) {output.C["GLMSTEP",i+1] <- tests$evaluations$GLMS.C@auc}
    if(is.null(tests$evaluations$GLMS.T)==F) {output.T["GLMSTEP",i+1] <- tests$evaluations$GLMS.T@auc}
    if(sum(tests$evaluations$TrainData$GLMSTEP) > 0) {output.LLC["GLMSTEP",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GLMSTEP)}
    if(sum(tests$evaluations$TestData$GLMSTEP) > 0) {output.LLT["GLMSTEP",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GLMSTEP)}

    if(is.null(tests$evaluations$GAM.C)==F) {output.C["GAM",i+1] <- tests$evaluations$GAM.C@auc} 
    if(is.null(tests$evaluations$GAM.T)==F) {output.T["GAM",i+1] <- tests$evaluations$GAM.T@auc} 
    if(sum(tests$evaluations$TrainData$GAM) > 0) {output.LLC["GAM",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GAM)}
    if(sum(tests$evaluations$TestData$GAM) > 0) {output.LLT["GAM",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GAM)}

    if(is.null(tests$evaluations$GAMS.C)==F) {output.C["GAMSTEP",i+1] <- tests$evaluations$GAMS.C@auc}
    if(is.null(tests$evaluations$GAMS.T)==F) {output.T["GAMSTEP",i+1] <- tests$evaluations$GAMS.T@auc}
    if(sum(tests$evaluations$TrainData$GAMSTEP) > 0) {output.LLC["GAMSTEP",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GAMSTEP)}
    if(sum(tests$evaluations$TestData$GAMSTEP) > 0) {output.LLT["GAMSTEP",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GAMSTEP)}

    if(is.null(tests$evaluations$MGCV.C)==F) {output.C["MGCV",i+1] <- tests$evaluations$MGCV.C@auc} 
    if(is.null(tests$evaluations$MGCV.T)==F) {output.T["MGCV",i+1] <- tests$evaluations$MGCV.T@auc} 
    if(sum(tests$evaluations$TrainData$MGCV) > 0) {output.LLC["MGCV",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MGCV)}
    if(sum(tests$evaluations$TestData$MGCV) > 0) {output.LLT["MGCV",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MGCV)}

    if(is.null(tests$evaluations$MGCVF.C)==F) {output.C["MGCVFIX",i+1] <- tests$evaluations$MGCVF.C@auc} 
    if(is.null(tests$evaluations$MGCVF.T)==F) {output.T["MGCVFIX",i+1] <- tests$evaluations$MGCVF.T@auc} 
    if(sum(tests$evaluations$TrainData$MGCVFIX) > 0) {output.LLC["MGCVFIX",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MGCVFIX)}
    if(sum(tests$evaluations$TestData$MGCVFIX) > 0) {output.LLT["MGCVFIX",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MGCVFIX)}

    if(is.null(tests$evaluations$EARTH.C)==F) {output.C["EARTH",i+1] <- tests$evaluations$EARTH.C@auc} 
    if(is.null(tests$evaluations$EARTH.T)==F) {output.T["EARTH",i+1] <- tests$evaluations$EARTH.T@auc} 
    if(sum(tests$evaluations$TrainData$EARTH) > 0) {output.LLC["EARTH",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$EARTH)}
    if(sum(tests$evaluations$TestData$EARTH) > 0) {output.LLT["EARTH",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$EARTH)}

    if(is.null(tests$evaluations$RPART.C)==F) {output.C["RPART",i+1] <- tests$evaluations$RPART.C@auc}
    if(is.null(tests$evaluations$RPART.T)==F) {output.T["RPART",i+1] <- tests$evaluations$RPART.T@auc}
    if(sum(tests$evaluations$TrainData$RPART) > 0) {output.LLC["RPART",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$RPART)}
    if(sum(tests$evaluations$TestData$RPART) > 0) {output.LLT["RPART",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$RPART)}

    if(is.null(tests$evaluations$NNET.C)==F) {output.C["NNET",i+1] <- tests$evaluations$NNET.C@auc} 
    if(is.null(tests$evaluations$NNET.T)==F) {output.T["NNET",i+1] <- tests$evaluations$NNET.T@auc} 
    if(sum(tests$evaluations$TrainData$NNET) > 0) {output.LLC["NNET",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$NNET)}
    if(sum(tests$evaluations$TestData$NNET) > 0) {output.LLT["NNET",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$NNET)}

    if(is.null(tests$evaluations$FDA.C)==F) {output.C["FDA",i+1] <- tests$evaluations$FDA.C@auc}
    if(is.null(tests$evaluations$FDA.T)==F) {output.T["FDA",i+1] <- tests$evaluations$FDA.T@auc}
    if(sum(tests$evaluations$TrainData$FDA) > 0) {output.LLC["FDA",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$FDA)}
    if(sum(tests$evaluations$TestData$FDA) > 0) {output.LLT["FDA",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$FDA)}

    if(is.null(tests$evaluations$SVM.C)==F) {output.C["SVM",i+1] <- tests$evaluations$SVM.C@auc}
    if(is.null(tests$evaluations$SVM.T)==F) {output.T["SVM",i+1] <- tests$evaluations$SVM.T@auc}
    if(sum(tests$evaluations$TrainData$SVM) > 0) {output.LLC["SVM",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$SVM)}
    if(sum(tests$evaluations$TestData$SVM) > 0) {output.LLT["SVM",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$SVM)}

    if(is.null(tests$evaluations$SVME.C)==F) {output.C["SVME",i+1] <- tests$evaluations$SVME.C@auc}
    if(is.null(tests$evaluations$SVME.T)==F) {output.T["SVME",i+1] <- tests$evaluations$SVME.T@auc}
    if(sum(tests$evaluations$TrainData$SVME) > 0) {output.LLC["SVME",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$SVME)}
    if(sum(tests$evaluations$TestData$SVME) > 0) {output.LLT["SVME",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$SVME)}

    if(is.null(tests$evaluations$GLMNET.C)==F) {output.C["GLMNET",i+1] <- tests$evaluations$GLMNET.C@auc}
    if(is.null(tests$evaluations$GLMNET.T)==F) {output.T["GLMNET",i+1] <- tests$evaluations$GLMNET.T@auc}
    if(sum(tests$evaluations$TrainData$GLMNET) > 0) {output.LLC["GLMNET",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GLMNET)}
    if(sum(tests$evaluations$TestData$GLMNET) > 0) {output.LLT["GLMNET",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GLMNET)}

    if(is.null(tests$evaluations$BIOCLIM.O.C)==F) {output.C["BIOCLIM.O",i+1] <- tests$evaluations$BIOCLIM.O.C@auc}
    if(is.null(tests$evaluations$BIOCLIM.O.T)==F) {output.T["BIOCLIM.O",i+1] <- tests$evaluations$BIOCLIM.O.T@auc}
    if(sum(tests$evaluations$TrainData$BIOCLIM.O) > 0) {output.LLC["BIOCLIM.O",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$BIOCLIM.O)}
    if(sum(tests$evaluations$TestData$BIOCLIM.O) > 0) {output.LLT["BIOCLIM.O",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$BIOCLIM.O)}

    if(is.null(tests$evaluations$BIOCLIM.C)==F) {output.C["BIOCLIM",i+1] <- tests$evaluations$BIOCLIM.C@auc}
    if(is.null(tests$evaluations$BIOCLIM.T)==F) {output.T["BIOCLIM",i+1] <- tests$evaluations$BIOCLIM.T@auc}
    if(sum(tests$evaluations$TrainData$BIOCLIM) > 0) {output.LLC["BIOCLIM",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$BIOCLIM)}
    if(sum(tests$evaluations$TestData$BIOCLIM) > 0) {output.LLT["BIOCLIM",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$BIOCLIM)}

    if(is.null(tests$evaluations$DOMAIN.C)==F) {output.C["DOMAIN",i+1] <- tests$evaluations$DOMAIN.C@auc}
    if(is.null(tests$evaluations$DOMAIN.T)==F) {output.T["DOMAIN",i+1] <- tests$evaluations$DOMAIN.T@auc}
    if(sum(tests$evaluations$TrainData$DOMAIN) > 0) {output.LLC["DOMAIN",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$DOMAIN)}
    if(sum(tests$evaluations$TestData$DOMAIN) > 0) {output.LLT["DOMAIN",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$DOMAIN)}

    if(is.null(tests$evaluations$MAHAL.C)==F) {output.C["MAHAL",i+1] <- tests$evaluations$MAHAL.C@auc}
    if(is.null(tests$evaluations$MAHAL.T)==F) {output.T["MAHAL",i+1] <- tests$evaluations$MAHAL.T@auc}
    if(sum(tests$evaluations$TrainData$MAHAL) > 0) {output.LLC["MAHAL",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MAHAL)}
    if(sum(tests$evaluations$TestData$MAHAL) > 0) {output.LLT["MAHAL",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MAHAL)}

    if(is.null(tests$evaluations$MAHAL01.C)==F) {output.C["MAHAL01",i+1] <- tests$evaluations$MAHAL01.C@auc}
    if(is.null(tests$evaluations$MAHAL01.T)==F) {output.T["MAHAL01",i+1] <- tests$evaluations$MAHAL01.T@auc}
    if(sum(tests$evaluations$TrainData$MAHAL01) > 0) {output.LLC["MAHAL01",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MAHAL01)}
    if(sum(tests$evaluations$TestData$MAHAL01) > 0) {output.LLT["MAHAL01",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MAHAL01)}

    if(is.null(tests$evaluations$ENSEMBLE.C)==F) {output.C["ENSEMBLE",i+1] <- tests$evaluations$ENSEMBLE.C@auc}
    if(is.null(tests$evaluations$ENSEMBLE.T)==F) {output.T["ENSEMBLE",i+1] <- tests$evaluations$ENSEMBLE.T@auc}
    if(sum(tests$evaluations$TrainData$ENSEMBLE) > 0) {output.LLC["ENSEMBLE",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$ENSEMBLE)}
    if(sum(tests$evaluations$TestData$ENSEMBLE) > 0) {output.LLT["ENSEMBLE",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$ENSEMBLE)}

    }


   if (variables.alone == T) {

## Models only using the focal variable

    output1.C <- output.C
    output1.T <- output.T
    output1.LLC <- output.LLC
    output1.LLT <- output.LLT

    colnames(output1.C) <- colnames(output1.T) <- colnames(output1.LLC) <- colnames(output1.LLT) <- c("all_vars", paste("only_", var.names, sep=""))

    for (i in 1:nv) {
        var.f <- var.names[i]
        var.f2 <- c("pb", var.f)
        cat(paste("\n", "3.", i, ". Only using variable: ", var.f, "\n\n", sep = ""))
        TrainData2 <- TrainData1[, which(names(TrainData1) %in% var.f2), drop=F]
        TestData2 <- TestData1[, which(names(TestData1) %in% var.f2), drop=F]

        factors2 <- NULL
        if (is.null(factors) == F) {
            factors2 <- factors[which(factors == var.f)]
            if (length(factors2) == 0) {factors2 <- NULL}
        }
        dummy.vars2 <- NULL
        if (is.null(dummy.vars) == F) {
            dummy.vars2 <- dummy.vars[which(dummy.vars == var.f)]
            if (length(dummy.vars2) == 0) {dummy.vars2 <- NULL}
        }

        var.f3 <- var.names[which(var.names != var.f)]

        if (is.null(layer.drops) == T) {
            layer.drops3 <- var.f3
        }else{
            layer.drops3 <- c(layer.drops, var.f)
        }

        x3 <- raster::dropLayer(x, which(names(x) %in% layer.drops3))
        x3 <- raster::stack(x3)
        if (raster::nlayers(x3) == 0) {
            MAXENT <- 0
            MAXLIKE <- 0
        }

        if(var.f %in% factors) {
            EARTH <- 0
            GLMNET <- 0
        }

    tests <- ensemble.calibrate.models(x=x3, 
        p=p1, a=a1, an=an, excludep=excludep, 
        k=k, pt=p2, at=a2, SSB.reduce=SSB.reduce, CIRCLES.d=CIRCLES.d,
        TrainData=TrainData2, TestData=TestData2, 
        PLOTS=FALSE, evaluations.keep=T,  
        VIF=VIF, COR=COR,
        formulae.defaults=T, maxit=maxit,
        ENSEMBLE.tune=ENSEMBLE.tune,
        ENSEMBLE.best=ENSEMBLE.best, ENSEMBLE.min=ENSEMBLE.min,
        ENSEMBLE.exponent=ENSEMBLE.exponent,
        input.weights=input.weights,
        MAXENT=MAXENT, MAXNET=MAXNET, MAXLIKE=MAXLIKE, GBM=GBM, GBMSTEP=GBMSTEP, RF=RF, CF=CF,
        GLM=GLM, GLMSTEP=GLMSTEP, GAM=GAM, GAMSTEP=GAMSTEP, MGCV=MGCV, MGCVFIX=MGCVFIX, 
        EARTH=EARTH, RPART=RPART, NNET=NNET, FDA=FDA, SVM=SVM, SVME=SVME, GLMNET=GLMNET,
        BIOCLIM.O=BIOCLIM.O, BIOCLIM=BIOCLIM, DOMAIN=DOMAIN, 
        MAHAL=MAHAL, MAHAL01=MAHAL01,
        PROBIT=PROBIT,
        Yweights=Yweights, 
        factors=factors2, dummy.vars=dummy.vars2,
        MAXENT.a=MAXENT.a, MAXENT.path=MAXENT.path,
        MAXNET.classes=MAXNET.classes, MAXNET.clamp=MAXNET.clamp, MAXNET.type=MAXNET.type,
        MAXLIKE.formula=NULL, MAXLIKE.method=MAXLIKE.method,
        GBM.formula=NULL, GBM.n.trees=GBM.n.trees,
        GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
        GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction,
        GBMSTEP.step.size=GBMSTEP.step.size,
        RF.formula=NULL, RF.ntree=RF.ntree, 
        CF.formula=NULL, CF.ntree=CF.ntree, 
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
        GLMNET.nlambda=GLMNET.nlambda,
        BIOCLIM.O.fraction=BIOCLIM.O.fraction,
        MAHAL.shape=MAHAL.shape)

    if(is.null(tests$evaluations$MAXENT.C)==F) {output1.C["MAXENT",i+1] <- tests$evaluations$MAXENT.C@auc}
    if(is.null(tests$evaluations$MAXENT.T)==F) {output1.T["MAXENT",i+1] <- tests$evaluations$MAXENT.T@auc}
    if(sum(tests$evaluations$TrainData$MAXENT) > 0) {output1.LLC["MAXENT",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MAXENT)}
    if(sum(tests$evaluations$TestData$MAXENT) > 0) {output1.LLT["MAXENT",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MAXENT)}

    if(is.null(tests$evaluations$MAXNET.C)==F) {output1.C["MAXNET",i+1] <- tests$evaluations$MAXNET.C@auc}
    if(is.null(tests$evaluations$MAXNET.T)==F) {output1.T["MAXNET",i+1] <- tests$evaluations$MAXNET.T@auc}
    if(sum(tests$evaluations$TrainData$MAXNET) > 0) {output1.LLC["MAXNET",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MAXNET)}
    if(sum(tests$evaluations$TestData$MAXNET) > 0) {output1.LLT["MAXNET",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MAXNET)}

    if(is.null(tests$evaluations$MAXLIKE.C)==F) {output1.C["MAXLIKE",i+1] <- tests$evaluations$MAXLIKE.C@auc}
    if(is.null(tests$evaluations$MAXLIKE.T)==F) {output1.T["MAXLIKE",i+1] <- tests$evaluations$MAXLIKE.T@auc}
    if(sum(tests$evaluations$TrainData$MAXLIKE) > 0) {output1.LLC["MAXLIKE",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MAXLIKE)}
    if(sum(tests$evaluations$TestData$MAXLIKE) > 0) {output1.LLT["MAXLIKE",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MAXLIKE)}

    if(is.null(tests$evaluations$GBM.C)==F) {output1.C["GBM",i+1] <- tests$evaluations$GBM.C@auc} 
    if(is.null(tests$evaluations$GBM.T)==F) {output1.T["GBM",i+1] <- tests$evaluations$GBM.T@auc} 
    if(sum(tests$evaluations$TrainData$GBM) > 0) {output1.LLC["GBM",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GBM)}
    if(sum(tests$evaluations$TestData$GBM) > 0) {output1.LLT["GBM",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GBM)}

    if(is.null(tests$evaluations$GBMSTEP.C)==F) {output1.C["GBMSTEP",i+1] <- tests$evaluations$GBMSTEP.C@auc} 
    if(is.null(tests$evaluations$GBMSTEP.T)==F) {output1.T["GBMSTEP",i+1] <- tests$evaluations$GBMSTEP.T@auc} 
    if(sum(tests$evaluations$TrainData$GBMSTEP) > 0) {output1.LLC["GBMSTEP",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GBMSTEP)}
    if(sum(tests$evaluations$TestData$GBMSTEP) > 0) {output1.LLT["GBMSTEP",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GBMSTEP)}

    if(is.null(tests$evaluations$RF.C)==F) {output1.C["RF",i+1] <- tests$evaluations$RF.C@auc}
    if(is.null(tests$evaluations$RF.T)==F) {output1.T["RF",i+1] <- tests$evaluations$RF.T@auc}
    if(sum(tests$evaluations$TrainData$RF) > 0) {output1.LLC["RF",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$RF)}
    if(sum(tests$evaluations$TestData$RF) > 0) {output1.LLT["RF",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$RF)}

    if(is.null(tests$evaluations$CF.C)==F) {output1.C["CF",i+1] <- tests$evaluations$CF.C@auc}
    if(is.null(tests$evaluations$CF.T)==F) {output1.T["CF",i+1] <- tests$evaluations$CF.T@auc}
    if(sum(tests$evaluations$TrainData$CF) > 0) {output1.LLC["CF",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$CF)}
    if(sum(tests$evaluations$TestData$CF) > 0) {output1.LLT["CF",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$CF)}

    if(is.null(tests$evaluations$GLM.C)==F) {output1.C["GLM",i+1] <- tests$evaluations$GLM.C@auc} 
    if(is.null(tests$evaluations$GLM.T)==F) {output1.T["GLM",i+1] <- tests$evaluations$GLM.T@auc} 
    if(sum(tests$evaluations$TrainData$GLM) > 0) {output1.LLC["GLM",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GLM)}
    if(sum(tests$evaluations$TestData$GLM) > 0) {output1.LLT["GLM",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GLM)}

    if(is.null(tests$evaluations$GLMS.C)==F) {output1.C["GLMSTEP",i+1] <- tests$evaluations$GLMS.C@auc}
    if(is.null(tests$evaluations$GLMS.T)==F) {output1.T["GLMSTEP",i+1] <- tests$evaluations$GLMS.T@auc}
    if(sum(tests$evaluations$TrainData$GLMSTEP) > 0) {output1.LLC["GLMSTEP",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GLMSTEP)}
    if(sum(tests$evaluations$TestData$GLMSTEP) > 0) {output1.LLT["GLMSTEP",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GLMSTEP)}

    if(is.null(tests$evaluations$GAM.C)==F) {output1.C["GAM",i+1] <- tests$evaluations$GAM.C@auc} 
    if(is.null(tests$evaluations$GAM.T)==F) {output1.T["GAM",i+1] <- tests$evaluations$GAM.T@auc} 
    if(sum(tests$evaluations$TrainData$GAM) > 0) {output1.LLC["GAM",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GAM)}
    if(sum(tests$evaluations$TestData$GAM) > 0) {output1.LLT["GAM",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GAM)}

    if(is.null(tests$evaluations$GAMS.C)==F) {output1.C["GAMSTEP",i+1] <- tests$evaluations$GAMS.C@auc}
    if(is.null(tests$evaluations$GAMS.T)==F) {output1.T["GAMSTEP",i+1] <- tests$evaluations$GAMS.T@auc}
    if(sum(tests$evaluations$TrainData$GAMSTEP) > 0) {output1.LLC["GAMSTEP",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GAMSTEP)}
    if(sum(tests$evaluations$TestData$GAMSTEP) > 0) {output1.LLT["GAMSTEP",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GAMSTEP)}

    if(is.null(tests$evaluations$MGCV.C)==F) {output1.C["MGCV",i+1] <- tests$evaluations$MGCV.C@auc} 
    if(is.null(tests$evaluations$MGCV.T)==F) {output1.T["MGCV",i+1] <- tests$evaluations$MGCV.T@auc} 
    if(sum(tests$evaluations$TrainData$MGCV) > 0) {output1.LLC["MGCV",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MGCV)}
    if(sum(tests$evaluations$TestData$MGCV) > 0) {output1.LLT["MGCV",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MGCV)}

    if(is.null(tests$evaluations$MGCVF.C)==F) {output1.C["MGCVFIX",i+1] <- tests$evaluations$MGCVF.C@auc} 
    if(is.null(tests$evaluations$MGCVF.T)==F) {output1.T["MGCVFIX",i+1] <- tests$evaluations$MGCVF.T@auc} 
    if(sum(tests$evaluations$TrainData$MGCVFIX) > 0) {output1.LLC["MGCVFIX",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MGCVFIX)}
    if(sum(tests$evaluations$TestData$MGCVFIX) > 0) {output1.LLT["MGCVFIX",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MGCVFIX)}

    if(is.null(tests$evaluations$EARTH.C)==F) {output1.C["EARTH",i+1] <- tests$evaluations$EARTH.C@auc} 
    if(is.null(tests$evaluations$EARTH.T)==F) {output1.T["EARTH",i+1] <- tests$evaluations$EARTH.T@auc} 
    if(sum(tests$evaluations$TrainData$EARTH) > 0) {output1.LLC["EARTH",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$EARTH)}
    if(sum(tests$evaluations$TestData$EARTH) > 0) {output1.LLT["EARTH",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$EARTH)}

    if(is.null(tests$evaluations$RPART.C)==F) {output1.C["RPART",i+1] <- tests$evaluations$RPART.C@auc}
    if(is.null(tests$evaluations$RPART.T)==F) {output1.T["RPART",i+1] <- tests$evaluations$RPART.T@auc}
    if(sum(tests$evaluations$TrainData$RPART) > 0) {output1.LLC["RPART",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$RPART)}
    if(sum(tests$evaluations$TestData$RPART) > 0) {output1.LLT["RPART",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$RPART)}

    if(is.null(tests$evaluations$NNET.C)==F) {output1.C["NNET",i+1] <- tests$evaluations$NNET.C@auc} 
    if(is.null(tests$evaluations$NNET.T)==F) {output1.T["NNET",i+1] <- tests$evaluations$NNET.T@auc} 
    if(sum(tests$evaluations$TrainData$NNET) > 0) {output1.LLC["NNET",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$NNET)}
    if(sum(tests$evaluations$TestData$NNET) > 0) {output1.LLT["NNET",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$NNET)}

    if(is.null(tests$evaluations$FDA.C)==F) {output1.C["FDA",i+1] <- tests$evaluations$FDA.C@auc}
    if(is.null(tests$evaluations$FDA.T)==F) {output1.T["FDA",i+1] <- tests$evaluations$FDA.T@auc}
    if(sum(tests$evaluations$TrainData$FDA) > 0) {output1.LLC["FDA",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$FDA)}
    if(sum(tests$evaluations$TestData$FDA) > 0) {output1.LLT["FDA",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$FDA)}

    if(is.null(tests$evaluations$SVM.C)==F) {output1.C["SVM",i+1] <- tests$evaluations$SVM.C@auc}
    if(is.null(tests$evaluations$SVM.T)==F) {output1.T["SVM",i+1] <- tests$evaluations$SVM.T@auc}
    if(sum(tests$evaluations$TrainData$SVM) > 0) {output1.LLC["SVM",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$SVM)}
    if(sum(tests$evaluations$TestData$SVM) > 0) {output1.LLT["SVM",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$SVM)}

    if(is.null(tests$evaluations$SVME.C)==F) {output1.C["SVME",i+1] <- tests$evaluations$SVME.C@auc}
    if(is.null(tests$evaluations$SVME.T)==F) {output1.T["SVME",i+1] <- tests$evaluations$SVME.T@auc}
    if(sum(tests$evaluations$TrainData$SVME) > 0) {output1.LLC["SVME",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$SVME)}
    if(sum(tests$evaluations$TestData$SVME) > 0) {output1.LLT["SVME",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$SVME)}

    if(is.null(tests$evaluations$GLMNET.C)==F) {output1.C["GLMNET",i+1] <- tests$evaluations$GLMNET.C@auc}
    if(is.null(tests$evaluations$GLMNET.T)==F) {output1.T["GLMNET",i+1] <- tests$evaluations$GLMNET.T@auc}
    if(sum(tests$evaluations$TrainData$GLMNET) > 0) {output1.LLC["GLMNET",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GLMNET)}
    if(sum(tests$evaluations$TestData$GLMNET) > 0) {output1.LLT["GLMNET",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GLMNET)}

    if(is.null(tests$evaluations$BIOCLIM.O.C)==F) {output1.C["BIOCLIM.O",i+1] <- tests$evaluations$BIOCLIM.O.C@auc}
    if(is.null(tests$evaluations$BIOCLIM.O.T)==F) {output1.T["BIOCLIM.O",i+1] <- tests$evaluations$BIOCLIM.O.T@auc}
    if(sum(tests$evaluations$TrainData$BIOCLIM.O) > 0) {output1.LLC["BIOCLIM.O",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$BIOCLIM.O)}
    if(sum(tests$evaluations$TestData$BIOCLIM.O) > 0) {output1.LLT["BIOCLIM.O",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$BIOCLIM.O)}

    if(is.null(tests$evaluations$BIOCLIM.C)==F) {output1.C["BIOCLIM",i+1] <- tests$evaluations$BIOCLIM.C@auc}
    if(is.null(tests$evaluations$BIOCLIM.T)==F) {output1.T["BIOCLIM",i+1] <- tests$evaluations$BIOCLIM.T@auc}
    if(sum(tests$evaluations$TrainData$BIOCLIM) > 0) {output1.LLC["BIOCLIM",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$BIOCLIM)}
    if(sum(tests$evaluations$TestData$BIOCLIM) > 0) {output1.LLT["BIOCLIM",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$BIOCLIM)}

    if(is.null(tests$evaluations$DOMAIN.C)==F) {output1.C["DOMAIN",i+1] <- tests$evaluations$DOMAIN.C@auc}
    if(is.null(tests$evaluations$DOMAIN.T)==F) {output1.T["DOMAIN",i+1] <- tests$evaluations$DOMAIN.T@auc}
    if(sum(tests$evaluations$TrainData$DOMAIN) > 0) {output1.LLC["DOMAIN",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$DOMAIN)}
    if(sum(tests$evaluations$TestData$DOMAIN) > 0) {output1.LLT["DOMAIN",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$DOMAIN)}

    if(is.null(tests$evaluations$MAHAL.C)==F) {output1.C["MAHAL",i+1] <- tests$evaluations$MAHAL.C@auc}
    if(is.null(tests$evaluations$MAHAL.T)==F) {output1.T["MAHAL",i+1] <- tests$evaluations$MAHAL.T@auc}
    if(sum(tests$evaluations$TrainData$MAHAL) > 0) {output1.LLC["MAHAL",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MAHAL)}
    if(sum(tests$evaluations$TestData$MAHAL) > 0) {output1.LLT["MAHAL",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MAHAL)}

    if(is.null(tests$evaluations$MAHAL01.C)==F) {output1.C["MAHAL01",i+1] <- tests$evaluations$MAHAL01.C@auc}
    if(is.null(tests$evaluations$MAHAL01.T)==F) {output1.T["MAHAL01",i+1] <- tests$evaluations$MAHAL01.T@auc}
    if(sum(tests$evaluations$TrainData$MAHAL01) > 0) {output1.LLC["MAHAL01",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MAHAL01)}
    if(sum(tests$evaluations$TestData$MAHAL01) > 0) {output1.LLT["MAHAL01",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MAHAL01)}

    if(is.null(tests$evaluations$ENSEMBLE.C)==F) {output1.C["ENSEMBLE",i+1] <- tests$evaluations$ENSEMBLE.C@auc}
    if(is.null(tests$evaluations$ENSEMBLE.T)==F) {output1.T["ENSEMBLE",i+1] <- tests$evaluations$ENSEMBLE.T@auc}
    if(sum(tests$evaluations$TrainData$ENSEMBLE) > 0) {output1.LLC["ENSEMBLE",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$ENSEMBLE)}
    if(sum(tests$evaluations$TestData$ENSEMBLE) > 0) {output1.LLT["ENSEMBLE",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$ENSEMBLE)}

    }

#
# end of variables alone loop
    }

## Arrange the tables

    output.C <- 100*output.C
    if (difference == T) {
        for (i in 1:nv) {
            output.C[,i+1] <- output.C[,i+1] - output.C[,1]
        }
    }
    output.C <- output.C[order(output.C[,"all_vars"], decreasing=T),]
    cat(paste("\n", "AUC for calibration data (as percentage)", "\n\n", sep = ""))
    if (difference == T) {
        cat(paste("\n", "Results for variables show change from full models", sep = ""))
        cat(paste("\n", "Note that positive differences indicate that the model without the variable", sep = ""))
        cat(paste("\n", "has higher AUC than the model with all the variables",  "\n\n", sep = ""))
    }
    print (output.C)

    output.T <- 100*output.T
    if (difference == T) {
        for (i in 1:nv) {
            output.T[,i+1] <- output.T[,i+1] - output.T[,1]
        }
    }
    output.T <- output.T[order(output.T[,"all_vars"], decreasing=T),]
    cat(paste("\n", "AUC for testing data (as percentage)", "\n\n", sep = ""))
    if (difference == T) {
        cat(paste("\n", "Results for variables show change from full models", sep = ""))
        cat(paste("\n", "Note that positive differences indicate that the model without the variable", sep = ""))
        cat(paste("\n", "has higher AUC than the model with all the variables",  "\n\n", sep = ""))
    }
    print (output.T)
#
# base null model on predictions of prevalence

    cat(paste("\n", "Null deviance for calibration data: ", preval.dev.cal, sep=""))
    cat(paste("\n", "Null deviance for worst possible predictions of calibration data: ", null.dev.cal, sep=""))
    cat(paste("\n", "Residual deviance for calibration data",  "\n", sep = ""))

    percentage.LLC <- output.LLC

    if (difference == T) {
        for (i in 1:nv) {
            output.LLC[,i+1] <- output.LLC[,1] - output.LLC[,i+1]
        }
        cat(paste("\n", "Results for variables show change from full models", sep = ""))
    }
    output.LLC <- output.LLC[order(output.LLC[,"all_vars"], decreasing=F),]
    cat(paste("\n\n"))
    print (output.LLC)

    cat(paste("\n", "Null deviance for testing data: ", preval.dev.test, sep=""))
    cat(paste("\n", "Null deviance for worst possible predictions of testing data: ", null.dev.test, sep=""))
    cat(paste("\n", "Residual deviance for testing data",  "\n", sep = ""))

    percentage.LLT <- output.LLT

    if (difference == T) {
        for (i in 1:nv) {
            output.LLT[,i+1] <-  output.LLT[,1] - output.LLT[,i+1]
        }
        cat(paste("\n", "Results for variables show change from full models", sep = ""))
    }

    output.LLT <- output.LLT[order(output.LLT[,"all_vars"], decreasing=F),]
    cat(paste("\n\n"))
    print (output.LLT)

    for (i in 1:(1+nv)) {
        percentage.LLC[, i] <- percentage.LLC[, i] - as.numeric(null.dev.cal)
        percentage.LLC[, i] <- -100 * percentage.LLC[, i]
        percentage.LLC[, i] <- percentage.LLC[, i] / as.numeric(null.dev.cal)
        percentage.LLC[, i] <- round(percentage.LLC[, i], 2) 

        percentage.LLT[, i] <- percentage.LLT[, i] - as.numeric(null.dev.test)
        percentage.LLT[, i] <- -100 * percentage.LLT[, i]
        percentage.LLT[, i] <- percentage.LLT[, i] / as.numeric(null.dev.test)
        percentage.LLT[, i] <- round(percentage.LLT[, i], 2) 
    }

    if (difference == T) {
        for (i in 1:nv) {
            percentage.LLC[,i+1] <- percentage.LLC[,1] - percentage.LLC[,i+1]
            percentage.LLT[,i+1] <- percentage.LLT[,1] - percentage.LLT[,i+1]
        }
        cat(paste("\n", "Results for variables show change from full models", sep = ""))
    }


    cat(paste("\n", "Percentage explained for calibration data",  "\n\n", sep = ""))
    percentage.LLC <- percentage.LLC[order(percentage.LLC[,"all_vars"], decreasing=F),]
    print(percentage.LLC)

    cat(paste("\n", "Percentage explained for testing data",  "\n\n", sep = ""))
    percentage.LLT <- percentage.LLT[order(percentage.LLT[,"all_vars"], decreasing=F),]
    print (percentage.LLT)

## Models with one variable only

    if (variables.alone == T) {

    output1.C <- 100*output1.C
    if (difference == T) {
        for (i in 1:nv) {
            output1.C[,i+1] <- output1.C[,i+1] - output1.C[,1]
        }
    }
    output1.C <- output1.C[order(output1.C[,"all_vars"], decreasing=T),]
    cat(paste("\n", "AUC for calibration data (as percentage)", "\n\n", sep = ""))
    if (difference == T) {
        cat(paste("\n", "Results for variables show change from full models", sep = ""))
        cat(paste("\n", "Note that positive differences indicate that the model only with the variable", sep = ""))
        cat(paste("\n", "has higher AUC than the model with all the variables",  "\n\n", sep = ""))
    }
    print (output1.C)

    output1.T <- 100*output1.T
    if (difference == T) {
        for (i in 1:nv) {
            output1.T[,i+1] <- output1.T[,i+1] - output1.T[,1]
        }
    }
    output1.T <- output1.T[order(output1.T[,"all_vars"], decreasing=F),]
    cat(paste("\n", "AUC for testing data (as percentage)", "\n\n", sep = ""))
    if (difference == T) {
        cat(paste("\n", "Results for variables show change from full models", sep = ""))
        cat(paste("\n", "Note that positive differences indicate that the model only with the variable", sep = ""))
        cat(paste("\n", "has higher AUC than the model with all the variables",  "\n\n", sep = ""))
    }
    print (output1.T)
#
# base null model on predictions of prevalence

    cat(paste("\n", "Null deviance for calibration data: ", preval.dev.cal, sep=""))
    cat(paste("\n", "Null deviance for worst possible predictions of calibration data: ", null.dev.cal, sep=""))
    cat(paste("\n", "Residual deviance for calibration data",  "\n", sep = ""))

    percentage1.LLC <- output1.LLC

    if (difference == T) {
        for (i in 1:nv) {
            output1.LLC[,i+1] <- output1.LLC[,1] - output1.LLC[,i+1]
        }
        cat(paste("\n", "Results for variables show change from full models", sep = ""))
    }
    output1.LLC <- output1.LLC[order(output1.LLC[,"all_vars"], decreasing=F),]
    cat(paste("\n\n"))
    print (output1.LLC)

    cat(paste("\n", "Null deviance for testing data: ", preval.dev.test, sep=""))
    cat(paste("\n", "Null deviance for worst possible predictions of testing data: ", null.dev.test, sep=""))
    cat(paste("\n", "Residual deviance for testing data",  "\n", sep = ""))

    percentage1.LLT <- output1.LLT

    if (difference == T) {
        for (i in 1:nv) {
            output1.LLT[,i+1] <-  output1.LLT[,1] - output1.LLT[,i+1]
        }
        cat(paste("\n", "Results for variables show change from full models", sep = ""))
    }

    output1.LLT <- output1.LLT[order(output1.LLT[,"all_vars"], decreasing=F),]
    cat(paste("\n\n"))
    print (output1.LLT)

    for (i in 1:(1+nv)) {
        percentage1.LLC[, i] <- percentage1.LLC[, i] - as.numeric(null.dev.cal)
        percentage1.LLC[, i] <- -100 * percentage1.LLC[, i]
        percentage1.LLC[, i] <- percentage1.LLC[, i] / as.numeric(null.dev.cal)
        percentage1.LLC[, i] <- round(percentage1.LLC[, i], 2) 

        percentage1.LLT[, i] <- percentage1.LLT[, i] - as.numeric(null.dev.test)
        percentage1.LLT[, i] <- -100 * percentage1.LLT[, i]
        percentage1.LLT[, i] <- percentage1.LLT[, i] / as.numeric(null.dev.test)
        percentage1.LLT[, i] <- round(percentage1.LLT[, i], 2) 
    }

    if (difference == T) {
        for (i in 1:nv) {
            percentage1.LLC[,i+1] <- percentage1.LLC[,1] - percentage1.LLC[,i+1]
            percentage1.LLT[,i+1] <- percentage1.LLT[,1] - percentage1.LLT[,i+1]
        }
        cat(paste("\n", "Results for variables show change from full models", sep = ""))
    }


    cat(paste("\n", "Percentage explained for calibration data",  "\n\n", sep = ""))
    percentage1.LLC <- percentage1.LLC[order(percentage1.LLC[,"all_vars"], decreasing=F),]
    print(percentage1.LLC)

    cat(paste("\n", "Percentage explained for testing data",  "\n\n", sep = ""))
    percentage1.LLT <- percentage1.LLT[order(percentage1.LLT[,"all_vars"], decreasing=F),]
    print (percentage1.LLT)

#
# end of variables alone loop
    }

#
    if (SINK==T && OLD.SINK==F) {sink(file=NULL, append=T)}
#
    if (variables.alone == F) {
        return(list(AUC.calibration=output.C, AUC.testing=output.T, 
            residual.deviance.calibration=output.LLC, residual.deviance.testing=output.LLT, 
            percentage.deviance.calibration=percentage.LLC, percentage.deviance.testing=percentage.LLT,
            call=match.call() ))
    }else{
        return(list(AUC.calibration=output.C, AUC.testing=output.T, 
            residual.deviance.calibration=output.LLC, residual.deviance.testing=output.LLT, 
            percentage.deviance.calibration=percentage.LLC, percentage.deviance.testing=percentage.LLT,
            AUC.single.calibration=output1.C, AUC.single.testing=output1.T, 
            residual.single.deviance.calibration=output1.LLC, residual.single.deviance.testing=output1.LLT, 
            percentage.single.deviance.calibration=percentage1.LLC, percentage.single.deviance.testing=percentage1.LLT,
            call=match.call() ))
    }    

}

