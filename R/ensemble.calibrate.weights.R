`ensemble.calibrate.weights` <- function(
    x=NULL, p=NULL, TrainTestData=NULL,
    a=NULL, an=1000, 
    get.block=FALSE, block.default=TRUE, get.subblocks=FALSE, 
    SSB.reduce=FALSE, CIRCLES.d=100000,
    excludep=FALSE, target.groups=FALSE, k=4, 
    VIF=FALSE, COR=FALSE,
    SINK=FALSE, PLOTS=FALSE, CATCH.OFF=FALSE,
    data.keep=FALSE,
    species.name = "Species001",
    threshold.method="spec_sens", threshold.sensitivity=0.9, threshold.PresenceAbsence=FALSE,
    ENSEMBLE.tune=FALSE, 
    ENSEMBLE.best=0, ENSEMBLE.min=0.7, ENSEMBLE.exponent=1, ENSEMBLE.weight.min=0.05,
    input.weights=NULL,
    MAXENT=1, MAXNET=1, MAXLIKE=1, GBM=1, GBMSTEP=1, RF=1, CF=1,
    GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=1, MGCV=1, MGCVFIX=0,
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, GLMNET=1,
    BIOCLIM.O=0, BIOCLIM=1, DOMAIN=1, MAHAL=1, MAHAL01=1,
    PROBIT=FALSE,
    Yweights="BIOMOD", 
    layer.drops=NULL, factors=NULL, dummy.vars=NULL,
    formulae.defaults=TRUE, maxit=100,
    MAXENT.a=NULL, MAXENT.an=10000, 
    MAXENT.path=paste(getwd(), "/models/maxent_", species.name,  sep=""), 
    MAXNET.classes="default", MAXNET.clamp=FALSE, MAXNET.type="cloglog",
    MAXLIKE.formula=NULL, MAXLIKE.method="BFGS",
    GBM.formula=NULL, GBM.n.trees=2001, 
    GBMSTEP.gbm.x=2:(length(var.names)+1), GBMSTEP.tree.complexity=5, GBMSTEP.learning.rate=0.005, 
    GBMSTEP.bag.fraction=0.5, GBMSTEP.step.size=100, 
    RF.formula=NULL, RF.ntree=751, RF.mtry=floor(sqrt(length(var.names))), 
    CF.formula=NULL, CF.ntree=751, CF.mtry=floor(sqrt(length(var.names))), 
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

    if (is.list(k) == F) {
        k <- as.integer(k)
        k.listed <- FALSE
        if (k < 2) {
            cat(paste("\n", "NOTE: parameter k was set to be smaller than 2", sep = ""))
            cat(paste("\n", "default value of 4 therefore set for parameter k", "\n", sep = ""))
            k <- 4
        }
    }else{
        k.listed <- TRUE
        k.list <- k
        k <- max(k.list$groupa)
    }

    x.was.terra <- FALSE
    
# new option using data.frame
    if (is.null(TrainTestData == TRUE)) {
        if(inherits(x, "RasterStack") == FALSE && inherits(x, "SpatRaster") == FALSE) {stop("x is not a RasterStack object")}
# use the raster format for preparing and checking the data
        if(inherits(x, "SpatRaster")) {
            x <- raster::stack(x)
            x.was.terra <- TRUE
        }
        if(is.null(p) == T) {stop("presence locations are missing (parameter p)")}
        if(is.null(p) == F) {
            p <- data.frame(p)
            names(p) <- c("x", "y")
        }
        if(is.null(a) == F) {
            a <- data.frame(a)
            names(a) <- c("x", "y")
        }
        if(is.null(MAXENT.a) == F) {
            MAXENT.a <- data.frame(MAXENT.a)
            names(MAXENT.a) <- c("x", "y")
        }
        if (is.null(layer.drops) == F) {
            layer.drops <- as.character(layer.drops)
            if (inherits(x, "RasterStack")) {
                x <- raster::dropLayer(x, which(names(x) %in% layer.drops))
                x <- raster::stack(x)
            }
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
    } # no data.frame

    if (is.null(TrainTestData) == F) {
        TrainTestData <- data.frame(TrainTestData)
        if (names(TrainTestData)[1] !="pb") {stop("first column for TrainTestData should be 'pb' containing presence (1) and absence (0) data")}
        if (inherits(x, "RasterStack")) { 
            if (raster::nlayers(x) != (ncol(TrainTestData)-1)) {
                cat(paste("\n", "WARNING: different number of explanatory variables in rasterStack and TrainTestData", sep = ""))
            }
        }
    }
#
    output.rownames <- c("MAXENT", "MAXNET", "MAXLIKE", "GBM", "GBMSTEP", "RF", "CF", 
        "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "GLMNET", 
        "BIOCLIM.O", "BIOCLIM", "DOMAIN", "MAHAL", "MAHAL01", "ENSEMBLE")
#
#    if(length(ENSEMBLE.exponent) > 1 || length(ENSEMBLE.best) > 1 || length(ENSEMBLE.min) > 1) {ENSEMBLE.tune <- TRUE}
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
    TrainData.all <- vector("list", k)
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
# also obtain p and a, possibly via target group sampling
# only do this when no data.frame

    if (is.null(TrainTestData) == TRUE) {

    tests <- ensemble.calibrate.models(x=x, 
        p=p, a=a, an=an, pt=NULL, at=NULL, excludep=excludep, target.groups=target.groups,
        k=0, TrainData=NULL, 
        VIF=F, COR=F,
        PLOTS=PLOTS, evaluations.keep=T, models.keep=F,
        ENSEMBLE.tune=F,
        ENSEMBLE.exponent=1, ENSEMBLE.best=1, ENSEMBLE.min=0.7,
        MAXENT=0, MAXNET=0, MAXLIKE=0, GBM=0, GBMSTEP=0, RF=0, CF=0, GLM=0, GLMSTEP=0, 
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

# new option of k.list surpasses other options of assigning observations to k-folds

    if (k.listed == F) {
        if (get.subblocks == T) {get.block <- T}

        if (get.block == F) {
            groupp <- dismo::kfold(p.all, k=k)
            groupa <- dismo::kfold(a.all, k=k)
        }else{
            p2.all <- p.all
            a2.all <- a.all
            if (block.default == F) {
                p2.all[, 1] <- p.all[, 2]
                p2.all[, 2] <- p.all[, 1]
                a2.all[, 1] <- a.all[, 2]
                a2.all[, 2] <- a.all[, 1]
                cat(paste("\n", "non-default ENMeval::get.block with first split along longitude", "\n\n", sep = ""))
            }else{
                cat(paste("\n", "default ENMeval::get.block with first split along latitude", "\n\n", sep = ""))
            }
            blocks <- ENMeval::get.block(occ=p2.all, bg.coords=a2.all)
            groupp <- blocks$occ.grp
            groupa <- blocks$bg.grp
            k <- 4

    # 'subblocking' whereby get.block is applied to each previously determined block
            if (get.subblocks == T) {
                occ.old <- groupp
                backg.old <- groupa
                for (i in 1:4) {
                    occ.i <- p2.all[occ.old == i, ]
                    backg.i <- a2.all[backg.old == i, ]
                    block2 <- ENMeval::get.block(occ=occ.i, bg.coords=backg.i)
                    occ.new <- block2$occ.grp
                    backg.new <- block2$bg.grp
                    groupp[occ.old == i] <- occ.new
                    groupa[backg.old == i] <- backg.new
                }
            }
        }

    }else{
        groupp <- k.list$groupp
        groupa <- k.list$groupa
        if (length(groupp) != nrow(p.all)) {cat(paste("WARNING: groupp length (", length(groupp), ") different from number of presence observations (", nrow(p.all), ")", "\n", sep = ""))}
        if (length(groupa) != nrow(a.all)) {cat(paste("WARNING: groupa length (", length(groupa), ") different from number of background observations (", nrow(a.all), ")", "\n", sep = ""))}
    }

    } # no data.frame

    if (is.null(TrainTestData) == FALSE) {

    TrainTest.p <- TrainTestData[TrainTestData[, "pb"] == 1, ]
    TrainTest.a <- TrainTestData[TrainTestData[, "pb"] == 0, ]
    var.names <- names(TrainTestData)
    var.names <- var.names[which(var.names != "pb")]
    p.all <- NULL
    a.all <- NULL

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

    if (k.listed == F) {
        groupp <- dismo::kfold(TrainTest.p, k=k)
        groupa <- dismo::kfold(TrainTest.a, k=k)
    }else{
        groupp <- k.list$groupp
        groupa <- k.list$groupa
        if (length(groupp) != nrow(TrainTest.p)) {cat(paste("WARNING: groupp length (", length(groupp), ") different from number of presence observations (", nrow(TrainTest.p), ")", "\n", sep = ""))}
        if (length(groupa) != nrow(TrainTest.a)) {cat(paste("WARNING: groupa length (", length(groupa), ") different from number of background observations (", nrow(TrainTest.a), ")", "\n", sep = ""))}
    }
    } # data.frame

# Start cross-validations
    
    for (i in 1:k){
        cat(paste(species.name, " ", k, "-FOLD CROSS-VALIDATION RUN: ", i, "\n", sep = ""))

    if (is.null(TrainTestData) == TRUE) {

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
                MAXENT=MAXENT, MAXNET=MAXNET, MAXLIKE=MAXLIKE, GBM=GBM, GBMSTEP=GBMSTEP, RF=RF, CF=CF, 
                GLM=GLM, GLMSTEP=GLMSTEP, GAM=GAM, GAMSTEP=GAMSTEP, MGCV=MGCV, MGCVFIX=MGCVFIX, 
                EARTH=EARTH, RPART=RPART, NNET=NNET, FDA=FDA, SVM=SVM, SVME=SVME, GLMNET=GLMNET,
                BIOCLIM.O=BIOCLIM.O, BIOCLIM=BIOCLIM, DOMAIN=DOMAIN, MAHAL=MAHAL, MAHAL01=MAHAL01,
                PROBIT=PROBIT,  
                Yweights=Yweights, 
                factors=factors2, dummy.vars=dummy.vars2,
                maxit=maxit,
                MAXENT.a=MAXENT.a, MAXENT.path=MAXENT.path,
                MAXNET.clamp=MAXNET.clamp, MAXNET.type=MAXNET.type,
                MAXLIKE.formula=MAXLIKE.formula, MAXLIKE.method=MAXLIKE.method,
                GBM.formula=GBM.formula, GBM.n.trees=GBM.n.trees, 
                GBMSTEP.gbm.x=GBMSTEP.gbm.x, GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
                GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction,
                GBMSTEP.step.size=GBMSTEP.step.size, 
                RF.formula=RF.formula, RF.ntree=RF.ntree, RF.mtry=RF.mtry, 
                CF.formula=CF.formula, CF.ntree=CF.ntree, CF.mtry=CF.mtry, 
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
        
        } # no data.frame
        
    if (is.null(TrainTestData) == FALSE) {

            TrainTest.p1 <- TrainTest.p[groupp != i,]
            TrainTest.p2 <- TrainTest.p[groupp == i,]

            TrainTest.a1 <- TrainTest.a[groupa != i,]
            TrainTest.a2 <- TrainTest.a[groupa == i,]
                        
            TrainTest.train <- rbind(TrainTest.p1, TrainTest.a1)
            TrainTest.test <- rbind(TrainTest.p2, TrainTest.a2)            

            tests <- ensemble.calibrate.models(x=NULL, 
                TrainData=TrainTest.train, TestData=TrainTest.test,
                p=NULL, a=NULL, pt=NULL, at=NULL, SSB.reduce=SSB.reduce, CIRCLES.d=CIRCLES.d,
                VIF=VIF, COR=COR,
                threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence,
                PLOTS=PLOTS, CATCH.OFF=CATCH.OFF,
                evaluations.keep=T, models.keep=F,
                ENSEMBLE.tune=ENSEMBLE.tune,
                ENSEMBLE.best=ENSEMBLE.best, ENSEMBLE.min=ENSEMBLE.min, ENSEMBLE.exponent=ENSEMBLE.exponent, ENSEMBLE.weight.min=ENSEMBLE.weight.min,
                MAXENT=MAXENT, MAXNET=MAXNET, MAXLIKE=MAXLIKE, GBM=GBM, GBMSTEP=GBMSTEP, RF=RF, CF=CF, 
                GLM=GLM, GLMSTEP=GLMSTEP, GAM=GAM, GAMSTEP=GAMSTEP, MGCV=MGCV, MGCVFIX=MGCVFIX, 
                EARTH=EARTH, RPART=RPART, NNET=NNET, FDA=FDA, SVM=SVM, SVME=SVME, GLMNET=GLMNET,
                BIOCLIM.O=BIOCLIM.O, BIOCLIM=BIOCLIM, DOMAIN=DOMAIN, MAHAL=MAHAL, MAHAL01=MAHAL01,
                PROBIT=PROBIT,  
                Yweights=Yweights, 
                factors=factors2, dummy.vars=dummy.vars2,
                maxit=maxit,
                MAXENT.a=MAXENT.a, MAXENT.path=MAXENT.path,
                MAXNET.clamp=MAXNET.clamp, MAXNET.type=MAXNET.type,
                MAXLIKE.formula=MAXLIKE.formula, MAXLIKE.method=MAXLIKE.method,
                GBM.formula=GBM.formula, GBM.n.trees=GBM.n.trees, 
                GBMSTEP.gbm.x=GBMSTEP.gbm.x, GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
                GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction,
                GBMSTEP.step.size=GBMSTEP.step.size, 
                RF.formula=RF.formula, RF.ntree=RF.ntree, RF.mtry=RF.mtry, 
                CF.formula=CF.formula, CF.ntree=CF.ntree, CF.mtry=CF.mtry, 
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
                
        } # data.frame        

        dummy.vars.noDOMAIN <- tests$evaluations$dummy.vars.noDOMAIN

        if(is.null(tests$evaluations$MAXENT.T)==F) {output["MAXENT",i] <- tests$evaluations$MAXENT.T@auc}
        if(is.null(tests$evaluations$MAXNET.T)==F) {output["MAXNET",i] <- tests$evaluations$MAXNET.T@auc}
        if(is.null(tests$evaluations$MAXLIKE.T)==F) {output["MAXLIKE",i] <- tests$evaluations$MAXLIKE.T@auc}
        if(is.null(tests$evaluations$GBM.T)==F) {output["GBM",i] <- tests$evaluations$GBM.T@auc} 
        if(is.null(tests$evaluations$GBMSTEP.T)==F) {output["GBMSTEP",i] <- tests$evaluations$GBMSTEP.T@auc} 
        if(is.null(tests$evaluations$RF.T)==F) {output["RF",i] <- tests$evaluations$RF.T@auc}
        if(is.null(tests$evaluations$CF.T)==F) {output["CF",i] <- tests$evaluations$CF.T@auc}
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
            output["MAXNET",k+1+i] <- tests$evaluations$STRATEGY.weights["MAXNET"]
            output["MAXLIKE",k+1+i] <- tests$evaluations$STRATEGY.weights["MAXLIKE"]
            output["GBM",k+1+i] <- tests$evaluations$STRATEGY.weights["GBM"]
            output["GBMSTEP",k+1+i] <- tests$evaluations$STRATEGY.weights["GBMSTEP"]
            output["RF",k+1+i] <- tests$evaluations$STRATEGY.weights["RF"]
            output["CF",k+1+i] <- tests$evaluations$STRATEGY.weights["CF"]
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

        TrainData.all[[i]] <- tests$evaluations$TrainData
        TestData.all[[i]] <- tests$evaluations$TestData

        eval.tablei <- tests$eval.table
        eval.tablei$ALGO <- rownames(eval.tablei)
        eval.tablei$k <- rep(i, nrow(eval.tablei))
        if (i==1) {
            eval.table.all <- eval.tablei
        }else{
            eval.table.all <- rbind(eval.table.all, eval.tablei)
        }
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
    output.weightsT <- output[, "MEAN.T"]
    output.weights <- output.weights[names(output.weights) != "ENSEMBLE"]
    output.weightsT <- output.weightsT[names(output.weightsT) != "ENSEMBLE"]
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

    if (ENSEMBLE.tune == F) {
# do not modify if ENSEMBLE.tune is not required, simply use same
#        cat(paste("\n\n", "parameters for next weighting: ENSEMBLE.min=", ENSEMBLE.min1, ", ENSEMBLE.best=", ENSEMBLE.best1, " and ENSEMBLE.exponent=", ENSEMBLE.exponent1, "\n\n", sep = ""))
#        output.weights <- ensemble.weights(weights=output.weights, exponent=ENSEMBLE.exponent1, best=ENSEMBLE.best1, min.weight=ENSEMBLE.min1)
#        print(output.weights)
         output.weights <- tests$evaluations$ensemble.weights
        cat(paste("\n", "final = original weights for ensemble forecasting", "\n", sep = ""))
        print(output.weights)
    }else{
        cat(paste("\n", "input weights for ensemble modelling based on MEAN column", "\n",  sep = ""))
        print(output.weights)
# if possible, select best models (models with highest weights)
        if (ENSEMBLE.best1 > 0) {
            cat(paste("\n\n", "parameters for next weighting: ENSEMBLE.min=0, ENSEMBLE.best=", ENSEMBLE.best1, " and ENSEMBLE.exponent=1", "\n\n", sep = ""))
            output.weights <- ensemble.weights(weights=output.weights, exponent=1, best=ENSEMBLE.best1, min.weight=0)
            print(output.weights)
        }
        output.weightsT <- ensemble.weights(weights=output.weightsT, exponent=ENSEMBLE.exponent1, best=ENSEMBLE.best1, min.weight=ENSEMBLE.min1)
# remove models with low input weights (mainly to reduce the number of models for final calibrations and mapping)
        cat(paste("\n", "Minimum input weight is ", ENSEMBLE.weight.min, "\n", sep=""))
        output.weights2 <- output.weights
        output.weights2T <- output.weightsT
        while(min(output.weights2) < ENSEMBLE.weight.min) {
            output.weights2 <- output.weights2[-which.min(output.weights2)]
            output.weights2 <- ensemble.weights(weights=output.weights2, exponent=1, best=0, min.weight=0)
        }
        while(min(output.weights2T) < ENSEMBLE.weight.min) {
            output.weights2T <- output.weights2T[-which.min(output.weights2T)]
            output.weights2T <- ensemble.weights(weights=output.weights2T, exponent=1, best=0, min.weight=0)
        }
        output.weights[] <- 0
        output.weightsT[] <- 0
        for (i in 1:length(output.weights2)) {output.weights[which(names(output.weights) == names(output.weights2)[i])] <- output.weights2[i]}
        for (i in 1:length(output.weights2T)) {output.weightsT[which(names(output.weightsT) == names(output.weights2T)[i])] <- output.weights2T[i]}
        cat(paste("\n", "final suggested weights for ensemble forecasting", "\n", sep = ""))
        print(output.weights)
    }

# test with suggested final weights
# was no need to repeat for no-tuning since here we already have those results, but possibly good to recheck
    output2 <- numeric(length=k+1)
    names(output2)[1:k] <- paste("T_", c(1:k), sep="")
    names(output2)[k+1] <- c("MEAN.T")
    for (i in 1:k) {
        TrainData <- TrainData.all[[i]]
        TrainData[,"ENSEMBLE"] <- output.weights["MAXENT"]*TrainData[,"MAXENT"] + output.weights["MAXNET"]*TrainData[,"MAXNET"] +
            output.weights["MAXLIKE"]*TrainData[,"MAXLIKE"] + output.weights["GBM"]*TrainData[,"GBM"] +
            output.weights["GBMSTEP"]*TrainData[,"GBMSTEP"] + output.weights["RF"]*TrainData[,"RF"] + output.weights["CF"]*TrainData[,"CF"] + output.weights["GLM"]*TrainData[,"GLM"] +
            output.weights["GLMSTEP"]*TrainData[,"GLMSTEP"] + output.weights["GAM"]*TrainData[,"GAM"] + output.weights["GAMSTEP"]*TrainData[,"GAMSTEP"] +
            output.weights["MGCV"]*TrainData[,"MGCV"] + output.weights["MGCVFIX"]*TrainData[,"MGCVFIX"] + output.weights["EARTH"]*TrainData[,"EARTH"] +
            output.weights["RPART"]*TrainData[,"RPART"] + output.weights["NNET"]*TrainData[,"NNET"] + output.weights["FDA"]*TrainData[,"FDA"] +
            output.weights["SVM"]*TrainData[,"SVM"] + output.weights["SVME"]*TrainData[,"SVME"] + output.weights["GLMNET"]*TrainData[,"GLMNET"] + 
            output.weights["BIOCLIM.O"]*TrainData[,"BIOCLIM.O"] + output.weights["BIOCLIM"]*TrainData[,"BIOCLIM"] +
            output.weights["DOMAIN"]*TrainData[,"DOMAIN"] + output.weights["MAHAL"]*TrainData[,"MAHAL"]+ output.weights["MAHAL01"]*TrainData[,"MAHAL01"]
        TestData <- TestData.all[[i]]
        TestData[,"ENSEMBLE"] <- output.weights["MAXENT"]*TestData[,"MAXENT"] + output.weights["MAXNET"]*TestData[,"MAXNET"] +
            output.weights["MAXLIKE"]*TestData[,"MAXLIKE"] + output.weights["GBM"]*TestData[,"GBM"] +
            output.weights["GBMSTEP"]*TestData[,"GBMSTEP"] + output.weights["RF"]*TestData[,"RF"] + output.weights["CF"]*TestData[,"CF"] + output.weights["GLM"]*TestData[,"GLM"] +
            output.weights["GLMSTEP"]*TestData[,"GLMSTEP"] + output.weights["GAM"]*TestData[,"GAM"] + output.weights["GAMSTEP"]*TestData[,"GAMSTEP"] +
            output.weights["MGCV"]*TestData[,"MGCV"] + output.weights["MGCVFIX"]*TestData[,"MGCVFIX"] + output.weights["EARTH"]*TestData[,"EARTH"] +
            output.weights["RPART"]*TestData[,"RPART"] + output.weights["NNET"]*TestData[,"NNET"] + output.weights["FDA"]*TestData[,"FDA"] +
            output.weights["SVM"]*TestData[,"SVM"] + output.weights["SVME"]*TestData[,"SVME"] + output.weights["GLMNET"]*TestData[,"GLMNET"] + 
            output.weights["BIOCLIM.O"]*TestData[,"BIOCLIM.O"] + output.weights["BIOCLIM"]*TestData[,"BIOCLIM"] +
            output.weights["DOMAIN"]*TestData[,"DOMAIN"] + output.weights["MAHAL"]*TestData[,"MAHAL"]+ output.weights["MAHAL01"]*TestData[,"MAHAL01"]
            eval1 <- evalT <- NULL
        TrainPres <- as.numeric(TrainData[TrainData[,"pb"]==1, "ENSEMBLE"])
        TrainAbs <- as.numeric(TrainData[TrainData[,"pb"]==0, "ENSEMBLE"])
        evalT <- dismo::evaluate(p=TrainPres, a=TrainAbs)
        TestPres <- as.numeric(TestData[TestData[,"pb"]==1, "ENSEMBLE"])
        TestAbs <- as.numeric(TestData[TestData[,"pb"]==0, "ENSEMBLE"])
        eval1 <- dismo::evaluate(p=TestPres, a=TestAbs)
        output2[i] <- eval1@auc
        cat(paste("\n", "Results with final weights for ensemble.evaluate for k = ", i, "\n\n", sep = ""))
        eval3 <- ensemble.evaluate(eval=eval1, eval.train=evalT)
        print(eval3)
        if (i==1) {
            eval.table.final <- data.frame(t(eval3))
        }else{
            eval.table.final <- data.frame(rbind(eval.table.final, t(eval3)))
        }
    }
    output2[k+1] <- mean(output2[1:k])
    cat(paste("\n", "AUC for ensemble models based on suggested input weights (using presence and background data sets generated for ", k, "-fold cross-validations)",  "\n", sep = ""))
    print(output2)

    cat(paste("\n", "(Results with input weights inferred from MEAN.T column with similar procedures",  "\n", sep = ""))
    cat(paste("parameters for weighting of MEAN.T: ENSEMBLE.min=", ENSEMBLE.min1, ", ENSEMBLE.best=", ENSEMBLE.best1, " and ENSEMBLE.exponent=", ENSEMBLE.exponent1, "\n\n", sep = ""))
    cat(paste("Final suggested weights with this alternative procedure", "\n", sep = ""))
    print(output.weightsT)

    output3 <- numeric(length=k+1)
    names(output3)[1:k] <- paste("T_", c(1:k), sep="")
    names(output3)[k+1] <- c("MEAN.T")
    for (i in 1:k) {
        TrainData <- TrainData.all[[i]]
        TrainData[,"ENSEMBLE"] <- output.weightsT["MAXENT"]*TrainData[,"MAXENT"] + output.weightsT["MAXNET"]*TrainData[,"MAXNET"] +
            output.weightsT["MAXLIKE"]*TrainData[,"MAXLIKE"] + output.weightsT["GBM"]*TrainData[,"GBM"] +
            output.weightsT["GBMSTEP"]*TrainData[,"GBMSTEP"] + output.weightsT["RF"]*TrainData[,"RF"] + output.weightsT["CF"]*TrainData[,"CF"] + output.weightsT["GLM"]*TrainData[,"GLM"] +
            output.weightsT["GLMSTEP"]*TrainData[,"GLMSTEP"] + output.weightsT["GAM"]*TrainData[,"GAM"] + output.weightsT["GAMSTEP"]*TrainData[,"GAMSTEP"] +
            output.weightsT["MGCV"]*TrainData[,"MGCV"] + output.weightsT["MGCVFIX"]*TrainData[,"MGCVFIX"] + output.weightsT["EARTH"]*TrainData[,"EARTH"] +
            output.weightsT["RPART"]*TrainData[,"RPART"] + output.weightsT["NNET"]*TrainData[,"NNET"] + output.weightsT["FDA"]*TrainData[,"FDA"] +
            output.weightsT["SVM"]*TrainData[,"SVM"] + output.weightsT["SVME"]*TrainData[,"SVME"] + output.weightsT["GLMNET"]*TrainData[,"GLMNET"] + 
            output.weightsT["BIOCLIM.O"]*TrainData[,"BIOCLIM.O"] + output.weightsT["BIOCLIM"]*TrainData[,"BIOCLIM"] +
            output.weightsT["DOMAIN"]*TrainData[,"DOMAIN"] + output.weightsT["MAHAL"]*TrainData[,"MAHAL"]+ output.weightsT["MAHAL01"]*TrainData[,"MAHAL01"]
        TestData <- TestData.all[[i]]
        TestData[,"ENSEMBLE"] <- output.weightsT["MAXENT"]*TestData[,"MAXENT"] + output.weightsT["MAXNET"]*TestData[,"MAXNET"] +
            output.weightsT["MAXLIKE"]*TestData[,"MAXLIKE"] + output.weightsT["GBM"]*TestData[,"GBM"] +
            output.weightsT["GBMSTEP"]*TestData[,"GBMSTEP"] + output.weightsT["RF"]*TestData[,"RF"] + output.weightsT["CF"]*TestData[,"CF"] + output.weightsT["GLM"]*TestData[,"GLM"] +
            output.weightsT["GLMSTEP"]*TestData[,"GLMSTEP"] + output.weightsT["GAM"]*TestData[,"GAM"] + output.weightsT["GAMSTEP"]*TestData[,"GAMSTEP"] +
            output.weightsT["MGCV"]*TestData[,"MGCV"] + output.weightsT["MGCVFIX"]*TestData[,"MGCVFIX"] + output.weightsT["EARTH"]*TestData[,"EARTH"] +
            output.weightsT["RPART"]*TestData[,"RPART"] + output.weightsT["NNET"]*TestData[,"NNET"] + output.weightsT["FDA"]*TestData[,"FDA"] +
            output.weightsT["SVM"]*TestData[,"SVM"] + output.weightsT["SVME"]*TestData[,"SVME"] + output.weightsT["GLMNET"]*TestData[,"GLMNET"] + 
            output.weightsT["BIOCLIM.O"]*TestData[,"BIOCLIM.O"] + output.weightsT["BIOCLIM"]*TestData[,"BIOCLIM"] +
            output.weightsT["DOMAIN"]*TestData[,"DOMAIN"] + output.weightsT["MAHAL"]*TestData[,"MAHAL"]+ output.weightsT["MAHAL01"]*TestData[,"MAHAL01"]
            eval1 <- eval2 <- NULL
            TestPres <- as.numeric(TestData[TestData[,"pb"]==1, "ENSEMBLE"])
            TestAbs <- as.numeric(TestData[TestData[,"pb"]==0, "ENSEMBLE"])
            eval1 <- dismo::evaluate(p=TestPres, a=TestAbs)
        TrainPres <- as.numeric(TrainData[TrainData[,"pb"]==1, "ENSEMBLE"])
        TrainAbs <- as.numeric(TrainData[TrainData[,"pb"]==0, "ENSEMBLE"])
        evalT <- dismo::evaluate(p=TrainPres, a=TrainAbs)
        TestPres <- as.numeric(TestData[TestData[,"pb"]==1, "ENSEMBLE"])
        TestAbs <- as.numeric(TestData[TestData[,"pb"]==0, "ENSEMBLE"])
        eval1 <- dismo::evaluate(p=TestPres, a=TestAbs)
        output3[i] <- eval1@auc
        cat(paste("\n", "Results with alternative final weights for ensemble.evaluate for k = ", i, "\n\n", sep = ""))
        eval3 <- ensemble.evaluate(eval=eval1, eval.train=evalT)
        print(eval3)
    }
    output3[k+1] <- mean(output3[1:k])
    cat(paste("\n", "AUC for ensemble models based on alternative input weights", "\n", sep = ""))
    print(output3)
    cat(paste(")", "\n", sep=""))

    if (SINK==T && OLD.SINK==F) {sink(file=NULL, append=T)}

    if (x.was.terra == TRUE) {x <- terra::rast(x)}
    
    if (data.keep == F) {
        cat(paste("\n\n"))
        return(list(AUC.table=output, table=output, output.weights=output.weights, AUC.with.suggested.weights=output2, 
            eval.table.all=eval.table.all, eval.table.final=eval.table.final,
            x=x, p=p.all, a=a.all, MAXENT.a=MAXENT.a, groupp=groupp, groupa=groupa,
            var.names=var.names, factors=factors2, dummy.vars=dummy.vars2, dummy.vars.noDOMAIN=dummy.vars.noDOMAIN,
            species.name=species.name, 
            call=match.call()))
    }else{
        cat(paste("\n\n"))
        return(list(data=TestData.all, 
            AUC.table=output, table=output, output.weights=output.weights, AUC.with.suggested.weights=output2, 
            data=TestData.all, eval.table.all=eval.table.all, eval.table.final=eval.table.final,
            x=x, p=p.all, a=a.all, MAXENT.a=MAXENT.a, groupp=groupp, groupa=groupa,
            var.names=var.names, factors=factors2, dummy.vars=dummy.vars2, dummy.vars.noDOMAIN=dummy.vars.noDOMAIN,
            species.name=species.name, 
            call=match.call()))
    }
}

