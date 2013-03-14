`ensemble.test.splits` <- function(
    x, p, a=NULL, an=1000, ext=NULL, k=5, 
    layer.drops=NULL, VIF=FALSE,
    digits=2, PLOTS=FALSE, SCRIPT=TRUE,
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
    if (k < 1) {stop("Parameter k smaller than 1")}
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
    if (is.null(a)==T) {a <- randomPoints(x, n=an, ext=ext)}
    groupp <- kfold(p, k=k)
    groupa <- kfold(a, k=k)
    output <- array(NA, dim=c(20, k+1))
    rownames(output) <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL", "GEODIST")
    colnames(output) <- c(1:k,"MEAN")
    for (i in 1:k){
        cat(paste("\n", "EVALUATION RUN: ", i, "\n\n", sep = ""))
        pc <- p[groupp != i,]
        pt <- p[groupp == i,]
        ac <- a[groupa != i,]
        at <- a[groupa == i,]
        assign("pc", pc, envir=.BiodiversityR)
        assign("pt", pt, envir=.BiodiversityR)
        assign("ac", ac, envir=.BiodiversityR)
        assign("at", at, envir=.BiodiversityR)
        tests <- ensemble.test(x=x, p=pc, a=ac, pt=pt, at=at, 
            VIF=VIF,
            PLOTS=PLOTS, evaluations.keep=T, models.keep=F,
            MAXENT=MAXENT, GBM=GBM, GBMSTEP=GBMSTEP, RF=RF, GLM=GLM, GLMSTEP=GLMSTEP, 
            GAM=GAM, GAMSTEP=GAMSTEP, MGCV=MGCV, MGCVFIX=MGCVFIX, EARTH=EARTH, RPART=RPART, 
            NNET=NNET, FDA=FDA, SVM=SVM, SVME=SVME, BIOCLIM=BIOCLIM, DOMAIN=DOMAIN, MAHAL=MAHAL,
            GEODIST=GEODIST, GEODIST.file.name=GEODIST.file.name, RASTER.format=RASTER.format,  
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
            GAMSTEP.steps=GAMSTEP.steps, GAMSTEP.scope=GAMSTEP.scope, GAMSTEP.pos=GAMSTEP.pos,
            MGCV.formula=MGCV.formula, MGCV.select=MGCV.select,
            MGCVFIX.formula=MGCVFIX.formula, 
            EARTH.formula=EARTH.formula, EARTH.glm=EARTH.glm,
            RPART.formula=RPART.formula, RPART.xval=RPART.xval, 
            NNET.formula=NNET.formula, NNET.size=NNET.size, NNET.decay=NNET.decay,
            FDA.formula=FDA.formula, SVM.formula=SVM.formula, SVME.formula=SVME.formula)
        if(is.null(tests$MAXENT.T)==F) {output["MAXENT",i] <- tests$MAXENT.T@auc}
        if(is.null(tests$GBM.T)==F) {output["GBM",i] <- tests$GBM.T@auc} 
        if(is.null(tests$GBMSTEP.T)==F) {output["GBMSTEP",i] <- tests$GBMSTEP.T@auc} 
        if(is.null(tests$RF.T)==F) {output["RF",i] <- tests$RF.T@auc}
        if(is.null(tests$GLM.T)==F) {output["GLM",i] <- tests$GLM.T@auc} 
        if(is.null(tests$GLMS.T)==F) {output["GLMSTEP",i] <- tests$GLMS.T@auc}
        if(is.null(tests$GAM.T)==F) {output["GAM",i] <- tests$GAM.T@auc} 
        if(is.null(tests$GAMS.T)==F) {output["GAMSTEP",i] <- tests$GAMS.T@auc}
        if(is.null(tests$MGCV.T)==F) {output["MGCV",i] <- tests$MGCV.T@auc} 
        if(is.null(tests$MGCVF.T)==F) {output["MGCVFIX",i] <- tests$MGCVF.T@auc} 
        if(is.null(tests$EARTH.T)==F) {output["EARTH",i] <- tests$EARTH.T@auc} 
        if(is.null(tests$RPART.T)==F) {output["RPART",i] <- tests$RPART.T@auc}
        if(is.null(tests$NNET.T)==F) {output["NNET",i] <- tests$NNET.T@auc} 
        if(is.null(tests$FDA.T)==F) {output["FDA",i] <- tests$FDA.T@auc}
        if(is.null(tests$SVM.T)==F) {output["SVM",i] <- tests$SVM.T@auc}
        if(is.null(tests$SVME.T)==F) {output["SVME",i] <- tests$SVME.T@auc}
        if(is.null(tests$BIOCLIM.T)==F) {output["BIOCLIM",i] <- tests$BIOCLIM.T@auc}
        if(is.null(tests$DOMAIN.T)==F) {output["DOMAIN",i] <- tests$DOMAIN.T@auc}
        if(is.null(tests$MAHAL.T)==F) {output["MAHAL",i] <- tests$MAHAL.T@auc}
        if(is.null(tests$GEODIST.T)==F) {output["GEODIST",i] <- tests$GEODIST.T@auc}
    }
    output[,k+1] <- rowMeans(output[,1:k], na.rm=T)
    output[is.na(output[,"MEAN"]),] <- 0  
    output <- 100*output
    output <- round(output, digits=digits)
    if (SCRIPT==TRUE) {
        weights <- output[,"MEAN"]
# no GEODIST
        weights <- weights[-20]
# AUC of 0.5 indicates random predictions
        weights <- weights - 50
        weights[weights < 0] <- 0
    }
    output <- output[order(output[,"MEAN"], decreasing=T),]
    cat(paste("\n", "Results (as percentage) of ensemble.test.splits sorted by average AUC ",  "\n\n", sep = ""))
    print(output)
    if (SCRIPT==TRUE) {
        cat(paste("\n", "Suggested parameters for ensemble.grd function", "\n", sep = ""))
        cat(paste("(these are based on average AUC minus 0.5)", "\n\n", sep=""))
        print(weights)
    }
    remove(pc, envir=.BiodiversityR)
    remove(pt, envir=.BiodiversityR)
    remove(ac, envir=.BiodiversityR)
    remove(at, envir=.BiodiversityR)
    return(output)
}




