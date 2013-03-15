`ensemble.drop1` <- function(
    x, p, a=NULL, an=1000, ext=NULL, k=0, pt=NULL, at=NULL,
    layer.drops=NULL, VIF=FALSE,
    digits=2, difference=FALSE,
    input.weights=NULL,
    MAXENT=1, GBM=1, GBMSTEP=1, RF=1, GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=1, MGCV=1, MGCVFIX=0, 
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, BIOCLIM=1, DOMAIN=1, MAHAL=1, 
    Yweights="BIOMOD", factors=NULL, dummy.vars=NULL,
    maxit=100,
    GBM.n.trees=2000,
    GBMSTEP.tree.complexity=5, GBMSTEP.learning.rate=0.005, 
    GBMSTEP.bag.fraction=0.5, GBMSTEP.step.size=100,
    RF.ntree=750, RF.mtry=log(nlayers(x)), 
    GLM.family=binomial(link="logit"), 
    GLMSTEP.steps=1000, GLMSTEP.k=2,
    GAM.family=binomial(link="logit"), 
    GAMSTEP.steps=1000,
    MGCV.select=FALSE, 
    EARTH.glm=list(family=binomial(link="logit"), maxit=maxit),
    RPART.xval=50, 
    NNET.size=8, NNET.decay=0.01
  )
{
    .BiodiversityR <- new.env()
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
# no point for GEODIST
    GEODIST <- -1
    if (MAXENT > 0) {
	jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
	if (!file.exists(jar)) {stop('maxent program is missing: ', jar, '\nPlease download it here: http://www.cs.princeton.edu/~schapire/maxent/')}
    }
    if (GBM > 0) {
        if (! require(gbm)) {stop("Please install the gbm package")}
    }
    if (GBMSTEP > 0) {
        if (! require(gbm)) {stop("Please install the gbm package")}
    }
    if (RF > 0) {
        if (! require(randomForest)) {stop("Please install the randomForest package")}
    }
    if (GLMSTEP > 0) {
        if (! require(MASS)) {stop("Please install the MASS package")}
    }
    if (GAM > 0) {
        if (! require(gam)) {stop("Please install the gam package")}
        detach(package:gam)
    }
    if (GAMSTEP > 0) {
        if (! require(gam)) {stop("Please install the gam package")}
        detach(package:gam)   
    }
    if (MGCV > 0) {
        if (! require(mgcv)) {stop("Please install the mgcv package")}
        detach(package:mgcv)   
    }
    if (MGCVFIX > 0) {
        if (! require(mgcv)) {stop("Please install the mgcv package")}
        detach(package:mgcv)   
    }
    if (EARTH > 0) {
        if (! require(earth)) {stop("Please install the earth package")}
    }
    if (RPART > 0) {
        if (! require(rpart)) {stop("Please install the rpart package")}
    }
    if (NNET > 0) {
        if (! require(nnet)) {stop("Please install the nnet package")}
    }
    if (FDA > 0) {
        if (! require(mda)) {stop("Please install the mda package")}
    }
    if (SVM > 0) {
        if (! require(kernlab)) {stop("Please install the kernlab package")}
    }
    if (SVME > 0) {
        if (! require(e1071)) {stop("Please install the e1071 package")}
#         get the probabilities from svm
            predict.svme <- function(model, newdata) {
                p <- predict(model, newdata, probability=T)
                return(attr(p, "probabilities")[,1])
             }
    }
    if (is.null(a)==T) {
        a <- randomPoints(x, n=an, ext=ext)
    }
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
    assign("p", p, envir=.BiodiversityR)
    assign("pt", pt, envir=.BiodiversityR)
    assign("a", a, envir=.BiodiversityR)
    assign("at", at, envir=.BiodiversityR)
    vars <- names(x)
    nv <- length(vars)
    output <- array(NA, dim=c(19, nv+1))
    rownames(output) <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL")
    colnames(output) <- c("all_vars", paste("without_",vars))
# first fit with all variables
    tests <- ensemble.test(x=x, p=p, a=a, pt=pt, at=at, 
        PLOTS=FALSE, evaluations.keep=T, models.keep=F,
        layer.drops=NULL, VIF=VIF,
        formulae.defaults=T, maxit=maxit,
        MAXENT=MAXENT, GBM=GBM, GBMSTEP=GBMSTEP, RF=RF, GLM=GLM, GLMSTEP=GLMSTEP, 
        GAM=GAM, GAMSTEP=GAMSTEP, MGCV=MGCV, MGCVFIX=MGCVFIX, EARTH=EARTH, RPART=RPART, 
        NNET=NNET, FDA=FDA, SVM=SVM, SVME=SVME, BIOCLIM=BIOCLIM, DOMAIN=DOMAIN, MAHAL=MAHAL,
        GEODIST=0, 
        Yweights=Yweights, factors=factors, dummy.vars=dummy.vars,
        GBM.formula=NULL, GBM.n.trees=GBM.n.trees,
        GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
        GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction,
        GBMSTEP.step.size=GBMSTEP.step.size,
        RF.formula=NULL, RF.ntree=RF.ntree, RF.mtry=RF.mtry, 
        GLM.formula=NULL, GLM.family=GLM.family, 
        GLMSTEP.k=GLMSTEP.k, GLMSTEP.steps=GLMSTEP.steps, STEP.formula=NULL, GLMSTEP.scope=NULL, 
        GAM.formula=NULL, GAM.family=GAM.family, 
        GAMSTEP.steps=GAMSTEP.steps, GAMSTEP.scope=NULL,
        MGCV.formula=NULL, MGCV.select=MGCV.select,
        MGCVFIX.formula=NULL, 
        EARTH.formula=NULL, EARTH.glm=EARTH.glm,
        RPART.formula=NULL, RPART.xval=RPART.xval, 
        NNET.formula=NULL, NNET.size=NNET.size, NNET.decay=NNET.decay,
        FDA.formula=NULL, SVM.formula=NULL, SVME.formula=NULL)
    if(is.null(tests$MAXENT.T)==F) {output["MAXENT",1] <- tests$MAXENT.T@auc}
    if(is.null(tests$GBM.T)==F) {output["GBM",1] <- tests$GBM.T@auc} 
    if(is.null(tests$GBMSTEP.T)==F) {output["GBMSTEP",1] <- tests$GBMSTEP.T@auc} 
    if(is.null(tests$RF.T)==F) {output["RF",1] <- tests$RF.T@auc}
    if(is.null(tests$GLM.T)==F) {output["GLM",1] <- tests$GLM.T@auc} 
    if(is.null(tests$GLMS.T)==F) {output["GLMSTEP",1] <- tests$GLMS.T@auc}
    if(is.null(tests$GAM.T)==F) {output["GAM",1] <- tests$GAM.T@auc} 
    if(is.null(tests$GAMS.T)==F) {output["GAMSTEP",1] <- tests$GAMS.T@auc}
    if(is.null(tests$MGCV.T)==F) {output["MGCV",1] <- tests$MGCV.T@auc} 
    if(is.null(tests$MGCVF.T)==F) {output["MGCVFIX",1] <- tests$MGCVF.T@auc} 
    if(is.null(tests$EARTH.T)==F) {output["EARTH",1] <- tests$EARTH.T@auc} 
    if(is.null(tests$RPART.T)==F) {output["RPART",1] <- tests$RPART.T@auc}
    if(is.null(tests$NNET.T)==F) {output["NNET",1] <- tests$NNET.T@auc} 
    if(is.null(tests$FDA.T)==F) {output["FDA",1] <- tests$FDA.T@auc}
    if(is.null(tests$SVM.T)==F) {output["SVM",1] <- tests$SVM.T@auc}
    if(is.null(tests$SVME.T)==F) {output["SVME",1] <- tests$SVME.T@auc}
    if(is.null(tests$BIOCLIM.T)==F) {output["BIOCLIM",1] <- tests$BIOCLIM.T@auc}
    if(is.null(tests$DOMAIN.T)==F) {output["DOMAIN",1] <- tests$DOMAIN.T@auc}
    if(is.null(tests$MAHAL.T)==F) {output["MAHAL",1] <- tests$MAHAL.T@auc}




# sequentially leave out the focal variable, then fit again

    vars <- names(x)
    nv <- length(vars)
    for (i in 1:nv) {
        var.f <- vars[i]
        cat(paste("\n", "Leaving out variable ", var.f, "\n\n", sep = ""))
    tests <- ensemble.test(x=x, p=p, a=a, pt=pt, at=at, 
        PLOTS=FALSE, evaluations.keep=T, models.keep=F,
        layer.drops=var.f, VIF=VIF,
        formulae.defaults=T, maxit=maxit,
        MAXENT=MAXENT, GBM=GBM, GBMSTEP=GBMSTEP, RF=RF, GLM=GLM, GLMSTEP=GLMSTEP, 
        GAM=GAM, GAMSTEP=GAMSTEP, MGCV=MGCV, MGCVFIX=MGCVFIX, EARTH=EARTH, RPART=RPART, 
        NNET=NNET, FDA=FDA, SVM=SVM, SVME=SVME, BIOCLIM=BIOCLIM, DOMAIN=DOMAIN, MAHAL=MAHAL,
        GEODIST=0,
        Yweights=Yweights, factors=factors, dummy.vars=dummy.vars,
        GBM.formula=NULL, GBM.n.trees=GBM.n.trees,
        GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
        GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction,
        GBMSTEP.step.size=GBMSTEP.step.size,
        RF.formula=NULL, RF.ntree=RF.ntree, RF.mtry=RF.mtry, 
        GLM.formula=NULL, GLM.family=GLM.family, 
        GLMSTEP.k=GLMSTEP.k, GLMSTEP.steps=GLMSTEP.steps, STEP.formula=NULL, GLMSTEP.scope=NULL, 
        GAM.formula=NULL, GAM.family=GAM.family, 
        GAMSTEP.steps=GAMSTEP.steps, GAMSTEP.scope=NULL,
        MGCV.formula=NULL, MGCV.select=MGCV.select,
        MGCVFIX.formula=NULL, 
        EARTH.formula=NULL, EARTH.glm=EARTH.glm,
        RPART.formula=NULL, RPART.xval=RPART.xval, 
        NNET.formula=NULL, NNET.size=NNET.size, NNET.decay=NNET.decay,
        FDA.formula=NULL, SVM.formula=NULL, SVME.formula=NULL)
    if(is.null(tests$MAXENT.T)==F) {output["MAXENT",i+1] <- tests$MAXENT.T@auc}
    if(is.null(tests$GBM.T)==F) {output["GBM",i+1] <- tests$GBM.T@auc} 
    if(is.null(tests$GBMSTEP.T)==F) {output["GBMSTEP",i+1] <- tests$GBMSTEP.T@auc} 
    if(is.null(tests$RF.T)==F) {output["RF",i+1] <- tests$RF.T@auc}
    if(is.null(tests$GLM.T)==F) {output["GLM",i+1] <- tests$GLM.T@auc} 
    if(is.null(tests$GLMS.T)==F) {output["GLMSTEP",i+1] <- tests$GLMS.T@auc}
    if(is.null(tests$GAM.T)==F) {output["GAM",i+1] <- tests$GAM.T@auc} 
    if(is.null(tests$GAMS.T)==F) {output["GAMSTEP",i+1] <- tests$GAMS.T@auc}
    if(is.null(tests$MGCV.T)==F) {output["MGCV",i+1] <- tests$MGCV.T@auc} 
    if(is.null(tests$MGCVF.T)==F) {output["MGCVFIX",i+1] <- tests$MGCVF.T@auc} 
    if(is.null(tests$EARTH.T)==F) {output["EARTH",i+1] <- tests$EARTH.T@auc} 
    if(is.null(tests$RPART.T)==F) {output["RPART",i+1] <- tests$RPART.T@auc}
    if(is.null(tests$NNET.T)==F) {output["NNET",i+1] <- tests$NNET.T@auc} 
    if(is.null(tests$FDA.T)==F) {output["FDA",i+1] <- tests$FDA.T@auc}
    if(is.null(tests$SVM.T)==F) {output["SVM",i+1] <- tests$SVM.T@auc}
    if(is.null(tests$SVME.T)==F) {output["SVME",i+1] <- tests$SVME.T@auc}
    if(is.null(tests$BIOCLIM.T)==F) {output["BIOCLIM",i+1] <- tests$BIOCLIM.T@auc}
    if(is.null(tests$DOMAIN.T)==F) {output["DOMAIN",i+1] <- tests$DOMAIN.T@auc}
    if(is.null(tests$MAHAL.T)==F) {output["MAHAL",i+1] <- tests$MAHAL.T@auc}

    }
    output <- 100*output
    output <- round(output, digits=digits)
    if (difference == T) {
        for (i in 1:nv) {
            output[,i+1] <- output[,i+1] - output[,1]
        }
    }
    cat(paste("\n", "Results (as percentage) of ensemble.drop1",  "\n\n", sep = ""))
    if (difference == T) {
        cat(paste("\n", "Note that positive differences indicate that the model without the variable", sep = ""))
        cat(paste("\n", "has higher AUC than the the model with all the variables",  "\n\n", sep = ""))
    }
    print (output)
    remove(p, envir=.BiodiversityR)
    remove(pt, envir=.BiodiversityR)
    remove(a, envir=.BiodiversityR)
    remove(at, envir=.BiodiversityR)
    return(output)

}