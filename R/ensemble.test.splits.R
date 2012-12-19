`ensemble.test.splits` <- function(
    x, p, a=NULL, an=1000, ext=NULL, k=5, digits=2, PLOTS=FALSE, SCRIPT=TRUE,
    MAXENT=1, GBM=1, GBMSTEP=1, RF=1, GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=1, 
    MGCV=1, EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, BIOCLIM=1, DOMAIN=1, MAHAL=1,   
    Yweights="BIOMOD", factors=NULL, formulae.defaults=TRUE,
    GBM.formula=NULL, GBM.n.trees=3000,
    GBMSTEP.gbm.x=2:(1+nlayers(x)), GBMSTEP.tree.complexity=5, GBMSTEP.learning.rate=0.01, GBMSTEP.bag.fraction=0.5,
    RF.formula=NULL, RF.ntree=750, RF.mtry=log(nlayers(x)), 
    GLM.formula=NULL, GLM.family=binomial(link="logit"), 
    GLMSTEP.steps=1000, STEP.formula=NULL, GLMSTEP.scope=NULL, GLMSTEP.k=2,
    GAM.formula=NULL, GAM.family=binomial(link="logit"), 
    GAMSTEP.steps=1000, GAMSTEP.scope=NULL,
    MGCV.formula=NULL,
    EARTH.formula=NULL, EARTH.glm=list(family=binomial(link="logit"), maxit=1000),
    RPART.formula=NULL, RPART.xval=50, 
    NNET.formula=NULL, NNET.size=8, NNET.decay=0.01,
    FDA.formula=NULL,
    SVM.formula=NULL  
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
    groupp <- kfold(p, k=k)
    groupa <- kfold(a, k=k)
    output <- array(NA, dim=c(17, k+1))
    rownames(output) <- c("MAXENT","GBM", "GBMSTEP","RF","GLM", "GLMSTEP", "GAM", "GAMSTEP",
        "MGCV","EARTH","RPART","NNET","FDA","SVM","BIOCLIM","DOMAIN","MAHAL")
    colnames(output) <- c(1:k,"MEAN")
    for (i in 1:k){
        cat(paste("\n", "EVALUATION RUN: ", i, "\n", "\n", sep = ""))
        pc <- p[groupp != i,]
        pt <- p[groupp == i,]
        ac <- a[groupa != i,]
        at <- a[groupa == i,]
        assign("pc", pc, envir=.BiodiversityR)
        assign("pt", pt, envir=.BiodiversityR)
        assign("ac", ac, envir=.BiodiversityR)
        assign("at", at, envir=.BiodiversityR)
        tests <- ensemble.test(x=x, p=pc, a=ac, pt=pt, at=at, 
            PLOTS=PLOTS, evaluations.keep=T,
            MAXENT=MAXENT, GBM=GBM, GBMSTEP=GBMSTEP, RF=RF, GLM=GLM, GLMSTEP=GLMSTEP, 
            GAM=GAM, GAMSTEP=GAMSTEP, MGCV=MGCV, EARTH=EARTH, RPART=RPART, 
            NNET=NNET, FDA=FDA, BIOCLIM=BIOCLIM, DOMAIN=DOMAIN, MAHAL=MAHAL,   
            Yweights=Yweights, factors=factors,
            GBM.formula=GBM.formula, GBM.n.trees=GBM.n.trees,
            GBMSTEP.gbm.x=GBMSTEP.gbm.x, GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
            GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction,
            RF.formula=RF.formula, RF.ntree=RF.ntree, RF.mtry=RF.mtry, 
            GLM.formula=GLM.formula, GLM.family=GLM.family, 
            GLMSTEP.k=GLMSTEP.k, GLMSTEP.steps=GLMSTEP.steps, STEP.formula=STEP.formula, GLMSTEP.scope=GLMSTEP.scope, 
            GAM.formula=GAM.formula, GAM.family=GAM.family, 
            GAMSTEP.steps=GAMSTEP.steps, GAMSTEP.scope=GAMSTEP.scope,
            MGCV.formula=MGCV.formula,
            EARTH.formula=EARTH.formula, EARTH.glm=EARTH.glm,
            RPART.formula=RPART.formula, RPART.xval=RPART.xval, 
            NNET.formula=NNET.formula, NNET.size=NNET.size, NNET.decay=NNET.decay,
            FDA.formula=FDA.formula, SVM.formula=SVM.formula)
        if(is.null(tests$MAXENT.T)==F) {output["MAXENT",i] <- tests$MAXENT.T@auc}
        if(is.null(tests$GBM.T)==F) {output["GBM",i] <- tests$GBM.T@auc} 
        if(is.null(tests$GBMSTEP.T)==F) {output["GBMSTEP",i] <- tests$GBMSTEP.T@auc} 
        if(is.null(tests$RF.T)==F) {output["RF",i] <- tests$RF.T@auc}
        if(is.null(tests$GLM.T)==F) {output["GLM",i] <- tests$GLM.T@auc} 
        if(is.null(tests$GLMS.T)==F) {output["GLMSTEP",i] <- tests$GLMS.T@auc}
        if(is.null(tests$GAM.T)==F) {output["GAM",i] <- tests$GAM.T@auc} 
        if(is.null(tests$GAMS.T)==F) {output["GAMSTEP",i] <- tests$GAMS.T@auc}
        if(is.null(tests$MGCV.T)==F) {output["MGCV",i] <- tests$MGCV.T@auc} 
        if(is.null(tests$EARTH.T)==F) {output["EARTH",i] <- tests$EARTH.T@auc} 
        if(is.null(tests$RPART.T)==F) {output["RPART",i] <- tests$RPART.T@auc}
        if(is.null(tests$NNET.T)==F) {output["NNET",i] <- tests$NNET.T@auc} 
        if(is.null(tests$FDA.T)==F) {output["FDA",i] <- tests$FDA.T@auc}
        if(is.null(tests$SVM.T)==F) {output["SVM",i] <- tests$SVM.T@auc}
        if(is.null(tests$BIOCLIM.T)==F) {output["BIOCLIM",i] <- tests$BIOCLIM.T@auc}
        if(is.null(tests$DOMAIN.T)==F) {output["DOMAIN",i] <- tests$DOMAIN.T@auc}
        if(is.null(tests$MAHAL.T)==F) {output["MAHAL",i] <- tests$MAHAL.T@auc}
    }
    output[,k+1] <- rowMeans(output[,1:k], na.rm=T)
    outp <- round(output, digits=4)
    outp[is.na(outp[,"MEAN"]),] <- 0  
    if (SCRIPT==TRUE) {
        cat(paste("\n", "Suggested parameters for ensemble.grd function", "\n", sep = ""))
        cat(paste("\n", "(these are based on average AUC)", "\n", "\n", sep=""))
        cat(paste("\n", "MAXENT=", outp[1, "MEAN"], ", GBM=", outp[2, "MEAN"], ", GBMSTEP=", outp[3, "MEAN"], ", RF=", outp[4, "MEAN"], 
            ", GLM=", outp[5, "MEAN"], ", GLMSTEP=", outp[6, "MEAN"],  
            ", GAM=", outp[7, "MEAN"], ", GAMSTEP=", outp[8, "MEAN"], ",", "\n",
            "MGCV=", outp[9, "MEAN"], ", EARTH=", outp[10, "MEAN"], ", RPART=", outp[11, "MEAN"],
            ", NNET=", outp[12, "MEAN"], ", FDA=", outp[13, "MEAN"], ", SVM=", outp[14, "MEAN"], ", BIOCLIM=", outp[15, "MEAN"],
            ", DOMAIN=", outp[16, "MEAN"], ", MAHAL=", outp[17, "MEAN"], "\n", "\n", "\n", sep=""))
    }
    output <- output[order(output[,"MEAN"], decreasing=T),]
    output <- 100*output
    output <- round(output, digits=digits)
    remove(pc, envir=.BiodiversityR)
    remove(pt, envir=.BiodiversityR)
    remove(ac, envir=.BiodiversityR)
    remove(at, envir=.BiodiversityR)
    return(output)
}




