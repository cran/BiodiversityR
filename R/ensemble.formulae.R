`ensemble.formulae` <- function(
    x, factors=NULL, dummy.vars=NULL
)
{
# in older version of raster used layerNames instead of names
    vars <- names(x)
    gamscope <- as.list(vars)
    names(gamscope) <- vars
    nv <- length(vars)
    nf <- length(factors)
    nd <- length(dummy.vars)
    if (is.null(factors) == F) {
        for (i in 1:nf) {
            if (any(vars==factors[i])==FALSE) {
                cat(paste("\n", "WARNING: categorical variable '", factors[i], "' not among grid layers", "\n", "\n", sep = ""))
            }
        }
    }
    if (is.null(dummy.vars) == F) {
        for (i in 1:nd) {
            if (any(vars==dummy.vars[i])==FALSE) {
               cat(paste("\n", "WARNING: dummy variable '", dummy.vars[i], "' not among grid layers", "\n", "\n", sep = ""))
            }
        }
    }
    results <- list(GBM.formula=NULL, RF.formula=NULL, GLM.formula=NULL, STEP.formula=NULL, GLMSTEP.scope=NULL,  
        GAM.formula=NULL, GAMSTEP.scope=NULL, MGCV.formula=NULL, 
        EARTH.formula=NULL, RPART.formula=NULL, NNET.formula=NULL,
        FDA.formula=NULL, SVM.formula=NULL)
    numpb <- paste("pb ~ ")
    catpb <- paste("as.factor(pb) ~ ")
    numvars <- NULL
    stepvars <- NULL
    for (i in 1:nv) {
        if (any(vars[i]==factors)==FALSE) {
            if(is.null(numvars)==T) {
                numvars <- paste(vars[i])
                stepvars <- numvars
            }else{
                numvars <- paste(numvars, "+", vars[i])
            }
        }
    }
    catvars <- NULL
    for (i in 1:nv) {
        if (any(vars[i]==factors)==TRUE) {
            if(is.null(numvars)==T) {
                catvars <- paste(vars[i])
                if (is.null(stepvars)==T) {stepvars <- catvars}
            }else{
                catvars <- paste(catvars, "+", vars[i], sep="")
            }
        }
    }
    explicitcatvars <- NULL
    for (i in 1:nv) {
        if (any(vars[i]==factors)==TRUE) {
            if(is.null(numvars)==T) {
                explicitcatvars <- paste("as.factor(", vars[i], ")", sep="")
            }else{
                explicitcatvars <- paste(explicitcatvars, " + as.factor(", vars[i], ")", sep="")
            }
            gamscope[[as.name(vars[i])]] <- as.formula(paste("~1 + as.factor(", vars[i], ")", sep=""))
        }
    }
    glmvars <- NULL
    for (i in 1:nv) {
        if (any(vars[i]==factors)==FALSE) {
            if(is.null(glmvars)==T) {
                if (any(vars[i]==dummy.vars)==FALSE) {
                    glmvars <- paste(vars[i], "+ I(", vars[i], "^2) + I(", vars[i], "^3)", sep="")
                }else{
                    glmvars <- paste(vars[i])
                }
            }else{
                if (any(vars[i]==dummy.vars)==FALSE) {
                    glmvars <- paste(glmvars, "+", vars[i], "+ I(", vars[i], "^2) + I(", vars[i], "^3)", sep="")
                }else{
                    glmvars <- paste(glmvars, "+", vars[i], sep="")
                }
            }
        }
    }
    gamvars <- NULL
    for (i in 1:nv) {
        if (any(vars[i]==factors)==FALSE) {
            if(is.null(gamvars)==T) {
                if (any(vars[i]==dummy.vars)==FALSE) {
                    gamvars <- paste("s(", vars[i], ",4)", sep="")
                }else{
                    gamvars <- paste(vars[i])
                }
            }else{
                if (any(vars[i]==dummy.vars)==FALSE) {
                    gamvars <- paste(gamvars, "+ s(", vars[i], ",4)", sep="")
                }else{
                    gamvars <- paste(gamvars, "+", vars[i], sep="")
                }
            }
            if (any(vars[i]==dummy.vars)==FALSE) {
                gamscope[[as.name(vars[i])]] <- as.formula(paste("~1 + ", vars[i], "+ s(", vars[i], ", 2) + s(", vars[i], ", 3) + s(", vars[i], ", 4)", sep=""))
            }else{
                gamscope[[as.name(vars[i])]] <- as.formula(paste("~1 + ", vars[i], sep=""))
            }  
        }
    }
    mgcvvars <- NULL
    for (i in 1:nv) {
        if (any(vars[i]==factors)==FALSE) {
            if(is.null(mgcvvars)==T) {
                if (any(vars[i]==dummy.vars)==FALSE) {
                    mgcvvars <- paste("s(", vars[i], ", k=4)", sep="")
                }else{
                    mgcvvars <- paste(vars[i])
                }
            }else{
                if (any(vars[i]==dummy.vars)==FALSE) {
                    mgcvvars <- paste(mgcvvars, "+ s(", vars[i], ", k=4)", sep="")
                }else{
                    mgcvvars <- paste(mgcvvars, "+", vars[i], sep="")
                }
            }
        }
    }
    results$GBM.formula <- as.formula(paste(numpb, numvars, catvars, sep=""))
    results$RF.formula <- as.formula(paste(numpb, numvars, catvars, sep=""))
    results$GLM.formula <- as.formula(paste(numpb, glmvars, catvars, sep=""))
    results$STEP.formula <- as.formula(paste(numpb, stepvars, sep=""))
    results$GLMSTEP.scope <- list(upper=as.formula(paste("~", glmvars, catvars, sep="")), lower=as.formula(paste("~1")))
    results$GAM.formula <- as.formula(paste(numpb, gamvars, explicitcatvars, sep=""))
    results$GAMSTEP.scope <- gamscope
    results$MGCV.formula <- as.formula(paste(numpb, mgcvvars, catvars, sep=""))
# no categorical variables for earth
    results$EARTH.formula <- as.formula(paste(catpb, numvars, sep=""))
    results$RPART.formula <- as.formula(paste(catpb, numvars, catvars, sep=""))
    results$NNET.formula <- as.formula(paste(catpb, numvars, explicitcatvars, sep=""))
    results$FDA.formula <- as.formula(paste(numpb, numvars, catvars, sep=""))
    results$SVM.formula <- as.formula(paste(numpb, numvars, catvars, sep=""))
    return(results)
}

