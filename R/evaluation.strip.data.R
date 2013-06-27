`evaluation.strip.data` <- function(
    x, ext=NULL, vars=names(x), factors=NULL, steps=50,
    modelnames=c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL")
)
{
    if (is.null(ext) == F) {x <- crop(x, y=ext, snap="in")}
# set minimum and maximum values
    for (i in 1:nlayers(x)) {
        x[[i]] <- setMinMax(x[[i]])
    }
    varslayers <- names(x)
    nvars <- length(vars)
    vars2 <- vars
    for (i in 1:nvars) {
        if(any(vars[i] == varslayers) == F) {
            cat(paste("\n", "Warning: ", vars[i], " not within rasterStack", "\n", sep=""))
            vars2 <- vars2[which(vars2 != vars[i])]
        }
    }
    vars <- vars2
    nvars <- length(vars)
    factors2 <- factors
    if (is.null(factors) == F) {
        factors <- as.character(factors)
        for (i in 1:length(factors)) {
            if(any(factors[i] == vars) == F) {
                cat(paste("\n", "Warning: ", factors[i], " not among variables", "\n", "\n", sep=""))
                factors2 <- factors2[which(factors2 != factors[i])]
            }
        }
    }
    factors <- factors2
    nnum <- nvars - length(factors)
    nrows <- nnum * steps
    out <- array(dim=c(nnum*steps, nvars+2), NA)
    dimnames(out)[[2]] <- c("focal.var","categorical",vars)
# for categorical variables first
    fixedlevel <- array(dim=c(nvars))
    for (i in 1:nvars) {
        if(any(vars[i] == factors) == T) {
            il <- which(varslayers == vars[i])       
            tabulation <- data.frame(freq(x[[il]]))
            NA.index <- !is.na(tabulation[,"value"])
            tabulation <- tabulation[NA.index,]           
            out2 <- array(dim=c(nrow(tabulation), nvars+2), NA)
            out2[, 1] <- rep(i, nrow(tabulation))
            out2[, 2] <- rep(1, nrow(tabulation))
            out2[, i+2] <- tabulation[,1]
            out <- rbind(out, out2)
            fixedlevel[i] <- tabulation[which.max(tabulation[,2]),1]
        }
    }
    for (i in 1:nvars) {
        if(any(vars[i] == factors) == T) {
            index <- is.na(out[,i+2])
            out[index,i+2] <- fixedlevel[i]
        }
    }
    nrows <- nrow(out)
# for numerical variables next
    for (i in 1:nvars) {
        if(any(vars[i] == factors) == F) {
            il <- which(varslayers == vars[i]) 
            out[,i+2] <- rep(cellStats(x[[il]], stat="mean"), nrows)
        }
    }
    j <- 0
    for (i in 1:nvars) {
        if(any(vars[i] == factors) == F) {
            il <- which(varslayers == vars[i]) 
            j <- j+1
            startpos <- (j-1)*steps+1
            endpos <- (j-1)*steps+steps
            out[startpos:endpos,1] <- rep(i, steps)
            out[startpos:endpos,2] <- rep(0, steps)
            minv <- minValue(x[[il]])
            maxv <- maxValue(x[[il]])
            out[startpos:endpos,i+2] <- seq(from=minv,to=maxv, length.out=steps)
         }
    }
    modelnames <- c(modelnames, "ENSEMBLE")
    nmodels <- length(modelnames)
    modelout <- array(dim=c(nrows, nmodels), 0.0)
    dimnames(modelout)[[2]] <- modelnames
    out <- cbind(out, modelout)
    out <- data.frame(out)
    for (i in 1:nvars) {
        if(any(vars[i] == factors) == T) {
            out[,i+2] <- as.factor(out[,i+2])
        }
    }
    return(out)
}




