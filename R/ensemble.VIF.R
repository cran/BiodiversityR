`ensemble.VIF` <- function(
    x=NULL, a=NULL, an=10000, 
    VIF.max=10, 
    layer.drops=NULL, factors=NULL, dummy.vars=NULL
)
{
    if(is.null(x) == T) {stop("value for parameter x is missing (RasterStack object)")}
    if(inherits(x,"RasterStack") == F) {stop("x is not a RasterStack object")}
    if(is.null(a) == F) {names(a) <- c("x", "y")}
#
    layer.drops <- as.character(layer.drops)
    factors <- as.character(factors)
    dummy.vars <- as.character(dummy.vars)
#
    vars <- names(x)
#
    if (length(layer.drops) > 0) {
        var.drops <- layer.drops
        nd <- length(layer.drops)
        for (i in 1:nd) {     
            if (any(vars == layer.drops[i]) == FALSE) {
                cat(paste("\n", "WARNING: variable to exclude '", layer.drops[i], "' not among grid layers", sep = ""))
            }else{
                cat(paste("\n", "NOTE: variable '", layer.drops[i], "' will not be included as explanatory variable", sep = ""))
                x <- raster::dropLayer(x, which(names(x) %in% c(layer.drops[i]) ))
                x <- raster::stack(x)
                vars <- names(x)
                if (length(factors) > 0) {
                    factors <- factors[factors != layer.drops[i]]
                }
                if (length(dummy.vars) > 0) {
                    dummy.vars <- dummy.vars[dummy.vars != layer.drops[i]]
                }
            }
        }
    }else{
        var.drops <- NULL
    }
#
    var.drops <- as.character(var.drops)
#
    vars <- names(x)
    if (length(factors) > 0) {
        var.drops <- c(var.drops, factors)
        for (i in 1:length(factors)) {vars <- vars[which(vars != factors[i])]}
    }
#
    nv <- length(vars)
#
    result <- data.frame(array(dim=c(nv, nv)))
    names(result) <- vars
    row.names(result) <- paste("step_", c(1:nv), sep="")
#
    if (is.null(a) == T) {a <- dismo::randomPoints(x[[1]], n=an, p=NULL, excludep=F)}
# presence locations will not be used, but are required for the ensemble.calibrate.models function
    p <- dismo::randomPoints(x[[1]], n=30, p=NULL, excludep=F)
#
    i <- 0
    VIF.result.max <- Inf

# for ensemble.calibrate.models, need to use NULL again for factors etc.
    if (length(var.drops) == 0) {var.drops <- NULL}
    if (length(dummy.vars) == 0) {dummy.vars <- NULL}

    cat(paste("\n", "Selection of explanatory variables based on the Variance Explanation Factor", sep = ""))
    cat(paste("\n", "Data obtained from  ", nrow(a), " point locations", sep = ""))
    cat(paste("\n", "If some variables have VIF > VIF.max, then the variable with largest VIF is excluded", sep = ""))
    cat(paste("\n", "The procedure stops when all variables have VIF <= VIF.max", "\n", sep = ""))

    while(VIF.result.max >= VIF.max) { 
        VIF.result <- ensemble.calibrate.models(x, p=p, a=a,
            layer.drops=var.drops, factors=NULL,
            VIF=T, AUC.weights=F, ENSEMBLE.tune=F, 
            MAXENT=0, MAXLIKE=0, GBM=0, GBMSTEP=0, RF=0, GLM=0, GLMSTEP=0, GAM=0, 
            GAMSTEP=0, MGCV=0, MGCVFIX=0,EARTH=0, RPART=0, NNET=0, FDA=0, 
            SVM=0, SVME=0, GLMNET=0,
            BIOCLIM.O=0, BIOCLIM=0, DOMAIN=0, MAHAL=0, MAHAL01=0)$VIF
        i <- i+1
        for (v in 1:length(VIF.result)) {result[i, which(names(result) == names(VIF.result)[v])] <- VIF.result[which(names(VIF.result) == names(VIF.result)[v])]}
        VIF.result.max <- VIF.result[1]
        var.drops <- c(var.drops, names(VIF.result)[1])     
    }

# remove last variable included
    if (length(var.drops) == 1) {
        var.drops <- as.character(NULL)
    }else{
        nvd <- length(var.drops)-1
        var.drops <- var.drops[1:nvd]
    }

    # include factors again as no information to exclude (drop)
    if (length(factors) > 0) {
        for (i in 1:length(factors)) {var.drops <- var.drops[which(var.drops != factors[i])]}
        vars.included <- c(names(VIF.result), factors)
    }else{
        vars.included <- names(VIF.result)
        factors <- NULL
    }

    dummy.vars <- as.character(dummy.vars)
    if (length(dummy.vars) > 0) {
        for (i in 1:length(var.drops)) {
            if (var.drops[i] %in% dummy.vars) {dummy.vars <- dummy.vars[which(dummy.vars != var.drops[i])]}
        }
    }else{
        dummy.vars <- NULL
    }

    if (length(var.drops) == 0) {var.drops <- NULL}

    result <- result[rowSums(result, na.rm=T) > 0, , drop=F]

    cat(paste("Summary of VIF selection process:", "\n", sep = ""))
    print(result)
    cat(paste("\n", sep = ""))

    cat(paste("Final selection of variables:", "\n", sep = ""))
    print(vars.included)
    cat(paste("\n", sep = ""))

    result <- result[rowSums(result, na.rm=T) > 0, ]
    return(list(stepwise.results=result, var.drops=var.drops, vars.included=vars.included,
        factors=factors, dummy.vars.included=dummy.vars, VIF.final=VIF.result))
}

