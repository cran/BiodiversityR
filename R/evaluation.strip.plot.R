`evaluation.strip.plot` <- function(
    data, 
    modelnames=c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP",
        "MGCV", "EARTH", "RPART", "NNET", "FDA", "SVM", "BIOCLIM", "DOMAIN", "MAHAL"),
    variable=NULL, model=NULL, ...
) 
{
    if(is.null(variable)==F) {
        v <- (which(colnames(data) == variable))
        f <- data[,1]==v
        vars <- max(data[,1])
        modelnames <- c(modelnames, "ENSEMBLE")
        nmodels <- length(modelnames)
#    models with data
        models <- 0
        for (j in 1:nmodels) {
            if (any(is.na(data[f, 2+vars+j])==F)) {models <- models + 1}
        }
        dim1 <- ceiling(sqrt(models))
        dim2 <- ceiling(models/dim1)
        par(mfrow=c(dim1,dim2))
        for (j in 1:models) {
            if (any(is.na(data[v, 2+vars+j])==F)) {
                plot(data[f,v+2], data[f, 2+vars+j], main=variable, xlab="", ylab=colnames(data)[2+vars+j],...)
            }
        }
        par(mfrow=c(1,1))
    }
    if(is.null(model)==F) {
        m <- colnames(data) == model
        vars <- max(data[,1])
        dim1 <- ceiling(sqrt(vars))
        dim2 <- ceiling(vars/dim1)
        par(mfrow=c(dim1,dim2))
        for (i in 1:vars) {
            f <- data[,1]==i
            plot(data[f,i+2], data[f, m], main=colnames(data)[i+2], xlab="", ylab=model,...)
        }
        par(mfrow=c(1,1))
    }
}




