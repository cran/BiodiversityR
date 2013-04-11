`evaluation.strip.plot2` <- function(
    data, 
    modelnames=c("MAXENT", "GBM", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP",
        "MGCV", "EARTH", "RPART", "NNET", "FDA", "SVM", "BIOCLIM", "DOMAIN", "MAHAL"),
    variable=NULL, model=NULL, 
    pch1=1, pch2=2, col1="black", col2="red", ...
) 
{
    data1 <- data[,which(data[,1]==0)]
    data2 <- data[,which(data[,1]==1)]
    if(is.null(variable)==F) {
        v <- (which(colnames(data) == variable))
        f <- data[,2]==v
        f1 <- data1[,2]==v
        f2 <- data2[,2]==v
        vars <- max(data[,2])
        modelnames <- c(modelnames, "ENSEMBLE")
        nmodels <- length(modelnames)
#    models with data
        models <- 0
        for (j in 1:nmodels) {
            if (any(is.na(data[f, 3+vars+j])==F)) {models <- models + 1}
        }
        dim1 <- ceiling(sqrt(models))
        dim2 <- ceiling(models/dim1)
        par(mfrow=c(dim1,dim2))
        for (j in 1:models) {
            if (any(is.na(data[v, 3+vars+j])==F)) {
                plot(data[f,v+3], data[f, 3+vars+j], type="n", main=variable, xlab="", ylab=colnames(data)[3+vars+j])
                points(data1[f1,v+3], data1[f1, 3+vars+j], type="o", pch=pch1, col=col1, ...)
                points(data2[f2,v+3], data2[f2, 3+vars+j], type="o", pch=pch2, col=col2, ...)
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
            f1 <- data1[,1]==i
            f2 <- data2[,1]==i
            plot(data[f,i+3], data[f, m], type="n", main=colnames(data)[i+3], xlab="", ylab=model)
            points(data1[f1,v+3], data1[f1, 3+vars+j], type="o", pch=pch1, col1=col1, ...)
            points(data2[f2,v+3], data2[f2, 3+vars+j], type="o", pch=pch1, col1=col1, ...)
        }
        par(mfrow=c(1,1))
    }
}




