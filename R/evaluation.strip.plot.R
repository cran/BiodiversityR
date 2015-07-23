`evaluation.strip.plot` <- function(
    data, 
    modelnames=c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL"),
    variable=NULL, model=NULL, ...
) 
{
    if(is.null(variable)==F) {
        v <- (which(names(data) == variable))
        v <- v-2
        f <- data[,1]==v
        vars <- max(data[,1])
# plot for all models
        if (is.null(model) ==T) {
            modelnames <- c(modelnames, "ENSEMBLE")
            nmodels <- length(modelnames)
# models with data
            models <- 0
            for (j in 1:nmodels) {
                if (any(is.na(data[f, 2+vars+j])==F)) {models <- models + 1}
            }
            dim1 <- max(1, ceiling(sqrt(models)))
            dim2 <- max(1, ceiling(models/dim1))
            graphics::par(mfrow=c(dim1,dim2))
            for (j in 1:models) {
                if (any(is.na(data[v, 2+vars+j])==F)) {
                    graphics::plot(data[f,v+2], data[f, 2+vars+j], main=variable, xlab="", ylab=names(data)[2+vars+j],...)
                }
            }
            graphics::par(mfrow=c(1,1))
        }else{
            m <- names(data) == model
            if (any(is.na(data[v, m])==F)) {
                graphics::plot(data[f,v+2], data[f, m], main=variable, xlab="", ylab=names(data)[m], ...)
            }
        }
    }
    if(is.null(model)==F && is.null(variable)==T) {
        m <- names(data) == model
# models with data
        if(is.na(sum(data[, m]))) { 
            cat(paste("NOTE: No data for model: ",  model, "\n", sep = ""))
        }else{
            vars <- max(data[,1])
            dim1 <- max(1, ceiling(sqrt(vars)))
            dim2 <- max(1, ceiling(vars/dim1))
            graphics::par(mfrow=c(dim1,dim2))
            for (i in 1:vars) {
                f <- data[,1]==i
                graphics::plot(data[f,i+2], data[f, m], main=names(data)[i+2], xlab="", ylab=model,...)
            }
            graphics::par(mfrow=c(1,1))
        }
    }
}

