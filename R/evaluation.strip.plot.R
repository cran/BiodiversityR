`evaluation.strip.plot` <- function(
    data, TrainData=NULL,
    variable.focal=NULL, model.focal=NULL, 
    dev.new.width=7, dev.new.height=7, ...
) 
{
    if (is.null(TrainData) == F) {
        TrainData <- TrainData[TrainData[, "pb"]==1, ]
        TrainData[, "pb"] <- as.factor(TrainData[, "pb"])
    }

    modelnames <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", 
        "MGCVFIX", "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "GLMNET", 
        "BIOCLIM.O", "BIOCLIM", "DOMAIN", "MAHAL", "MAHAL01", "ENSEMBLE")

    modelnames <- names(data)[which(names(data) %in% modelnames)]

    if(is.null(variable.focal)==F) {
        v <- which(names(data) == variable.focal)
        v <- v-2
        f <- data[,1]==v
        vars <- max(data[,1])
# plot for all model.focals
        if (is.null(model.focal) == T) {
            n.models <- length(modelnames)
# model.focals with data
            dim1 <- max(1, ceiling(sqrt(n.models)))
            dim2 <- max(1, ceiling(n.models/dim1))
            par.old <- graphics::par(no.readonly=T)
            if (dev.new.width > 0 && dev.new.height > 0) {grDevices::dev.new(width=dev.new.width, height=dev.new.height)}
            graphics::par(mfrow=c(dim1,dim2))
            for (j in 1:n.models) {
                 if (is.null(TrainData)==T  || is.factor(TrainData[, which(names(TrainData) == variable.focal)])==T) {
                    graphics::plot(data[f,v+2], data[f, 2+vars+j], main=variable.focal, xlab="", ylab=names(data)[2+vars+j], ylim=c(0, 1.25), ...)
                }else{
                    graphics::plot(data[f,v+2], data[f, 2+vars+j], main=variable.focal, xlab="", ylab=names(data)[2+vars+j], ylim=c(0, 1.25), ...)
                    graphics::boxplot(TrainData[, which(names(TrainData) == variable.focal)] ~ TrainData[,"pb"], add=T, horizontal=T)
                }
            }
            graphics::par(par.old)
        }else{
            m <- which(names(data) == model.focal)
            if (dev.new.width > 0 && dev.new.height > 0) {grDevices::dev.new(width=dev.new.width, height=dev.new.height)}
            if (is.null(TrainData)==T  || is.factor(TrainData[, which(names(TrainData) == variable.focal)])==T) {
                graphics::plot(data[f,v+2], data[f, m], main=variable.focal, xlab="", ylab=names(data)[2+vars+j], ylim=c(0, 1.25), ...)
            }else{
                graphics::plot(data[f,v+2], data[f, m], main=variable.focal, xlab="", ylab=names(data)[2+vars+j], ylim=c(0, 1.25), ...)
                graphics::boxplot(TrainData[, which(names(TrainData) == variable.focal)] ~ TrainData[,"pb"], add=T, horizontal=T)
            }
        }
    }
    if(is.null(model.focal)==F && is.null(variable.focal)==T) {
        m <- which(names(data) == model.focal)
# model.focals with data
        vars <- max(data[,1])
        dim1 <- max(1, ceiling(sqrt(vars)))
        dim2 <- max(1, ceiling(vars/dim1))
        if (dev.new.width > 0 && dev.new.height > 0) {grDevices::dev.new(width=dev.new.width, height=dev.new.height)}
        par.old <- graphics::par(no.readonly=T)
        graphics::par(mfrow=c(dim1,dim2))
        for (i in 1:vars) {
            f <- which(data[,1]==i)
            if (is.null(TrainData)==T  || is.factor(TrainData[, which(names(TrainData) == names(data)[i+2])])==T) {
                graphics::plot(data[f,i+2], data[f, m], main=names(data)[i+2], xlab="", ylab=model.focal, ylim=c(0, 1.25), ...)
            }else{
                graphics::plot(data[f,i+2], data[f, m], main=names(data)[i+2], xlab="", ylab=model.focal, ylim=c(0, 1.25), ...)
                varfocal <- names(data)[i+2]
                graphics::boxplot(TrainData[, which(names(TrainData) == varfocal)] ~ TrainData[,"pb"], add=T, horizontal=T)
            }
        }
        graphics::par(par.old)
    }
}

