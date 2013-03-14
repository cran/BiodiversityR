`evaluation.strip.newdata` <- function(
    x, ext=NULL, factors=NULL, steps=50, 
    modelnames=c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL"),
    xn=x
)
{
    data1 <- evaluation.strip.data(x=x, ext=ext, factors=factors, steps=steps, modelnames=modelnames)
    data2 <- evaluation.strip.data(x=xn, ext=ext, factors=factors, steps=steps, modelnames=modelnames)
    data3 <- rbind(data1, data2)
    column1 <- array(dim=c(nrow(data3), 1), 1)
    colnames(column1) <- "newdata"
    data3 <- cbind(column1, data3)
    data3[1:nrow(data1), 1] <- rep(0, nrow(data1))
    return(data3)
}






