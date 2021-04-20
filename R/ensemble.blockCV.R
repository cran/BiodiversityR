`ensemble.spatialBlock` <- function(
    x=NULL, p=NULL, 
    a=NULL, an=1000, EPSG=NULL,
    excludep=FALSE, target.groups=FALSE, k=4, 
    factors=NULL,
    theRange=NULL, return.object=FALSE, ...
)
{
# Function to assign presence and background data to spatially separated folds via blockCV::spatialBlock

    ensemble.data <- ensemble.calibrate.models(x=x, p=p, a=a, an=an, 
        SSB.reduce=FALSE, 
        excludep=excludep, target.groups=target.groups, k=0, 
        ENSEMBLE.tune=F,
        MAXENT=0, MAXNET=0, MAXLIKE=0, GBM=0, GBMSTEP=0, RF=0, CF=0, 
        GLM=0, GLMSTEP=0, GAM=0, GAMSTEP=0, MGCV=0, MGCVFIX=0, 
        EARTH=0, RPART=0, NNET=0, FDA=0, SVM=0, SVME=0, GLMNET=0,
        BIOCLIM.O=0, BIOCLIM=0, DOMAIN=0, MAHAL=0, MAHAL01=0,
        factors=factors,
        evaluations.keep=TRUE)

    p.new <- as.data.frame(ensemble.data$evaluations$p)
    a.new <- as.data.frame(ensemble.data$evaluations$a)
    names(a.new) <- names(p.new)

    PA.input <- data.frame(pb=c(rep(1, nrow(p.new)), rep(0, nrow(a.new))), rbind(p.new, a.new))
#    PA.Spatial <- sp::SpatialPointsDataFrame(PA.input[, c(2:3)], data=PA.input, proj4string=raster::crs(x))
    PA.Spatial <- sf::st_as_sf(PA.input, coords=names(a.new), crs=raster::crs(x))
    if (is.null(EPSG) == FALSE) {sf::st_crs(PA.Spatial) <- EPSG}

    sb1 <- blockCV::spatialBlock(speciesData=PA.Spatial, species="pb", theRange=theRange, k=k, ...)

    k <- list(p=p.new, a=a.new, groupp=sb1$foldID[PA.input$pb == 1], groupa=sb1$foldID[PA.input$pb == 0])

    if (return.object == F) {
        return(k)
    }else{
        results <- list(k=k, block.object=sb1, speciesData=PA.Spatial)
        return(results)
    }
}

`ensemble.envBlock` <- function(
    x=NULL, p=NULL, 
    a=NULL, an=1000, EPSG=NULL,
    excludep=FALSE, target.groups=FALSE, k=4, 
    factors=NULL,
    return.object=FALSE, ...
)
{
# Function to assign presence and background data to spatially separated folds via blockCV::envBlock

    ensemble.data <- ensemble.calibrate.models(x=x, p=p, a=a, an=an, 
        SSB.reduce=FALSE, 
        excludep=excludep, target.groups=target.groups, k=0, 
        ENSEMBLE.tune=F,
        MAXENT=0, MAXNET=0, MAXLIKE=0, GBM=0, GBMSTEP=0, RF=0, CF=0, 
        GLM=0, GLMSTEP=0, GAM=0, GAMSTEP=0, MGCV=0, MGCVFIX=0, 
        EARTH=0, RPART=0, NNET=0, FDA=0, SVM=0, SVME=0, GLMNET=0,
        BIOCLIM.O=0, BIOCLIM=0, DOMAIN=0, MAHAL=0, MAHAL01=0,
        factors=factors,
        evaluations.keep=TRUE)

    p.new <- as.data.frame(ensemble.data$evaluations$p)
    a.new <- as.data.frame(ensemble.data$evaluations$a)
    names(a.new) <- names(p.new)

    PA.input <- data.frame(pb=c(rep(1, nrow(p.new)), rep(0, nrow(a.new))), rbind(p.new, a.new))
    PA.Spatial <- sp::SpatialPointsDataFrame(PA.input[, c(2:3)], data=PA.input, proj4string=raster::crs(x))
    if (is.null(EPSG) == FALSE) {sf::st_crs(PA.Spatial) <- EPSG}

    eb1 <- blockCV::envBlock(rasterLayer=x, speciesData=PA.Spatial, species="pb", k=k, ...)
    k <- list(p=p.new, a=a.new, groupp=eb1$foldID[PA.input$pb == 1], groupa=eb1$foldID[PA.input$pb == 0])

    if (return.object == F) {
        return(k)
    }else{
        results <- list(k=k, block.object=eb1, speciesData=PA.Spatial)
        return(results)
    }
}

