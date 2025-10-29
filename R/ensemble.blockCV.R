`ensemble.spatialBlock` <- function(
    x=NULL, p=NULL, 
    a=NULL, an=1000,
    EPSG=NULL,
    excludep=FALSE, target.groups=FALSE, k=4,
    factors=NULL,
    size=NULL, return.object=FALSE, ...
)
{
# Function to assign presence and background data to spatially separated folds via blockCV's spatialBlock

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
    PA.Spatial <- sf::st_as_sf(PA.input, coords=names(a.new), crs=raster::crs(x))
    if (is.null(EPSG) == FALSE) {sf::st_crs(PA.Spatial) <- EPSG}

# Disabled August 2025 when blockCV was archived - reenabled October 2025
# modified October 2025 to use cv_spatial instead of blockCV::spatialBlock
#    sb1 <- blockCV::spatialBlock(speciesData=PA.Spatial, species="pb", theRange=theRange, k=k, ...)
    sb1 <- blockCV::cv_spatial(x=PA.Spatial, column="pb", size=size, k=k, ...)
    
# October 2025: remove locations not assigned to any fold
    PA.output <- PA.input[is.na(sb1$folds_ids) == FALSE, ]
    folds.output <- sb1$folds_ids[is.na(sb1$folds_ids) == FALSE]

    k <- list(p=PA.output[PA.output$pb == 1, -1], 
              a=PA.output[PA.output$pb == 0, -1], 
              groupp=folds.output[PA.output$pb == 1], 
              groupa=folds.output[PA.output$pb == 0])

# Disabled August 2025 - reenabled October 2025
#    k <- sb1 <- NULL
    
    if (return.object == FALSE) {
        return(k)
    }else{
        results <- list(k=k, block.object=sb1, 
                        speciesData=PA.Spatial,
                        speciesData.folded=PA.Spatial[is.na(sb1$folds_ids) == FALSE, ])
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
# Function to assign presence and background data to spatially separated folds via blockCV's envBlock

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

# disabled August 2025 when blockCV was archived - reenabled October 2025
# modified October 2025 to use cv_cluster instead of blockCV::envBlock
#    eb1 <- blockCV::envBlock(rasterLayer=x, speciesData=PA.Spatial, species="pb", k=k, ...)
    eb1 <- blockCV::cv_cluster(r=terra::rast(x), x=PA.Spatial, column="pb", k=k, ...)
    k <- list(p=p.new, a=a.new, groupp=eb1$foldID[PA.input$pb == 1], groupa=eb1$foldID[PA.input$pb == 0])
    
#    k <- eb1 <- NULL

    # October 2025: remove locations not assigned to any fold
    PA.output <- PA.input[is.na(eb1$folds_ids) == FALSE, ]
    folds.output <- eb1$folds_ids[is.na(eb1$folds_ids) == FALSE]
    
    k <- list(p=PA.output[PA.output$pb == 1, -1], 
              a=PA.output[PA.output$pb == 0, -1], 
              groupp=folds.output[PA.output$pb == 1], 
              groupa=folds.output[PA.output$pb == 0])
    
        
    if (return.object == F) {
        return(k)
    }else{
        results <- list(k=k, block.object=eb1, 
                        speciesData=PA.Spatial,
                        speciesData.folded=PA.Spatial[is.na(eb1$folds_ids) == FALSE, ])
        return(results)
    }
}

