
`ensemble.PET.season` <- function(
    PREC.stack=NULL, PET.stack=NULL, 
    filename=NULL, overwrite=TRUE,
    CATCH.OFF=FALSE, ...
)
{
    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if(inherits(PREC.stack, "RasterStack") == FALSE && inherits(PREC.stack, "SpatRaster") == FALSE) {stop("PREC.stack is not a RasterStack or SpatRaster object")}
    if(inherits(PET.stack, "RasterStack") == FALSE && inherits(PET.stack, "SpatRaster") == FALSE) {stop("PET.stack is not a RasterStack or SpatRaster object")}

    names(PREC.stack) <- paste0("PREC", 1:length(names(PREC.stack)))
    names(PET.stack) <- paste0("PET", 1:length(names(PET.stack)))    

    if (inherits(PREC.stack, "RasterStack")) {
        x <- raster::stack(c(PREC.stack, PET.stack))
    }else{
        x <- terra::rast(list(PREC.stack, PET.stack))
    }

    PET.season.object <- list(PREC.names=names(PREC.stack), PET.names=names(PET.stack))
 
    predict.PET.season <- function(object=PET.season.object, newdata=newdata) {
        PREC.names <- object$PREC.names
        PET.names <- object$PET.names        
        result <- array(nrow(newdata))
        for (i in 1:nrow(newdata)) {
            datai <- newdata[i,,drop=F]
            datai.PREC <- datai[, PREC.names]
            datai.PET <- datai[, PET.names]
            datai.BAL <- datai.PET - datai.PREC
            datai.DRY <- datai.BAL > 0
            period.DRY <- rep(0, length(datai.DRY))
            count.period <- 1
            for (j in 1:length(period.DRY)) {
                if (datai.DRY[j] == 1) {
                    period.DRY[j] <- count.period
                }else{
                    count.period <- count.period+1
                }
            }
            if (datai.DRY[1] == TRUE && datai.DRY[length(period.DRY)] == TRUE) {
                old.period <- period.DRY[length(period.DRY)]
                period.DRY[period.DRY == old.period] <- 1
            }
            unique.periods <- c(unique(period.DRY[period.DRY !=0]))
            bal.max <- 0
            if (length(unique.periods) > 0) {
                for (j in length(unique.periods)) {
                    bal.period <- sum(datai.BAL[period.DRY == unique.periods[j]])
                    if (bal.period > bal.max) {bal.max <- bal.period}
                }
            }
            result[i] <- -1 * bal.max        
        }
        return(result)
    }
#
# predict
    
    if (inherits(PREC.stack, "RasterStack")) {
    
        if (CATCH.OFF == F) {
        tryCatch(PET.season.raster <- raster::predict(object=x, model=PET.season.object, fun=predict.PET.season, na.rm=TRUE, 
               filename=filename, overwrite=overwrite, ...),
                   error= function(err) {print(paste("prediction of aridity deficit failed"))},
                   silent=F)
        }else{
            PET.season.raster <- raster::predict(object=x, model=PET.season.object, fun=predict.PET.season, na.rm=TRUE, 
               filename=filename, overwrite=overwrite, ...)
        }

    }else{        
        if (CATCH.OFF == F) {
            tryCatch(PET.season.raster <- terra::predict(object=x, model=PET.season.object, fun=predict.PET.season, na.rm=TRUE, 
                                                          filename=filename, overwrite=overwrite, ...),
                     error= function(err) {print(paste("prediction of aridity deficit failed"))},
                     silent=F)
        }else{
            PET.season.raster <- terra::predict(object=x, model=PET.season.object, fun=predict.PET.season, na.rm=TRUE, 
                                                 filename=filename, overwrite=overwrite, ...)
        }       
          
    }
  
    names(PET.season.raster) <- "PET.season"  
    return(PET.season.raster)  
}

