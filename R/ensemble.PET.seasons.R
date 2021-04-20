
`ensemble.PET.seasons` <- function(
    PREC.stack=NULL, PET.stack=NULL, 
    index=c("seasons", "start1", "length1", "start2", "length2", "start3", "length3"),
    filename=NULL, overwrite=TRUE,
    CATCH.OFF=FALSE, ...
)
{
    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if(inherits(PREC.stack, "RasterStack") == F) {stop("PREC.stack is not a RasterStack object")}
    if(inherits(PET.stack, "RasterStack") == F) {stop("PET.stack is not a RasterStack object")}
    
    names(PREC.stack) <- paste0("PREC", 1:length(names(PREC.stack)))
    names(PET.stack) <- paste0("PET", 1:length(names(PET.stack)))    

    x <- raster::stack(c(PREC.stack, PET.stack)) 

    PET.season.object <- list(PREC.names=names(PREC.stack), PET.names=names(PET.stack))

    indices <- c("seasons", "start1", "length1", "start2", "length2", "start3", "length3")
    index1 <- indices[match(index, indices)]
 
    predict.PET.seasons <- function(object=PET.season.object, newdata=newdata, index=index1) {
        PREC.names <- object$PREC.names
        PET.names <- object$PET.names        
        result <- array(nrow(newdata))
        for (i in 1:nrow(newdata)) {
            datai <- newdata[i,,drop=F]
            datai.PREC <- datai[, PREC.names]
            datai.PET <- datai[, PET.names]
            datai.BAL <- datai.PET - 2 * datai.PREC
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
            unique.periods <- sort(unique(period.DRY))
            unique.periods <- unique.periods[unique.periods != 0]
            for (j in 1:length(unique.periods)) {
                period.DRY[period.DRY == unique.periods[j]] <- -1*j
            }
            if (period.DRY[1] == 0) {period.DRY[1] <- abs(min(period.DRY))}
            if (period.DRY[1] == 0) {period.DRY[1] <- 1}
            for (j in 2:length(period.DRY)) {
                if (period.DRY[j] == 0) {period.DRY[j] <- abs(period.DRY[j-1])}
            }
            names(result) <- paste0("PE", c(1:length(result)))          
            if (i == 1) {
                result <- period.DRY
            }else{
                result <- rbind(result, period.DRY)
            }
        }
        result2 <- data.frame(result)
        result2$seasons <- apply(result, FUN="max", MARGIN=1)
        result2[result2$seasons == -1, "seasons"] <- 0
        result2$start3 <- result2$start2 <- result2$start1 <- NA
        result2$length3 <- result2$length2 <- result2$length1 <- NA
        for (i in 1:nrow(result2)) {
            resulti <- result[i, ]
            resulti.pos <- resulti
            resulti.pos[resulti.pos < 0] <- NA
            if (result2[i, "seasons"] !=0) {
                resulti.pos1 <- resulti.pos
                resulti.pos1[resulti.pos1 > 1] <- NA
                s1 <- which.min(resulti.pos1 == 1)
                if (s1 == 1) {
                    resulti.neg <- resulti
                    resulti.neg[resulti.neg > 0] <- 0      
                    resulti.neg[resulti.neg < 0] <- 1                    
                    resulti.neg <- data.frame(t(resulti.neg))              
                    e1 <- max(0, max.col(resulti.neg, ties.method="last"), na.rm=T)
                    if (e1 < length(resulti.neg)) {s1 <- e1+1}     
                }              
                result2[i, "start1"] <- s1
                result2[i, "length1"] <- sum(resulti.pos1 == 1, na.rm=T)
            }
            if (result2[i, "seasons"] > 1) {
                resulti.pos2 <- resulti.pos
                resulti.pos2[resulti.pos2 != 2] <- NA
                s2 <- which.min(resulti.pos2 == 2)         
                result2[i, "start2"] <- s2
                result2[i, "length2"] <- sum(resulti.pos2 == 2, na.rm=T)
            }
            if (result2[i, "seasons"] > 2) {
                resulti.pos3 <- resulti.pos
                resulti.pos3[resulti.pos3 != 3] <- NA
                s3 <- which.min(resulti.pos3 == 3)         
                result2[i, "start3"] <- s3
                result2[i, "length3"] <- sum(resulti.pos3 == 3, na.rm=T)
            }       
        }
        rownames(result2) <- NULL
        result3 <- result2[, index]
        return(result3)
    }
#
# predict
    if (CATCH.OFF == F) {
        tryCatch(PET.seasons.raster <- raster::predict(object=x, model=PET.season.object, fun=predict.PET.seasons, index=index1, na.rm=TRUE, 
               filename=filename, overwrite=overwrite, ...),
           error= function(err) {print(paste("prediction of", index1, "failed"))},
           silent=F)
    }else{
        PET.seasons.raster <- raster::predict(object=x, model=PET.season.object, fun=predict.PET.seasons, na.rm=TRUE, 
               filename=filename, overwrite=overwrite, ...)
    }
#    
    return(PET.seasons.raster)  
}


`ensemble.prec.season` <- function(
    PREC.stack=NULL, start.layer=NULL, length.layer=NULL,
    filename=NULL, overwrite=TRUE,
    CATCH.OFF=FALSE, ...
)
{
    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if(inherits(PREC.stack, "RasterStack") == F) {stop("PREC.stack is not a RasterStack object")}
    if(inherits(start.layer, "RasterLayer") == F) {stop("start.layer is not a RasterLayer object")}
    if(inherits(length.layer, "RasterLayer") == F) {stop("length.layer is not a RasterLayer object")}
    
    names(PREC.stack) <- paste0("PREC", 1:length(names(PREC.stack)))
    names(start.layer) <- "start"    
    names(length.layer) <- "length"  

    x <- raster::stack(c(PREC.stack, start.layer, length.layer)) 

    prec.season.object <- list(PREC.names=names(PREC.stack))
 
    predict.prec.season <- function(object=prec.season.object, newdata=newdata) {
        PREC.names <- object$PREC.names
        n.mts <- length(PREC.names)       
        result <- array(nrow(newdata))
        for (i in 1:nrow(newdata)) {
            datai <- newdata[i, , drop=F]
            datai.PREC <- datai[, PREC.names]
            mts <- seq(from=datai[, "start"], length=datai[, "length"])
            mts[mts > n.mts] <- mts[mts > n.mts] - n.mts
            result[i] <- sum(datai.PREC[mts])                    
        }
        return(result)
    }
#
# predict
    if (CATCH.OFF == F) {
        tryCatch(prec.season.raster <- raster::predict(object=x, model=prec.season.object, fun=predict.prec.season, na.rm=TRUE, 
               filename=filename, overwrite=overwrite, ...),
           error= function(err) {print(paste("prediction failed"))},
           silent=F)
    }else{
        prec.season.raster <- raster::predict(object=x, model=prec.season.object, fun=predict.prec.season, na.rm=TRUE, 
               filename=filename, overwrite=overwrite, ...)
    }
#    
    return(prec.season.raster)  
}


`ensemble.tmean.season` <- function(
    TMEAN.stack=NULL, start.layer=NULL, length.layer=NULL,
    filename=NULL, overwrite=TRUE,
    CATCH.OFF=FALSE, ...
)
{
    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if(inherits(TMEAN.stack, "RasterStack") == F) {stop("TMEAN.stack is not a RasterStack object")}
    if(inherits(start.layer, "RasterLayer") == F) {stop("start.layer is not a RasterLayer object")}
    if(inherits(length.layer, "RasterLayer") == F) {stop("length.layer is not a RasterLayer object")}
    
    names(TMEAN.stack) <- paste0("TMEAN", 1:length(names(TMEAN.stack)))
    names(start.layer) <- "start"    
    names(length.layer) <- "length"  

    x <- raster::stack(c(TMEAN.stack, start.layer, length.layer)) 

    tmean.season.object <- list(TMEAN.names=names(TMEAN.stack))
 
    predict.tmean.season <- function(object=tmean.season.object, newdata=newdata) {
        tmean.names <- object$TMEAN.names
        n.mts <- length(tmean.names)       
        result <- array(nrow(newdata))
        for (i in 1:nrow(newdata)) {
            datai <- newdata[i, , drop=F]
            datai.TMEAN <- datai[, tmean.names]
            mts <- seq(from=datai[, "start"], length=datai[, "length"])
            mts[mts > n.mts] <- mts[mts > n.mts] - n.mts
            result[i] <- sum(datai.TMEAN[mts])/length(mts)
        }
        return(result)
    }
#
# predict
    if (CATCH.OFF == F) {
        tryCatch(tmean.season.raster <- raster::predict(object=x, model=tmean.season.object, fun=predict.tmean.season, na.rm=TRUE, 
               filename=filename, overwrite=overwrite, ...),
           error= function(err) {print(paste("prediction failed"))},
           silent=F)
    }else{
        tmean.season.raster <- raster::predict(object=x, model=tmean.season.object, fun=predict.tmean.season, na.rm=TRUE, 
               filename=filename, overwrite=overwrite, ...)
    }
#    
    return(tmean.season.raster)  
}


`ensemble.season.suitability` <- function(
    season.raster=NULL, thresholds=NULL,
    filename=NULL, overwrite=TRUE,
    CATCH.OFF=FALSE, ...
)
{
    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if(inherits(season.raster, "RasterLayer") == F) {stop("season.raster is not a RasterLayer object")}
    
    suit.object <- list(thresholds=thresholds[order(thresholds)])
 
    predict.suit <- function(suit.object=suit.object, newdata=newdata) {
        VMIN <- as.numeric(suit.object$thresholds[1])
        VOPMN <- as.numeric(suit.object$thresholds[2])
        VOPMX <- as.numeric(suit.object$thresholds[3])
        VMAX <- as.numeric(suit.object$thresholds[4])     
        result <- array(nrow(newdata))
        for (i in 1:nrow(newdata)) {
            datai <- newdata[i, , drop=F]
            R.out <- 1
            if (datai <= VMIN) {R.out <- 0}
            if (datai > VMIN && datai < VOPMN) {
                R.out <- 1 - ((VOPMN - datai) / (VOPMN - VMIN))
            }
            if (datai >= VMAX) {R.out <- 0}  
            if (datai > VOPMX && datai < VMAX) {
                R.out <- 1 - ((datai - VOPMX) / (VMAX - VOPMX))
            }
            result[i] <- R.out                   
        }
        return(result)
    }
#
# predict
    if (CATCH.OFF == F) {
        tryCatch(suit.raster <- raster::predict(object=season.raster, model=suit.object, fun=predict.suit, na.rm=TRUE, 
               filename=filename, overwrite=overwrite, ...),
           error= function(err) {print(paste("prediction failed"))},
           silent=F)
    }else{
        suit.raster <- raster::predict(object=season.raster, model=suit.object, fun=predict.suit, na.rm=TRUE, 
               filename=filename, overwrite=overwrite, ...)
    }
#    
    return(suit.raster)  
}



