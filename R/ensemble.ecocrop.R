`ensemble.ecocrop.object` <- function(
    temp.thresholds, rain.thresholds, name="crop01", 
    temp.multiply=1, annual.temps=TRUE, transform=1
)
{
    temps <- as.numeric(temp.thresholds[order(temp.thresholds)])
    temps <- temps*temp.multiply
    names(temps) <- c("tminabs", "tminopt", "tmaxopt", "tmaxabs")
    rains <- rain.thresholds[order(rain.thresholds)]
    names(rains) <- c("pminabs", "pminopt", "pmaxopt", "pmaxabs")
    ecocrop.object <- list(name=name, temp.thresholds=temps, rain.thresholds=rains, annual.temps=annual.temps, transform=transform)
    return(ecocrop.object)
}


`ensemble.ecocrop` <- function(
    x=NULL, ecocrop.object=NULL, 
    RASTER.object.name=ecocrop.object$name, 
    RASTER.stack.name = "xTitle", RASTER.format = "GTiff", 
    RASTER.datatype = "INT2S", RASTER.NAflag = -32767, 
    CATCH.OFF = FALSE
) 
{
  .BiodiversityR <- new.env()
  if (is.null(x) == T) {
    stop("value for parameter x is missing (RasterStack or SpatRaster object)")
  }
  if (inherits(x, "RasterStack") == F  && inherits(x, "SpatRaster") == FALSE) {
    stop("x is not a RasterStack or SpatRaster object")
  }
  names(x)[which(names(x) == "bio01")] <- "bio1"
  names(x)[which(names(x) == "bio05")] <- "bio5"
  names(x)[which(names(x) == "bio06")] <- "bio6"
  if (is.null(ecocrop.object) == TRUE) {
    stop("value for parameter ecocrop.object is missing (hint: use the ensemble.ecocrop.object function)")
  }
  vars <- names(x)
  if (any(vars == "bio12") == F) {
    stop("Bioclimatic variable 'bio12' not provided with data")
  }
  if (ecocrop.object$annual.temps == F) {
    if (any(vars == "bio5") == F) {
      stop("Bioclimatic variable 'bio5' not provided with data")
    }
    if (any(vars == "bio6") == F) {
      stop("Bioclimatic variable 'bio6' not provided with data")
    }
  }
  else {
    if (any(vars == "bio1") == F) {
      stop("Bioclimatic variable 'bio1' not provided with data")
    }
  }
  predict.ecocrop <- function(object = ecocrop.object, newdata = newdata) {
    tminopt <- object$temp.thresholds["tminopt"]
    tminabs <- object$temp.thresholds["tminabs"]
    tmaxopt <- object$temp.thresholds["tmaxopt"]
    tmaxabs <- object$temp.thresholds["tmaxabs"]
    pminopt <- object$rain.thresholds["pminopt"]
    pminabs <- object$rain.thresholds["pminabs"]
    pmaxopt <- object$rain.thresholds["pmaxopt"]
    pmaxabs <- object$rain.thresholds["pmaxabs"]
    annual.temps <- object$annual.temps
    z <- object$transform
    # modified from array to numeric
    result <- numeric(nrow(newdata))
    
    Pall <- as.numeric(newdata[, "bio12"])     
    if (annual.temps == F) {
      TMIall <- as.numeric(newdata[, "bio6"])
      TMXall <- as.numeric(newdata[, "bio5"])
    }else{
      TMIall <- as.numeric(newdata[, "bio1"])
      TMXall <- TMIall
    }
    
    for (i in 1:length(Pall)) {
      #            datai <- newdata[i, , drop = F]
      P <- Pall[i]
      TMI <- TMIall[i]
      TMX <- TMXall[i]
      
      if (is.na(P)==FALSE && is.na(TMI)==FALSE && is.na(TMX)==FALSE) {
        
        PS1 <- PS2 <- PS3 <- 0
        if (P > pminabs && P < pminopt) {
          PS1 <- (P - pminabs)/(pminopt - pminabs)
        }
        if (P >= pminopt && P <= pmaxopt) {
          PS2 <- 1
        }
        if (P > pmaxopt && P < pmaxabs) {
          PS3 <- (pmaxabs - P)/(pmaxabs - pmaxopt)
        }
        PS <- max(c(PS1, PS2, PS3))
        
        TMI1 <- TMI2 <- 0
        if (TMI >= tminopt) {
          TMI1 <- 1
        }
        if (TMI > tminabs && TMI < tminopt) {
          TMI2 <- (TMI - tminabs)/(tminopt - tminabs)
        }
        TMIS <- max(c(TMI1, TMI2))
        
        
        TMX1 <- TMX2 <- 0
        if (TMX <= tmaxopt) {
          TMX1 <- 1
        }
        if (TMX > tmaxopt && TMX < tmaxabs) {
          TMX2 <- (tmaxabs - TMX)/(tmaxabs - tmaxopt)
        }
        TMXS <- max(c(TMX1, TMX2))
        
        SFINAL <- min(PS, TMIS, TMXS)
        result[i] <- SFINAL
      }else{
        result[i] <- NA
      }
    }  
    #        p <- as.numeric(result)
    result <- result^z
    result <- trunc(1000*result)
    return(result)
  }
  dir.create("ensembles", showWarnings = F)
  dir.create("ensembles/ecocrop", showWarnings = F)
  
  stack.title <- RASTER.stack.name
  rasterfull <- paste("ensembles/ecocrop/", RASTER.object.name, 
                      "_", stack.title, sep = "")
  
  if (inherits(x, "RasterStack") == TRUE) {
    
    if (CATCH.OFF == F) {
      tryCatch(ecocrop.raster <- raster::predict(object = x, 
                                                 model = ecocrop.object, fun = predict.ecocrop, na.rm = TRUE, 
                                                 filename = rasterfull, progress = "text", overwrite = TRUE, 
                                                 format = RASTER.format), error = function(err) {
                                                   print(paste("ecocrop prediction failed"))
                                                 }, silent = F)
    }else{
      ecocrop.raster <- raster::predict(object = x, model = ecocrop.object, 
                                        fun = predict.ecocrop, na.rm = TRUE, filename = rasterfull, 
                                        progress = "text", overwrite = TRUE, format = RASTER.format)
    }
    
  } # RasterStack
  
  if (inherits(x, "SpatRaster") == TRUE) {
    
    rasterfull <- paste0(rasterfull, ".tif") 
    
    if (CATCH.OFF == F) {
      tryCatch(ecocrop.raster <- terra::predict(object = x, 
                                                model = ecocrop.object, fun = predict.ecocrop, na.rm = TRUE, 
                                                filename = rasterfull,
                                                overwrite = TRUE), error = function(err) {
                                                  print(paste("ecocrop prediction failed"))
                                                }, silent = F)
    }else{
      ecocrop.raster <- terra::predict(object = x, model = ecocrop.object, 
                                       fun = predict.ecocrop, na.rm = TRUE, filename = rasterfull, 
                                       overwrite = TRUE)
    }
    
  } # SpatRaster   
  
  #    ecocrop.raster <- trunc(1000 * ecocrop.raster)
  cat(paste("\n", "raster layer created (probabilities multiplied by 1000)", 
            "\n", sep = ""))
  #    raster::setMinMax(ecocrop.raster)
  print(ecocrop.raster)
  #    raster::writeRaster(ecocrop.raster, filename = "working.grd", 
  #        overwrite = T)
  #    working.raster <- raster::raster("working.grd")
  #    names(working.raster) <- paste(RASTER.object.name, "_", 
  #        stack.title, "_ecocrop", sep = "")
  #    raster::writeRaster(ecocrop.raster, filename = rasterfull, 
  #        progress = "text", overwrite = TRUE, format = RASTER.format, 
  #        datatype = RASTER.datatype, NAflag = RASTER.NAflag)
  
  cat(paste("\n", "ecocrop raster provided in folder: ", 
            getwd(), "//ensembles//ecocrop", "\n", sep = ""))
  return(ecocrop.raster)
}
