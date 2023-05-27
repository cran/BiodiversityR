`ensemble.envirem.masterstack` <- function(
    x,
    precipstack, 
    tmaxstack, tminstack, 
    tmeanstack=NULL
)
{
    if (inherits(precipstack, "RasterBrick")) {precipstack <- raster::stack(precipstack)}
    if (inherits(precipstack, "RasterStack")) {
        precip.data <- data.frame(raster::extract(precipstack, y=x))
    }else if (inherits(precipstack, "SpatRaster")){
        precip.data <- data.frame(terra::extract(precipstack, y=x))
    }else{  
        stop("precipstack is not a RasterStack or SpatRaster object")
    }
    
    names(precip.data) <- paste0("precip_", 1:ncol(precip.data))

    if (inherits(tmaxstack, "RasterBrick")) {tmaxstack <- raster::stack(tmaxstack)}    
    if (inherits(tmaxstack, "RasterStack")) {
        tmax.data <- data.frame(raster::extract(tmaxstack, y=x))
    }else if (inherits(tmaxstack, "SpatRaster")){
        tmax.data <- data.frame(terra::extract(tmaxstack, y=x))
    }else{  
        stop("tmaxstack is not a RasterStack or SpatRaster object")
    }
    
    names(tmax.data) <- paste0("tmax_", 1:ncol(tmax.data))

    if (inherits(tminstack, "RasterBrick")) {tminstack <- raster::stack(tminstack)}        
    if (inherits(tminstack, "RasterStack")) {
        tmin.data <- data.frame(raster::extract(tminstack, y=x))
    }else if (inherits(tminstack, "SpatRaster")){
        tmin.data <- data.frame(terra::extract(tminstack, y=x))
    }else{  
        stop("tminstack is not a RasterStack or SpatRaster object")
    }  
    
    names(tmin.data) <- paste0("tmin_", 1:ncol(tmin.data))   
    
    if (is.null(tmeanstack) == FALSE) {
 
        if (inherits(tmeanstack, "RasterBrick")) {tmeanstack <- raster::stack(tmeanstack)}        
        if (inherits(tmeanstack, "RasterStack")) {
            tmean.data <- data.frame(raster::extract(tmeanstack, y=x))
        }else if (inherits(tmeanstack, "SpatRaster")){
            tmean.data <- data.frame(terra::extract(tmeanstack, y=x))
        }else{  
            stop("tmeanstack is not a RasterStack or SpatRaster object")
        }  
        names(tmean.data) <- paste0("tmean_", 1:ncol(tmean.data))   
        input.data <- cbind(precip.data, tmax.data, tmin.data, tmean.data)
        
    }else{
        input.data <- cbind(precip.data, tmax.data, tmin.data, tmean.data)
    }
    
    for (i in 1:ncol(input.data)) {
        rasteri <- raster::raster(matrix(input.data[, i]))
        if (i == 1) {
            masterstack <- raster::stack(rasteri)
        }else{
            masterstack <- raster::stack(c(masterstack, rasteri))
        }
    }
    
    names(masterstack) <- names(input.data)
    return(masterstack)
}

`ensemble.envirem.solradstack` <- function(
    x, solrad
)
{
    if (inherits(solrad, "RasterBrick")) {solrad <- raster::stack(solrad)}
    if (inherits(solrad, "RasterStack")) {
        input.data <- data.frame(raster::extract(solrad, y=x))
    }else if (inherits(solrad, "SpatRaster")){
        input.data <- data.frame(terra::extract(solrad, y=x))
    }else{  
        stop("solrad is not a RasterStack or SpatRaster object")
    }
    
    names(input.data) <- paste0("et_solrad_", 1:ncol(input.data))
    
    for (i in 1:ncol(input.data)) {
        rasteri <- raster::raster(matrix(input.data[, i]))
        if (i == 1) {
            solradout <- raster::stack(rasteri)
        }else{
            solradout <- raster::stack(c(solradout, rasteri))
        }
    }
    
    names(solradout) <- names(input.data)
    return(solradout)
    
}

`ensemble.envirem.run` <- function(
    masterstack, solradstack,
    var="all", ...
)
{
    envirem.out <- envirem::layerCreation(masterstack=masterstack,
                                          solradstack=solradstack,
                                          var=var, ...)
# Modified May 2023 with optional argument  
    return(as.data.frame(envirem.out, optional=FALSE))
}


