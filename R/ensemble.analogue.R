# Function to get information for future climate of the current location
# implemented through a function to generate the main parameters

# modified from the Bioclim methods of Robert J. Hijmans

if (require(raster) == T) {

setClass('Climsurf',
    representation (
        p='data.frame',
	target.values='numeric',
	norm.values='numeric'),	
    prototype (),
    validity = function(object)	{return(TRUE)}
)


if (!isGeneric("climsurface")) {
	setGeneric("climsurface", function(ref.location, future.stack, current.stack, 
            method='quantile', q.limits=c(0.25, 0.75), weights=rep(1, nlayers(current.stack)), ...)
	standardGeneric("climsurface"))
}


# climsurface generates the parameters needed to predict
setMethod('climsurface', signature(ref.location='data.frame', future.stack='RasterStack', current.stack='RasterStack'),
          function(ref.location, future.stack, current.stack, method, q.limits, weights)
          {
            target.values <- as.numeric(extract(future.stack, ref.location))
            names(target.values) <- names(future.stack)
            if (method == "quantile") {
              lower.interval <- data.frame(t(quantile(current.stack, q.limits[1])))
              upper.interval <- data.frame(t(quantile(current.stack, q.limits[2])))
              norm.values <- as.numeric(upper.interval - lower.interval)
              names(norm.values) <- names(current.stack)
            }
            if (method == "sd") {
              norm.values <- cellStats(current.stack, stat="sd")
            }
            if (method == "none") {
              norm.values <- rep(1, length=length(names(current.stack)))
              names(norm.values) <- names(current.stack)
            }
            if (length(weights) != nlayers(future.stack)) {
              paste("WARNING: length of weights is different from number of variables in current RasterStack", "\n", sep="") 
              weight.values <- rep(1, nlayers(current.stack))
            }
            
            # problem if some of the norm values are zero
            zero.norm.values <- which(norm.values == 0)
            if(length(zero.norm.values) > 0) {
              cat(paste("WARNING: some of the normalizing values were zero", "\n\n", sep=""))
              print(names(zero.norm.values))
              cat(paste("\n", "respective values were now set to one", "\n", sep=""))
              norm.values[names(norm.values) %in% names(zero.norm.values)] <- 1
            }
            weight.values <- weights
            if (length(names(weight.values)) == 0) {names(weight.values) <- names(current.stack)}
            weight.values <- weight.values / sum(weight.values)
            
            out <- list(p=ref.location, target.values=target.values, norm.values=norm.values, weight.values=weight.values, stack.name=future.stack@title)
            class(out) <- "Climsurf"
            return(out)
          }
)

# predict function generates the predictions
setMethod('predict', signature(object='Climsurf'), 
    function(object, current.data, z=2) {

	out <- data.frame(current.data)
	targetdata <-object$target.values
	normdata <- object$norm.values
	weightdata <- object$weight.values
	for (i in 1:ncol(out)) {
		out[,i] <- as.numeric(out[,i]) - as.numeric(targetdata[which(names(targetdata) == names(out)[i])])
		out[,i] <- abs(out[,i])
		out[,i] <- as.numeric(out[,i]) / as.numeric(normdata[which(names(normdata) == names(out)[i])])
		out[,i] <- (out[,i]) ^ z
		out[,i] <- as.numeric(out[,i]) * as.numeric(weightdata[which(names(weightdata) == names(out)[i])])
	}
	out2 <- rowSums(out)
	z2 <- 1/z
	out2 <- out2 ^ z2
	return(out2)
    }
)


# climate.analogues produces the grid layers and finds required number of analogues
# analogues gives number of analogue locations required

`climate.analogues` <- function(current.stack, model, analogues=1, z=2,
	filename=paste(row.names(model$p), "_", model$stack.name, sep=""), 
	overwrite=T, KML.out=T, KML.blur=10, 
	limits=c(1, 5, 20, 50), limit.colours=c('red', 'orange', 'blue', 'grey'), ...) 
{
	raster.out <- current.stack[[1]]
	raster.out <- raster::predict(object=current.stack, model=model)
	setMinMax(raster.out)
	names(raster.out) <- paste(row.names(model$p), "_", model$stack.name, sep="")
# file to be saved in analogues subfolder
        dir.create("analogues", showWarnings = F)
	filename <- paste(getwd(), "//analogues//", filename, sep="")
	writeRaster(raster.out, filename=filename, overwrite=overwrite, ...)
	names(raster.out) <- paste(row.names(model$p), "_", model$stack.name, sep="")
	print(raster.out)
#	
	limits.data <- data.frame(limits)
	limits.data <- data.frame(cbind(limits, limits))
	names(limits.data) <- c("threshold", "count")
	freqs <- freq(raster.out, digits=8)
	for (i in 1:length(limits)) {
		j <- 1
        	while (sum(freqs[1:j, 2]) < limits[i]) {j <- j+1}
        	limits.data[i, 1] <- freqs[j, 1]
        	limits.data[i, 2] <- sum(freqs[1:j, 2])
	}
	cat(paste("suggested breaks in colour scheme", "\n", sep=""))	
	print(limits.data)
	if (KML.out == T) {
		breaks1 <- c(minValue(raster.out), limits.data[,1])
		KML(raster.out, filename=filename, overwrite=overwrite, blur=KML.blur, col=limit.colours, breaks=breaks1)
	}
# get locations
	j <- 1
        while (sum(freqs[1:j, 2]) < analogues) {j <- j+1}
        threshold <- freqs[j, 1]
	index1 <- which(raster.out[,] <= threshold)
	pres1 <- xyFromCell(raster.out, index1)
#
	vars <- length(names(current.stack))
	output1 <- data.frame(array(dim=c(length(index1), 4+vars)))
	names(output1) <- c("model", "lon", "lat", "distance", names(current.stack))
	output1[, 1] <- rep(model$stack.name, nrow(output1))
	output1[, c(2:3)] <- pres1
	point.data1 <- extract(current.stack, pres1)
	output1[, c(5:(4+vars))] <- point.data1
	point.data2 <- extract(raster.out, pres1)
	output1[, 4] <- point.data2
	output1 <- output1[order(output1[,"distance"], decreasing=F),]
	return(output1)
}

}
