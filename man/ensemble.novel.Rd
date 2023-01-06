\name{ensemble.novel}
\alias{ensemble.novel}
\alias{ensemble.novel.object}

\title{
Mapping of novel environmental conditions (areas where some of the environmental conditions are outside the range of environmental conditions of a reference area).
}

\description{
Function \code{ensemble.novel} creates the map with novel conditions. Function \code{ensemble.novel.object} provides the reference values used by the prediction function used by \code{\link[raster]{predict}} . 
}

\usage{
ensemble.novel(x = NULL, novel.object = NULL,
    RASTER.object.name = novel.object$name, RASTER.stack.name = x@title,
    RASTER.format = "GTiff", RASTER.datatype = "INT2S", RASTER.NAflag = -32767, 
    CATCH.OFF = FALSE)

ensemble.novel.object(x = NULL, name = "reference1", mask.raster = NULL,
    quantiles = FALSE, probs = c(0.05, 0.95), factors = NULL)

}

\arguments{
  \item{x}{ RasterStack object (\code{\link[raster]{stack}}) containing all environmental layers for which novel conditions should be calculated. With \code{\link{ensemble.novel.object}}, x can also be a data.frame.  }
  \item{novel.object}{ Object listing minima and maxima for the environmental layers, used by the prediction function that is used internally by \code{\link[raster]{predict}}. This object is created with \code{\link{ensemble.novel.object}}. }
  \item{RASTER.object.name}{ First part of the names of the raster file that will be generated, expected to identify the area and time period for which ranges were calculated }
  \item{RASTER.stack.name}{ Last part of the names of the raster file that will be generated, expected to identify the predictor stack used }
  \item{RASTER.format}{ Format of the raster files that will be generated. See \code{\link[raster]{writeFormats}} and \code{\link[raster]{writeRaster}}.  }
  \item{RASTER.datatype}{ Format of the raster files that will be generated. See \code{\link[raster]{dataType}} and \code{\link[raster]{writeRaster}}.  }
  \item{RASTER.NAflag}{ Value that is used to store missing data. See \code{\link[raster]{writeRaster}}.  }
  \item{CATCH.OFF}{ Disable calls to function \code{\link[base]{tryCatch}}.  }

  \item{name}{ Name of the object, expect to expected to identify the area and time period for which ranges were calculated and where no novel conditions will be detected }
  \item{mask.raster}{ RasterLayer object (\code{\link[raster]{raster}}) that can be used to select the area for which reference values are obtained (see \code{\link[raster]{mask}})  }
  \item{quantiles}{ If \code{TRUE}, then replace minima and maxima with quantile values. See also \code{\link[stats]{quantile}} and \code{\link[raster]{quantile}})  }
  \item{probs}{ Numeric vector of probabilities [0, 1] as used by \code{\link[stats]{quantile}} and \code{\link[raster]{quantile}})  }
  \item{factors}{ vector that indicates which variables are factors; these variables will be ignored for novel conditions } 
}

\details{
Function \code{ensemble.novel} maps zones (coded '1') that are novel (outside the minimum-maximum range) relative to the range provided by function \code{ensemble.novel.object}. Values that are not novel (inside the range of minimum-maximum values) are coded '0'. In theory, the maps show the same areas that have negative Multivariate Environmental Similarity Surface (MESS) values ((\code{\link[dismo]{mess}}))
}

\value{
Function \code{ensemble.novel.object} returns a list with following objects:
\item{minima }{ minima of the reference environmental conditions }
\item{maxima }{ maxima of the reference environmental conditions }
\item{name }{ name for the reference area and time period }
}

\author{ Roeland Kindt (World Agroforestry Centre)}

\seealso{
\code{\link{ensemble.raster}}, \code{\link{ensemble.bioclim}} and \code{\link[BiodiversityR]{ensemble.bioclim.graph}} 
}

\examples{

\dontrun{
# get predictor variables
library(dismo)
predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
    pattern='grd', full.names=TRUE)
predictors <- stack(predictor.files)
predictors <- subset(predictors, subset=c("bio1", "bio5", "bio6", "bio7", "bio8", 
    "bio12", "bio16", "bio17"))
predictors
predictors@title <- "base"

# reference area to calculate environmental ranges
ext <- extent(-70, -50, -10, 10)
extent.values2 <- c(-70, -50, -10, 10)
predictors.current <- crop(predictors, y=ext)
predictors.current <- stack(predictors.current)

novel.test <- ensemble.novel.object(predictors.current, name="noveltest")
novel.test
novel.raster <- ensemble.novel(x=predictors, novel.object=novel.test)
novel.raster

plot(novel.raster)
# no novel conditions within reference area
rect(extent.values2[1], extent.values2[3], extent.values2[2], extent.values2[4])

# use novel conditions as a simple species suitability mapping method
# presence points
presence_file <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
pres <- read.table(presence_file, header=TRUE, sep=',')[,-1]
pres.data <- data.frame(extract(predictors, y=pres))

# ranges and maps
Bradypus.ranges1 <- ensemble.novel.object(pres.data, name="Bradypus", quantiles=F)
Bradypus.ranges1
Bradypus.novel1 <- ensemble.novel(x=predictors, novel.object=Bradypus.ranges1)
Bradypus.novel1

par.old <- graphics::par(no.readonly=T)
graphics::par(mfrow=c(1,2))

# suitable where there are no novel conditions
raster::plot(Bradypus.novel1, breaks=c(-0.1, 0, 1), col=c("green", "grey"), 
    main="Suitability mapping using minimum to maximum range")
points(pres[, 2] ~ pres[, 1], pch=1, col="red", cex=0.8)

# use 90 percent intervals similar to BIOCLIM methodology
Bradypus.ranges2 <- ensemble.novel.object(pres.data, name="BradypusQuantiles", quantiles=T)
Bradypus.ranges2
Bradypus.novel2 <- ensemble.novel(x=predictors, novel.object=Bradypus.ranges2)
Bradypus.novel2
raster::plot(Bradypus.novel2, breaks=c(-0.1, 0, 1), col=c("green", "grey"), 
    main="Suitability mapping using quantile range")
points(pres[, 2] ~ pres[, 1], pch=1, col="red", cex=0.8)

graphics::par(par.old)

# deal with novel factor levels through dummy variables
predictors <- stack(predictor.files)
biome.layer <- predictors[["biome"]]
biome.layer
ensemble.dummy.variables(xcat=biome.layer, most.frequent=0, freq.min=1,
    overwrite=TRUE)

predictors.dummy <- stack(predictor.files)
predictors.dummy <- subset(predictors.dummy, subset=c("biome_1", "biome_2",  "biome_3",  
    "biome_4", "biome_5", "biome_7",  "biome_8",  "biome_9", 
    "biome_10", "biome_12", "biome_13", "biome_14"))
predictors.dummy
predictors.dummy@title <- "base_dummy"

predictors.dummy.current <- crop(predictors.dummy, y=ext)
predictors.dummy.current <- stack(predictors.dummy.current)

novel.levels <- ensemble.novel.object(predictors.dummy.current, name="novellevels")
novel.levels
novel.levels.raster <- ensemble.novel(x=predictors.dummy, novel.object=novel.levels)
novel.levels.raster

novel.levels.quantiles <- ensemble.novel.object(predictors.dummy.current, quantiles=TRUE,
    name="novellevels_quantiles")
novel.levels.quantiles
novel.levels.quantiles.raster <- ensemble.novel(x=predictors.dummy, 
    novel.object=novel.levels.quantiles)
novel.levels.quantiles.raster

# difference in ranges for variables with low frequencies
background <- dismo::randomPoints(predictors.dummy.current, n=10000, p=NULL, excludep=F)
extract.data <- extract(predictors.dummy.current, y=background)
colSums(extract.data)/sum(extract.data)*100
novel.levels
novel.levels.quantiles

par.old <- graphics::par(no.readonly=T)
graphics::par(mfrow=c(1,2))
raster::plot(novel.levels.raster, breaks=c(-0.1, 0, 1), col=c("grey", "green"), 
    main="novel outside minimum to maximum range")
rect(extent.values2[1], extent.values2[3], extent.values2[2], extent.values2[4])
raster::plot(novel.levels.quantiles.raster, breaks=c(-0.1, 0, 1), col=c("grey", "green"), 
    main="novel outside quantile range")
rect(extent.values2[1], extent.values2[3], extent.values2[2], extent.values2[4])
graphics::par(par.old)

}
}

