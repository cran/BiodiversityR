\name{ensemble.ecocrop}
\alias{ensemble.ecocrop}
\alias{ensemble.ecocrop.object}

\title{
Mapping of novel environmental conditions (areas where some of the environmental conditions are outside the range of environmental conditions of a reference area).
}

\description{
Function \code{ensemble.novel} creates the map with novel conditions. Function \code{ensemble.novel.object} provides the reference values used by the prediction function used by \code{\link[raster]{predict}} . 
}

\usage{
ensemble.ecocrop(x = NULL, ecocrop.object = NULL, 
    RASTER.object.name = ecocrop.object$name, RASTER.stack.name = x@title,
    RASTER.format = "raster", RASTER.datatype = "INT2S", RASTER.NAflag = -32767,
    KML.out = TRUE, KML.blur = 10, KML.maxpixels = 100000,
    CATCH.OFF = FALSE)

ensemble.ecocrop.object(temp.thresholds, rain.thresholds, name = "crop01", 
    temp.multiply = 10, annual.temps = TRUE, transform = 1)

}

\arguments{
  \item{x}{ RasterStack object (\code{\link[raster]{stack}}) containing all environmental layers for which suitability should be calculated.  }
  \item{ecocrop.object}{ Object listing optimal and absolute minima and maxima for the rainfall and temperature values, used by the prediction function that is used internally by \code{\link[raster]{predict}}. This object is created with \code{\link{ensemble.ecocrop.object}}. }
  \item{RASTER.object.name}{ First part of the names of the raster file that will be generated, expected to identify the species or crop for which ranges were calculated }
  \item{RASTER.stack.name}{ Last part of the names of the raster file that will be generated, expected to identify the predictor stack used }
  \item{RASTER.format}{ Format of the raster files that will be generated. See \code{\link[raster]{writeFormats}} and \code{\link[raster]{writeRaster}}.  }
  \item{RASTER.datatype}{ Format of the raster files that will be generated. See \code{\link[raster]{dataType}} and \code{\link[raster]{writeRaster}}.  }
  \item{RASTER.NAflag}{ Value that is used to store missing data. See \code{\link[raster]{writeRaster}}.  }
  \item{KML.out}{ If \code{TRUE}, then kml files will be saved in a subfolder 'kml/zones'.  }
  \item{KML.maxpixels}{ Maximum number of pixels for the PNG image that will be displayed in Google Earth. See also \code{\link[raster]{KML}}.  }
  \item{KML.blur}{ Integer that results in increasing the size of the PNG image by \code{KML.blur^2}, which may help avoid blurring of isolated pixels. See also \code{\link[raster]{KML}}.  }
  \item{CATCH.OFF}{ Disable calls to function \code{\link[base]{tryCatch}}.  }

  \item{temp.thresholds}{ Optimal and absolute thresholds for temperatures. These will be sorted as: absolute minimum temperature, optimal minimum temperature, optimal maximum temperature and absolute maximum temperature.  }
  \item{rain.thresholds}{ Optimal and absolute thresholds for annual rainfall. These will be sorted as: absolute minimum rainfall, optimal minimum rainfall, optimal maximum rainfall and absolute maximum rainfall.  }
  \item{name}{ Name of the object, expect to expected to identify the species or crop }
  \item{temp.multiply}{ Multiplier for temperature values. Default of 10 is to be used with raster layers where temperature was multiplied by 10 such as Worldclim or AFRICLIM.   }
  \item{annual.temps}{ If \code{TRUE} then temperature limits are assumed to apply to mean annual temperature (bioclimatic variable bio1). If \code{FALSE} then minimum temperature limits are assumed to apply to the temperature of the coldest month (bioclimatic variable bio6) and maximum temperature limits are assumed to apply to the temperature of the hottest month (bioclimatic variable bio5). See also \code{\link[dismo]{biovars}}.   } 
  \item{transform}{ Exponent used to transform probability values obtained from interpolating between optimal and absolute limits. Exponent of 2 results in squaring probabilities, for example input probabilities of 0.5 transformed to 0.5^2 = 0.25.  }
}

\details{
Function \code{ensemble.ecocrop} maps suitability for a species or crop based on optimal and absolute temperature and rainfall limits. Where both temperature and rainfall are within the optimal limits, suitability of 1000 is calculated. Where both temperature and rainfall are outside the absolute limits, suitability of 0 is calculated. In situations where temperature or rainfall is in between the optimal and absolute limits, then suitability is interpolated between 0 and 1000, and the lowest suitability from temperature and rainfall is calculated. Setting very wide rainfall limits will simulate the effect of irrigation, i.e. where suitability only depends on temperature limits. 

For a large range of crop and plant species, optimal and absolute limits are available from the FAO ecocrop database (\url{http://ecocrop.fao.org/ecocrop}), hence the name of the function. A different implementation of suitability mapping based on ecocrop limits is available from \code{\link[dismo]{ecocrop}}. Ecocrop thresholds for several species are available from: \code{\link[dismo]{getCrop}}

}

\value{
Function \code{ensemble.ecocrop.object} returns a list with following objects:
\item{name }{ name for the crop or species }
\item{temp.thresholds }{ optimal and absolute minimum and maximum temperature limits }
\item{rain.thresholds }{ optimal and absolute minimum and maximum annual rainfall limits }
\item{annual.temps }{ logical indicating whether temperature limits apply to annual temperatures }
\item{transform }{ exponent to transform suitability values }
}

\author{ Roeland Kindt (World Agroforestry Centre)}

\seealso{
\code{\link[dismo]{biovars}}
}

\examples{

\dontrun{
#test with Brazil nut (limits from FAO ecocrop)
#temperature: (12) 20-36 (40)
#annnual rainfall: (1400) 2400-2800 (3500)

# get predictor variables
library(dismo)
predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
    pattern='grd', full.names=TRUE)
predictors <- stack(predictor.files)
# subset based on Variance Inflation Factors
predictors <- subset(predictors, subset=c("bio5", "bio6", "bio12"))
predictors
predictors@title <- "base"

Brazil.ecocrop <- ensemble.ecocrop.object(temp.thresholds=c(20, 36, 12, 40), 
    rain.thresholds=c(2400, 2800, 1400, 3500), 
    annual.temps=FALSE, name="Bertholletia_excelsa")
Brazil.ecocrop
ensemble.ecocrop(predictors, ecocrop.object=Brazil.ecocrop)

dev.new()
par.old <- graphics::par(no.readonly=T)
graphics::par(mfrow=c(1,2))


rasterfull1 <- paste("ensembles//ecocrop//Bertholletia_excelsa_base.grd", sep="")
rasterfull1 <- raster(rasterfull1)
# raster file saved probabilities as integer values between 0 and 1000
rasterfull1 <- rasterfull1/1000
raster::plot(rasterfull1, main="Ecocrop suitability")

GBIFloc <- gbif(genus="Bertholletia", species="excelsa", geo=TRUE)
GBIFpres <- cbind(GBIFloc$lon, GBIFloc$lat)
GBIFpres <- GBIFpres[complete.cases(GBIFpres), ]
GBIFpres <- GBIFpres[duplicated(GBIFpres) == FALSE, ]
point.suitability <- extract(rasterfull1, y=GBIFpres)
point.suitability[is.na(point.suitability)] <- -1

GBIFpres.optimal <- GBIFpres[point.suitability == 1, ]
GBIFpres.suboptimal <- GBIFpres[point.suitability < 1 & point.suitability > 0, ]
GBIFpres.not <- GBIFpres[point.suitability == 0, ]

raster::plot(rasterfull1, main="GBIF locations", 
    sub="blue: optimal, cyan: suboptimal, red: not suitable")
bg.legend <- c("blue", "cyan", "red")

points(GBIFpres.suboptimal, pch=21, cex=1.2, bg=bg.legend[2])
points(GBIFpres.optimal, pch=21, cex=1.2, bg=bg.legend[1])
points(GBIFpres.not, pch=21, cex=1.2, bg=bg.legend[3])

graphics::par(par.old)
}

}

