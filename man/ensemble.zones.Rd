\name{ensemble.zones}
\alias{ensemble.zones}
\alias{ensemble.centroids}

\title{
Mapping of environmental zones based on the Mahalanobis distance from centroids in environmental space.
}

\description{
Function \code{ensemble.zones} maps the zone of each raster cell within a presence map based on the minimum Mahalanobis distance (via \code{\link[stats]{mahalanobis}}) to different centroids. Function \code{ensemble.centroids} defines centroids within a presence map based on Principal Components Analysis (via \code{\link[vegan]{rda}}) and K-means clustering (via \code{\link[stats]{kmeans}}). 
}

\usage{
ensemble.zones(presence.raster = NULL, centroid.object = NULL, 
    x = NULL, ext = NULL,
    RASTER.species.name = centroid.object$name, RASTER.stack.name = x@title, 
    RASTER.format = "raster", RASTER.datatype = "INT1S", RASTER.NAflag = -127, 
    KML.out = FALSE, KML.maxpixels = 100000, KML.blur = 10)

ensemble.centroids(presence.raster = NULL, x = NULL, categories.raster = NULL,
    an = 10000, ext = NULL, name = "Species001", 
    pca.var = 0.95, centers = 0, use.silhouette = TRUE, 
    plotit = FALSE, dev.new.width = 7, dev.new.height = 7)

}

\arguments{
  \item{presence.raster}{ RasterLayer object (\code{\link[raster]{raster}}) documenting presence (coded 1) of an organism }
  \item{centroid.object}{ Object listing values for centroids and covariance to be used with the \code{\link[stats]{mahalanobis}} distance (used internally by the prediction function called from \code{\link[raster]{predict}}). }
  \item{x}{ RasterStack object (\code{\link[raster]{stack}}) containing all environmental layers that correspond to explanatory variables }
  \item{ext}{ an Extent object to limit the predictions and selection of background points to a sub-region of \code{presence.raster} and \code{x}, typically provided as c(lonmin, lonmax, latmin, latmax). See also \code{\link[dismo]{randomPoints}} and \code{\link[raster]{extent}}. }
  \item{RASTER.species.name}{ First part of the names of the raster file that will be generated, expected to identify the modelled species (or organism) }
  \item{RASTER.stack.name}{ Last part of the names of the raster file that will be generated, expected to identify the predictor stack used }
  \item{RASTER.format}{ Format of the raster files that will be generated. See \code{\link[raster]{writeFormats}} and \code{\link[raster]{writeRaster}}.  }
  \item{RASTER.datatype}{ Format of the raster files that will be generated. See \code{\link[raster]{dataType}} and \code{\link[raster]{writeRaster}}.  }
  \item{RASTER.NAflag}{ Value that is used to store missing data. See \code{\link[raster]{writeRaster}}.  }
  \item{KML.out}{ If \code{TRUE}, then kml files will be saved in a subfolder 'kml/zones'.  }
  \item{KML.maxpixels}{ Maximum number of pixels for the PNG image that will be displayed in Google Earth. See also \code{\link[raster]{KML}}.  }
  \item{KML.blur}{ Integer that results in increasing the size of the PNG image by \code{KML.blur^2}, which may help avoid blurring of isolated pixels. See also \code{\link[raster]{KML}}.  }

  \item{categories.raster}{ RasterLayer object (\code{\link[raster]{raster}}) documenting predefined zones such as vegetation types. In case this object is provided, then centroids will be calculated for each zone. }
  \item{an}{ Number of presence points to be used for Principal Components Analysis (via \code{\link[vegan]{rda}}); see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{name}{ Name for the centroid object, for example identifying the species and area for which centroids are calculated }
  \item{pca.var}{ Minimum number of axes based on the fraction of variance explained (default value of 0.95 indicates that at least 95 percent of variance will be explained on the selected number of axes). Axes and coordinates are obtained from Principal Components Analysis (\code{\link[vegan]{scores}}).  }
  \item{centers}{ Number of centers (clusters) to be used for K-means clustering (\code{\link[stats]{kmeans}}). In case a value smaller than 1 is provided, function \code{\link[vegan]{cascadeKM}} is called to determine the optimal number of centers via the Calinski-Harabasz criterion.  }
  \item{use.silhouette}{ If \code{TRUE}, then centroid values are only based on presence points that have silhouette values (\code{\link[cluster]{silhouette}}) larger than 0.  }
  \item{plotit}{ If \code{TRUE}, then a plot is provided that shows the locations of centroids in geographical and environmental space. Plotting in geographical space is based on determination of the presence location (analogue) with smallest Mahalanobis distance to the centroid in environmental space.  }
  \item{dev.new.width}{ Width for new graphics device (\code{\link[grDevices]{dev.new}}). If < 0, then no new graphics device is opened.  }
  \item{dev.new.height}{ Heigth for new graphics device (\code{\link[grDevices]{dev.new}}). If < 0, then no new graphics device is opened.  }
}

\details{
Function \code{ensemble.zones} maps the zone of each raster cell of a predefined presence map, whereby the zone is defined as the centroid with the smallest Mahalanobis distance. The function returns a RasterLayer object (\code{\link[raster]{raster}}) and possibly a KML layer.

Function \code{ensemble.centroid} provides the centroid locations in environmental space and a covariance matrix (\code{\link[stats]{cov}}) to be used with \code{\link[stats]{mahalanobis}}. Also provided is information on the analogue presence location that is closest to the centroid in environmental space. 

}

\value{
Function \code{ensemble.centroid} returns a list with following objects:
\item{centroids }{ Location of centroids in environmental space }
\item{centroid.analogs }{ Location of best analogs to centroids in environmental space }
\item{cov.mahal }{ Covariance matrix }
}

\author{ Roeland Kindt (World Agroforestry Centre)}

\seealso{
\code{\link{ensemble.raster}}
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

# choose background points
background <- randomPoints(predictors, n=1000, extf=1.00)

# predicted presence from GLM
ensemble.calibrate.step1 <- ensemble.calibrate.models(
    x=predictors, p=pres, a=background,
    species.name="Bradypus",
    MAXENT=0, MAXLIKE=0, GBM=0, GBMSTEP=0, RF=0, GLM=1, GLMSTEP=0, 
    GAM=0, GAMSTEP=0, MGCV=0, MGCVFIX=0,
    EARTH=0, RPART=0, NNET=0, FDA=0, SVM=0, SVME=0, GLMNET=0,
    BIOCLIM.O=0, BIOCLIM=0, DOMAIN=0, MAHAL=0, MAHAL01=0,
    Yweights="BIOMOD",
    models.keep=TRUE)

ensemble.raster.results <- ensemble.raster(xn=predictors, 
    models.list=ensemble.calibrate.step1$models, 
    RASTER.species.name="Bradypus", RASTER.stack.name="base")

# get presence map as for example created with ensemble.raster in subfolder 'ensemble/presence'
# presence values are values equal to 1
presence.file <- paste("ensembles//presence//Bradypus_base.grd", sep="")
presence.raster <- raster(presence.file)

# let cascadeKM decide on the number of clusters
dev.new()
centroids <- ensemble.centroids(presence.raster=presence.raster, 
    x=predictors, an=1000, plotit=T)
ensemble.zones(presence.raster=presence.raster, centroid.object=centroids, 
    x=predictors, RASTER.species.name="Bradypus", KML.out=T)

dev.new()
zones.file <- paste("ensembles//zones//Bradypus_base.grd", sep="")
zones.raster <- raster(zones.file)
max.zones <- maxValue(zones.raster)
plot(zones.raster, breaks=c(0, c(1:max.zones)), 
    col = grDevices::rainbow(n=max.zones), main="zones")
ensemble.zones(presence.raster=presence.raster, centroid.object=centroids, 
    x=predictors, RASTER.species.name="Bradypus", KML.out=T)

# manually choose 6 zones
dev.new()
centroids6 <- ensemble.centroids(presence.raster=presence.raster, 
    x=predictors, an=1000, plotit=T, centers=6)
ensemble.zones(presence.raster=presence.raster, centroid.object=centroids6, 
    x=predictors, RASTER.species.name="Bradypus6", KML.out=T)

dev.new()
zones.file <- paste("ensembles//zones//Bradypus6_base.grd", sep="")
zones.raster <- raster(zones.file)
max.zones <- maxValue(zones.raster)
plot(zones.raster, breaks=c(0, c(1:max.zones)), 
    col = grDevices::rainbow(n=max.zones), main="six zones")

}
}

