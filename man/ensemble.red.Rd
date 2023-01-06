\name{ensemble.red}
\alias{ensemble.red}
\alias{ensemble.chull.create}
\alias{ensemble.chull.apply}
\alias{ensemble.chull.buffer.distances}
\alias{ensemble.chull.MSDM}

\title{
Area of Occupancy (AOO) and Extent of Occurrence (EOO) via the \pkg{red} library.
}

\description{
Function \code{ensemble.red} is a wrapper function for estimation of AOO and EOO computed for redlisting of species based on IUCN criteria (\url{https://www.iucnredlist.org/about/regional}). Function \code{ensemble.chull.create} creates a mask layer based on a convex hull around known presence locations, inspired by \code{mcp} argument of the \code{\link[red]{map.sdm}} function. 
}

\usage{

ensemble.red(x)

ensemble.chull.create(x.pres = NULL, p = NULL, buffer.width = 0.2,
    buffer.maxmins = FALSE, lonlat.dist = FALSE,
    poly.only = FALSE,
    RASTER.format = "GTiff", RASTER.datatype = "INT1U", RASTER.NAflag = 255,
    overwrite = TRUE, ...)

ensemble.chull.apply(x.spec = NULL, mask.layer=NULL, keep.old=T,
    RASTER.format="GTiff", RASTER.datatype="INT1U", RASTER.NAflag = 255,
    overwrite=TRUE, ...)

ensemble.chull.buffer.distances(p = NULL, 
    buffer.maxmins = FALSE, lonlat.dist = FALSE)

ensemble.chull.MSDM(p = NULL, a = NULL, species.name = NULL,
    suit.file = NULL, suit.divide = 1000, MSDM.dir = NULL,
    method = "BMCP", threshold = "spec_sens", 
    buffer = "species_specific")

}

\arguments{
  \item{x}{ RasterLayer object (\code{\link[raster]{raster}}), representing 'count' suitability layers (available from the 'count' and 'consensuscount' subdirectories of the 'ensembles' directory)  }

  \item{x.pres}{ RasterLayer object (\code{\link[raster]{raster}}), representing 'presence' suitability layers (available from the 'presence' and 'consensuspresence' subdirectories of the 'ensembles' directory)  }
  \item{p}{ known presence locations, available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{buffer.width}{ multiplier to create buffer (via \code{\link[rgeos]{gBuffer}}) by multiplying the maximum distance among the presence locations (calculated via \code{\link[raster]{pointDistance}})  }
  \item{buffer.maxmins}{ Calculate the buffer width based on the two neighbouring locations that are furthest apart (maximum of minimum distances from each location).  }
  \item{lonlat.dist}{ Estimate the distance in km for longitude latitude data.  }
  \item{poly.only}{ Only return the polygon with the convex hull, but do not create the mask layer.  }
  \item{RASTER.format}{ Format of the raster files that will be generated. See \code{\link[raster]{writeFormats}} and \code{\link[raster]{writeRaster}}.  }
  \item{RASTER.datatype}{ Format of the raster files that will be generated. See \code{\link[raster]{dataType}} and \code{\link[raster]{writeRaster}}.  }
  \item{RASTER.NAflag}{ Value that is used to store missing data. See \code{\link[raster]{writeRaster}}.  }
  \item{overwrite}{ Overwrite existing raster files. See \code{\link[raster]{writeRaster}}.  }
  \item{...}{ Additional arguments for \code{\link[raster]{writeRaster}}.  }

  \item{x.spec}{ RasterLayer object (\code{\link[raster]{raster}}), representing any suitability layer for the species under investigation)  }
  \item{mask.layer}{ RasterLayer object (\code{\link[raster]{raster}}), representing the mask based on the convex hull around known presence locations. The function will replace all values in \code{x.spec} to zero where corresponding values in the mask.layer are zero.  }
  \item{keep.old}{ keep a copy of the RasterLayer before the mask is applied.  }

  \item{a}{ absence of background locations, available in 2-column (lon, lat) dataframe.  }
  \item{species.name}{ name of the species, ideally without spaces.  }
  \item{suit.file}{ file with raster data corresponding to suitability values of the focal species.  }
  \item{suit.divide}{ number by which values in the suitability raster should be divided to result in probabilities (BiodiversityR saves data as 1000 * suitability, hence these values need to be divided by 1000).  }
  \item{MSDM.dir}{ name of the directory where input and processed raster files will be saved.  }
  \item{method}{ method for MSDM_Posteriori function from c("OBR", "PRES", "LQ", "MCP", "BMCP").  }
  \item{threshold}{ threshold for MSDM_Posteriori function from c("kappa", "spec_sens", "no_omission", "prevalence", "equal_sens_spec",
    "sensitivty").  }
  \item{buffer}{ buffer for MSDM_Posteriori function.  }


}

\details{
Function \code{ensemble.red} calculates AOO (\code{\link[red]{aoo}}) and EOO (\code{\link[red]{aoo}}) statistics calculated for areas with different consensus levels on species presence (1 model predicting presence, 2 models predicting presence, ...). In case that these statistics are within IUCN criteria for Critically Endangered (CR), Endangered (EN) or Vulnerable (VU), then this information is added in columns documenting the types of AOO and EOO.

Function \code{ensemble.chull.create} first creates a convex hull around known presence locations. Next, a buffer is created around the convex hull where the width of this buffer is calculated as the maximum distance among presence locations (\code{\link[raster]{pointDistance}}) multiplied by argument \code{buffer.width}. Finally, the mask is created by including all polygons of predicted species presence that are partially covered by the convex hull and its buffer.
}

\value{
Function \code{ensemble.red} returns an array with AOO and EOO 
Function \code{ensemble.chull.create} creates a mask layer based on a convex hull around known presence locations.
Function \code{ensemble.chull.MSDM} prepares the input data and script for the MSDM_Posteriori function of the MSDM package.
}

\references{
Cardoso P. 2017. red - an R package to facilitate species red list assessments according to the IUCN criteria. Biodiversity Data Journal 5:e20530. \doi{10.3897/BDJ.5.e20530}

Mendes, P.; Velazco S.J.E.; Andrade, A.F.A.; De Marco, P. (2020) Dealing with overprediction in species distribution models: how adding distance constraints can improve model accuracy, Ecological Modelling, in press. \doi{10.1016/j.ecolmodel.2020.109180}

Kindt R. 2018. Ensemble species distribution modelling with transformed suitability values. Environmental Modelling & Software 100: 136-145. \doi{10.1016/j.envsoft.2017.11.009}
}

\author{ Roeland Kindt (World Agroforestry Centre)}

\seealso{
\code{\link{ensemble.batch}}
}

\examples{

\dontrun{

## Not run: 
# based on examples in the dismo package

# get predictor variables
library(dismo)
predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
    pattern='grd', full.names=TRUE)
predictors <- stack(predictor.files)
# subset based on Variance Inflation Factors
predictors <- subset(predictors, subset=c("bio5", "bio6", 
    "bio16", "bio17"))
predictors
predictors@title <- "red"

# presence points
presence_file <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
pres <- read.table(presence_file, header=TRUE, sep=',')

# fit 4 ensemble models (could take some time!)
# (examples for the red package use 100 models)
ensembles <- ensemble.batch(x=predictors, 
    xn=c(predictors),
    species.presence=pres, 
    thin.km=100,
    k.splits=4, k.test=0, 
    n.ensembles=4, 
    SINK=TRUE, 
    ENSEMBLE.best=10, ENSEMBLE.exponent=c(1, 2, 3), 
    ENSEMBLE.min=0.6,
    MAXENT=0, MAXNET=1, MAXLIKE=1, GBM=1, GBMSTEP=0, RF=1, CF=0,
    GLM=0, GLMSTEP=1, GAM=1, GAMSTEP=0, MGCV=1, MGCVFIX=0, 
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, 
    BIOCLIM.O=0, BIOCLIM=1, DOMAIN=0, MAHAL=0, MAHAL01=0,
    PROBIT=TRUE,
    Yweights="BIOMOD",
    formulae.defaults=TRUE)

# first application of ensemble.red before applying the convex hull mask
# AOO and EOO are determined for each count level
library(red)
count.file <- paste(getwd(), 
    "/ensembles/consensuscount/Bradypus variegatus_red.tif", sep="")
count.raster <- raster(count.file)
ensemble.red(count.raster)

# do not predict presence in polygons completely outside convex hull
# of known presence locations
pres.file <- paste(getwd(), 
    "/ensembles/consensuspresence/Bradypus variegatus_red.tif", sep="")
pres.raster <- raster(pres.file)
pres1 <- pres[, -1]
chull.created <- ensemble.chull.create(x.pres=pres.raster, p=pres1)

mask.raster <- chull.created$mask.layer
plot(mask.raster, col=c("black", "green"))
mask.poly <- chull.created$convex.hull

pres.chull <- ensemble.chull.apply(pres.raster, mask=mask.raster, keep.old=T)

par.old <- graphics::par(no.readonly=T)
par(mfrow=c(1,2))
plot(pres.raster, breaks=c(-1, 0, 1), col=c("grey", "green"),
    main="before convex hull")
points(pres1, col="blue")

# load new
pres.file.new <- paste(getwd(), 
    "/ensembles/chull/Bradypus variegatus_red.tif", sep="")
pres.raster.new <- raster(pres.file.new)
plot(pres.raster.new, breaks=c(-1, 0, 1), col=c("grey", "green"),
    main="after convex hull")
plot(mask.poly, add=T, border="blue")

# create a smaller hull (0.05 * largest distance)
chull.created <- ensemble.chull.create(x.pres=pres.raster, p=pres1, 
    buffer.width=0.05, lonlat.dist=TRUE)
mask.raster <- chull.created$mask.layer
mask.poly <- chull.created$convex.hull
pres.chull <- ensemble.chull.apply(pres.raster, mask=mask.raster, keep.old=T)

par(mfrow=c(1,2))
plot(pres.raster, breaks=c(-1, 0, 1), col=c("grey", "green"),
    main="before convex hull")
points(pres1, col="blue")
pres.raster.new <- raster(pres.file.new)
plot(pres.raster.new, breaks=c(-1, 0, 1), col=c("grey", "green"),
    main="after convex hull")
plot(mask.poly, add=T, border="blue")

# create a hull based on the distance to the location with the farthest neighbour
chull.created <- ensemble.chull.create(x.pres=pres.raster, p=pres1, 
    buffer.maxmins=TRUE, buffer.width=0.9, lonlat.dist=TRUE)
mask.raster <- chull.created$mask.layer
mask.poly <- chull.created$convex.hull
pres.chull <- ensemble.chull.apply(pres.raster, mask=mask.raster, keep.old=T)

par(mfrow=c(1,2))
plot(pres.raster, breaks=c(-1, 0, 1), col=c("grey", "green"),
    main="before convex hull")
points(pres1, col="blue")
pres.raster.new <- raster(pres.file.new)
plot(pres.raster.new, breaks=c(-1, 0, 1), col=c("grey", "green"),
    main="after convex hull")
plot(mask.poly, add=T, border="blue")

par.old <- graphics::par(no.readonly=T)

# how distances were derived
# maximum distance between observations
ensemble.chull.buffer.distances(pres1, lonlat.dist=TRUE)   
# the closest neigbhour that is farthest away from each observation
# this is the distance calculated by MSDM_posteriori for buffer="species_specific"
ensemble.chull.buffer.distances(pres1, buffer.maxmins=TRUE, lonlat.dist=TRUE) 

}
}

