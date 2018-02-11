\name{ensemble.spatialThin}
\alias{ensemble.spatialThin}

\title{
Spatial thinning of presence point locations using the highly accurate geodesic estimates from the geosphere package
}

\description{
The function creates a randomly selected subset of point locations where the shortest distance (geodesic) is above a predefined minimum. The geodesic is calculated more accurately (via \code{\link[geosphere]{distGeo}}) than in the \code{spThin} or \code{red} packages.
}

\usage{
ensemble.spatialThin(x, thin.km=0.1, 
    runs=100, verbose=FALSE, return.matrix=FALSE)
}

\arguments{
  \item{x}{ Point locations provided in 2-column (lon, lat) format.  }
  \item{thin.km}{ Threshold for minimum distance (km) in final point location data set.  }
  \item{runs}{ Number of runs to maximize the retained number of point locations.  }
  \item{verbose}{ Provide some details on each run.  }
  \item{return.matrix}{ Include the distance matrix (calculated via \code{\link[geosphere]{distGeo}})) in the results.  }
}

\details{
Locations with distances smaller than the threshold distance are randomly removed from the data set until no distance is smaller than the threshold. The function uses a similar algorithm as functions in the \code{spThin} or \code{red} packages, but the geodesic is more accurately calculated via \code{\link[geosphere]{distGeo}}.

With several runs (default of 100 as in the \code{red} package or some \code{spThin} examples), the (first) data set with the maximum number of records is retained.
}

\value{
The function returns a spatially thinned point location data set. 

}

\references{
Aiello-Lammens ME, Boria RA, Radosavljevic A, Vilela B and Anderson RP. 2015. spThin: an R package for spatial thinning of species occurrence records for use in ecological niche models. Ecography 38: 541-545
}


\author{ Roeland Kindt (World Agroforestry Centre) }


\seealso{
\code{\link{ensemble.batch}}
}

\examples{
\dontrun{
# get predictor variables, only needed for plotting
library(dismo)
predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
    pattern='grd', full.names=TRUE)
predictors <- stack(predictor.files)
# subset based on Variance Inflation Factors
predictors <- subset(predictors, subset=c("bio5", "bio6", 
    "bio16", "bio17", "biome"))
predictors
predictors@title <- "base"

# presence points
presence_file <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
pres <- read.table(presence_file, header=TRUE, sep=',')[, -1]

# number of locations
nrow(pres)

par.old <- graphics::par(no.readonly=T)
par(mfrow=c(2,2))

pres.thin1 <- ensemble.spatialThin(pres, thin.km=100, runs=10, verbose=T)
plot(predictors[[1]], main="10 runs", ext=extent(SpatialPoints(pres.thin1)))
points(pres.thin1, pch=20, col="red")

pres.thin2 <- ensemble.spatialThin(pres, thin.km=100, runs=10, verbose=T)
plot(predictors[[1]], main="10 runs", ext=extent(SpatialPoints(pres.thin2)))
points(pres.thin2, pch=20, col="red")

pres.thin3 <- ensemble.spatialThin(pres, thin.km=100, runs=100, verbose=T)
plot(predictors[[1]], main="100 runs", ext=extent(SpatialPoints(pres.thin3)))
points(pres.thin3, pch=20, col="red")

pres.thin4 <- ensemble.spatialThin(pres, thin.km=100, runs=100, verbose=T)
plot(predictors[[1]], main="100 runs", ext=extent(SpatialPoints(pres.thin4)))
points(pres.thin4, pch=20, col="red")

graphics::par(par.old)

}

}