\name{ensemble.spatialThin}
\alias{ensemble.spatialThin}
\alias{ensemble.spatialThin.quant}
\alias{ensemble.environmentalThin}
\alias{ensemble.outlierThin}

\title{
Thinning of presence point coordinates in geographical or environmental space
}

\description{
Function \code{ensemble.spatialThin} creates a randomly selected subset of point coordinates where the shortest distance (geodesic) is above a predefined minimum. The geodesic is calculated more accurately (via \code{\link[geosphere]{distGeo}}) than in the \code{spThin} or \code{red} packages.
}

\usage{
ensemble.spatialThin(x, thin.km = 0.1, 
    runs = 100, verbose = FALSE, 
    return.notRetained = FALSE)

ensemble.spatialThin.quant(x, thin.km = 0.1, 
    runs = 100, verbose = FALSE, 
    LON.length = 21, LAT.length = 21)

ensemble.environmentalThin(x, predictors.stack = NULL, thin.n = 50,
    runs = 100, pca.var = 0.95, verbose = FALSE,
    return.notRetained = FALSE)

ensemble.outlierThin(x, predictors.stack = NULL, k = 10,
    quant = 0.95, pca.var = 0.95, 
    return.outliers = FALSE)
}

\arguments{
  \item{x}{ Point locations provided in 2-column (lon, lat) format.  }
  \item{thin.km}{ Threshold for minimum distance (km) in final point location data set.  }
  \item{runs}{ Number of runs to maximize the retained number of point coordinates.  }
  \item{verbose}{ Provide some details on each run.  }
  \item{return.notRetained}{ Return in an additional data set the point coordinates that were thinned out.  }

  \item{LON.length}{ Number of quantile limits to be calculated from longitudes; see also \code{\link[stats]{quantile}}  }
  \item{LAT.length}{ Number of quantile limits to be calculated from latitudes; see also \code{\link[stats]{quantile}}  }

  \item{predictors.stack}{ RasterStack object (\code{\link[raster]{stack}}) containing environmental layers that define the environmental space of point observations.  }
  \item{thin.n}{ Target number of environmentally thinned points.  }
  \item{pca.var}{ Minimum number of axes based on the fraction of variance explained (default value of 0.95 indicates that at least 95 percent of variance will be explained on the selected number of axes). Axes and coordinates are obtained from Principal Components Analysis (\code{\link[vegan]{scores}}).  }

  \item{k}{ The number of neighbours for the Local Outlier Factor analysis; see \code{\link[Rlof]{lof}}  }
  \item{quant}{ The quantile probability above with local outlier factors are classified as outliers; see also \code{\link[stats]{quantile}}  }
  \item{return.outliers}{ Return in an additional data set the point coordinates that were flagged as outliers.  }
}

\details{
Locations with distances smaller than the threshold distance are randomly removed from the data set until no distance is smaller than the threshold. The function uses a similar algorithm as functions in the \code{spThin} or \code{red} packages, but the geodesic is more accurately calculated via \code{\link[geosphere]{distGeo}}.

With several runs (default of 100 as in the \code{red} package or some \code{spThin} examples), the (first) data set with the maximum number of records is retained.

Function \code{ensemble.spatialThin.quant} was designed to be used with large data sets where the size of the object with pairwise geographical distances could create memory problems. With this function, spatial thinning is only done within geographical areas defined by quantile limits of geographical coordinates.

Function \code{ensemble.environmentalThin} performs an analysis in environmental space similar to the analysis in geographical space by \code{ensemble.spatialThin}. However, the target number of retained point coordinates needs to be defined by the user. Coordinates are obtained in environmental space by a principal components analysis (function \code{\link[vegan]{rda}}). Internally, first points are randomly selected from the pair with the smallest environmental distance until the selected target number of retained point coordinates is reached. From the retained point coordinates, the minimum environmental distance is determined. In a second step (more similar to spatial thinning), locations are randomly removed from all pairs that have a distance larger than the minimum distance calculated in step 1.

Function \code{ensemble.outlierThin} selects point coordinates that are less likely to be local outliers based on a Local Outlier Factor analysis (\code{\link[Rlof]{lof}}). Since LOF does not result in strict classification of outliers, a user-defined quantile probability is used to identify outliers.

}

\value{
The function returns a spatially or environmentally thinned point location data set. 

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
plot(predictors[[1]], main="5 runs", ext=extent(SpatialPoints(pres.thin1)))
points(pres, pch=20, col="black")
points(pres.thin1, pch=20, col="red")

pres.thin2 <- ensemble.spatialThin(pres, thin.km=100, runs=10, verbose=T)
plot(predictors[[1]], main="5 runs (after fresh start)", ext=extent(SpatialPoints(pres.thin2)))
points(pres, pch=20, col="black")
points(pres.thin2, pch=20, col="red")

pres.thin3 <- ensemble.spatialThin(pres, thin.km=100, runs=100, verbose=T)
plot(predictors[[1]], main="100 runs", ext=extent(SpatialPoints(pres.thin3)))
points(pres, pch=20, col="black")
points(pres.thin3, pch=20, col="red")

pres.thin4 <- ensemble.spatialThin(pres, thin.km=100, runs=100, verbose=T)
plot(predictors[[1]], main="100 runs (after fresh start)", ext=extent(SpatialPoints(pres.thin4)))
points(pres, pch=20, col="black")
points(pres.thin4, pch=20, col="red")

graphics::par(par.old)

## thinning in environmental space

env.thin <- ensemble.environmentalThin(pres, predictors.stack=predictors, thin.n=60,
    return.notRetained=T)
pres.env1 <- env.thin$retained
pres.env2 <- env.thin$not.retained

# plot in geographical space
par.old <- graphics::par(no.readonly=T)
par(mfrow=c(1, 2))

plot(predictors[[1]], main="black = not retained", ext=extent(SpatialPoints(pres.thin3)))
points(pres.env2, pch=20, col="black")
points(pres.env1, pch=20, col="red")

# plot in environmental space
background.data <- data.frame(raster::extract(predictors, pres))
rda.result <- vegan::rda(X=background.data, scale=T)
# select number of axes
ax <- 2
while ( (sum(vegan::eigenvals(rda.result)[c(1:ax)])/
    sum(vegan::eigenvals(rda.result))) < 0.95 ) {ax <- ax+1}
rda.scores <- data.frame(vegan::scores(rda.result, display="sites", scaling=1, choices=c(1:ax)))
rownames(rda.scores) <- rownames(pres)
points.in <- rda.scores[which(rownames(rda.scores) \%in\% rownames(pres.env1)), c(1:2)]
points.out <- rda.scores[which(rownames(rda.scores) \%in\% rownames(pres.env2)), c(1:2)]
plot(points.out, main="black = not retained", pch=20, col="black", 
    xlim=range(rda.scores[, 1]), ylim=range(rda.scores[, 2]))
points(points.in, pch=20, col="red")

graphics::par(par.old)

## removing outliers
out.thin <- ensemble.outlierThin(pres, predictors.stack=predictors, k=10,
    return.outliers=T)
pres.out1 <- out.thin$inliers
pres.out2 <- out.thin$outliers

# plot in geographical space
par.old <- graphics::par(no.readonly=T)
par(mfrow=c(1, 2))

plot(predictors[[1]], main="black = outliers", ext=extent(SpatialPoints(pres.thin3)))
points(pres.out2, pch=20, col="black")
points(pres.out1, pch=20, col="red")

# plot in environmental space
background.data <- data.frame(raster::extract(predictors, pres))
rda.result <- vegan::rda(X=background.data, scale=T)
# select number of axes
ax <- 2
while ( (sum(vegan::eigenvals(rda.result)[c(1:ax)])/
    sum(vegan::eigenvals(rda.result))) < 0.95 ) {ax <- ax+1}
rda.scores <- data.frame(vegan::scores(rda.result, display="sites", scaling=1, choices=c(1:ax)))
rownames(rda.scores) <- rownames(pres)
points.in <- rda.scores[which(rownames(rda.scores) \%in\% rownames(pres.out1)), c(1:2)]
points.out <- rda.scores[which(rownames(rda.scores) \%in\% rownames(pres.out2)), c(1:2)]
plot(points.out, main="black = outliers", pch=20, col="black", 
    xlim=range(rda.scores[, 1]), ylim=range(rda.scores[, 2]))
points(points.in, pch=20, col="red")

graphics::par(par.old)

}

}