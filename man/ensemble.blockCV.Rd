\name{ensemble.spatialBlock}
\alias{ensemble.spatialBlock}
\alias{ensemble.envBlock}

\title{
Spatially or environmentally separated folds for cross-validation via blockCV::spatialBlock or blockCV::envBlock
}

\description{
The functions internally calls blockCV::spatialBlock and blockCV::envBlock. Syntax is very similar to that of BiodiversityR::ensemble.calibrate.weights.
}

\usage{
ensemble.spatialBlock(x = NULL, p = NULL, 
    a = NULL, an = 1000, EPSG=NULL,
    excludep = FALSE, target.groups = FALSE, k = 4,
    factors = NULL,
    theRange = NULL, return.object = FALSE, ...)

ensemble.envBlock(x = NULL, p = NULL, 
    a = NULL, an = 1000, EPSG=NULL,
    excludep = FALSE, target.groups = FALSE, k = 4,
    factors = NULL,
    return.object = FALSE, ...)
}

\arguments{
  \item{x}{ RasterStack object (\code{\link[raster]{stack}}) containing all layers that correspond to explanatory variables  }
  \item{p}{ presence points used for calibrating the suitability models, typically available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{a}{ background points used for calibrating the suitability models, typically available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{an}{ number of background points for calibration to be selected with \code{\link[dismo]{randomPoints}} in case argument \code{a} is missing}
  \item{EPSG}{ EPSG number (\url{https://spatialreference.org/}) to be assigned internally to the coordinate reference system of the locations via \code{\link[sf]{st_crs}}. Although the function internally first assigns the coordinate reference from the RasterStack \code{x} via \code{\link[raster]{crs}}, this method fails in some situations as in the example shown below. In such cases, manually assigning the EPSG could resolve this problem. }
  \item{excludep}{ parameter that indicates (if \code{TRUE}) that presence points will be excluded from the background points; see also \code{\link[dismo]{randomPoints}} }
  \item{target.groups}{ Parameter that indicates (if \code{TRUE}) that the provided background points (argument \code{a}) represent presence points from a target group sensu Phillips et al. 2009 (these are species that are all collected or observed using the same methods or equipment). Setting the parameter to \code{TRUE} results in selecting the centres of cells of the target groups as background points, while avoiding to select the same cells twice. Via argument \code{excludep}, it is possible to filter out cells with presence observations (argument \code{p}).  }
  \item{k}{ Integer value. The number of desired folds for cross-validation. The default is \code{k = 4}. The interpretation of the argument is exactly the same as in \code{\link[BiodiversityR]{ensemble.calibrate.models}} and \code{\link[dismo]{kfold}}.  }
  \item{factors}{ vector that indicates which variables are factors; see also \code{\link[dismo]{prepareData}}  }
  \item{theRange}{ Numeric value of the specified range by which blocks are created and training/testing data are separated. This distance should be in metres. See also \code{\link[blockCV]{spatialBlock}}.  }
  \item{return.object}{ If \code{TRUE}, then also return ('block.object') the complete result of \code{\link[blockCV]{spatialBlock}} or \code{\link[blockCV]{envBlock}}. In addtion (if \code{TRUE}), return the species data ('speciesData') that was created for \code{blockCV}. To visualize these results, see below or from \code{\link[blockCV]{foldExplorer}}.  }
  \item{...}{ Other arguments to pass to \code{\link[blockCV]{spatialBlock}} or \code{\link[blockCV]{envBlock}}, such \code{numLimit} (The minimum number of points in each fold for training-presence, training-absence, testing-presence and testing-absence) and \code{iteration} (The number of attempts to create folds that fulfil the \code{numLimit} requirement).   }

}

\details{
The functions internally call \code{\link[blockCV]{spatialBlock}} or \code{\link[blockCV]{envBlock}}. 

The result of the function includes a list (\code{k}) with following elements. This list can be directly imported into \code{\link[BiodiversityR]{ensemble.calibrate.weights}}, but only elements \code{groupp} and \code{groupa} will be used.

     - \code{p} : Presence locations, created by \code{\link[BiodiversityR]{ensemble.calibrate.models}} where points with missing data were excluded and possibly points were added for missing factor levels 

     - \code{a} : Background locations, created by \code{\link[BiodiversityR]{ensemble.calibrate.models}} where points with missing data were excluded and possibly points were added for missing factor levels 

     - \code{groupp} : k-fold identities for the presence locations

     - \code{groupa} : k-fold identities for the background locations

Optionally the function also returns elements \code{block.object} and \code{speciesData}. These can be used to visualize data with \code{\link[blockCV]{foldExplorer}}.

}

\value{
The function returns a list with the following elements:. 

\item{k}{ A list with data on folds that can be directly used by \code{\link[BiodiversityR]{ensemble.calibrate.weights}}. }
\item{block.object}{ the results of \code{\link[blockCV]{spatialBlock}} or \code{\link[blockCV]{envBlock}} }
\item{speciesData}{ a \code{\link[sp]{SpatialPointsDataFrame}} with species data}

}

\references{
Roberts et al., 2017. Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. Ecography. 40: 913-929. 
}


\author{ Roeland Kindt (World Agroforestry Centre) }


\examples{
\dontrun{

library(blockCV)
library(sf)

# get predictor variables
library(dismo)
predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
    pattern='grd', full.names=TRUE)
predictors <- stack(predictor.files)
# subset based on Variance Inflation Factors
predictors <- subset(predictors, subset=c("bio5", "bio6", 
    "bio16", "bio17"))
predictors
predictors@title <- "base"

# presence points
presence_file <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
pres <- read.table(presence_file, header=TRUE, sep=',')[, -1]

# choose background points
background <- randomPoints(predictors, n=1000, p=pres, excludep=T, extf=1.00)
background <- data.frame(background)
colnames(background)=c('lon', 'lat')

# spatial blocking with square blocks of 1000 km and minimum 20 points in each categor
# fails if EPSG is not assigned
block.data <- ensemble.spatialBlock(x=predictors, p=pres, a=background, 
    EPSG=NULL,
    showBlocks=F, theRange=1000000, k=4, numLimit=20, iteration=1000, return.object=T)

block.data <- ensemble.spatialBlock(x=predictors, p=pres, a=background, 
    EPSG=4326,
    showBlocks=F, theRange=1000000, k=4, numLimit=20, iteration=1000, return.object=T)

# explore the results
foldExplorer(blocks=block.data$block.object, rasterLayer=predictors, 
    speciesData=block.data$speciesData)

# apply in calibration of ensemble weights
# make sure that folds apply to subset of points
p.spatial <- block.data$k$p 
a.spatial <- block.data$k$a
k.spatial <- block.data$k

ensemble.w1 <- ensemble.calibrate.weights(x=predictors, 
    p=p.spatial, a=a.spatial, k=k.spatial, 
    species.name="Bradypus",
    SINK=FALSE, PROBIT=TRUE,
    MAXENT=0, MAXNET=1, MAXLIKE=1, GBM=1, GBMSTEP=0, RF=0, CF=1,
    GLM=1, GLMSTEP=0, GAM=1, GAMSTEP=0, MGCV=0, MGCVFIX=0, 
    EARTH=0, RPART=0, NNET=1, FDA=0, SVM=0, SVME=0, GLMNET=0,
    BIOCLIM.O=1, BIOCLIM=1, DOMAIN=0, MAHAL=0, MAHAL01=0,
    ENSEMBLE.tune=TRUE,
    ENSEMBLE.best=0, ENSEMBLE.exponent=c(1, 2, 3),
    ENSEMBLE.min=0.7,
    Yweights="BIOMOD", 
    formulae.defaults=TRUE)

# confirm that correct folds were used
all.equal(ensemble.w1$groupp, block.data$k$groupp)
all.equal(ensemble.w1$groupa, block.data$k$groupa)

# environmental blocking with minimum 5 points in each category
block.data2 <- ensemble.envBlock(x=predictors, p=pres, a=background, 
    factors="biome",
    k=4, numLimit=5, return.object=T)

# explore the results
foldExplorer(blocks=block.data2$block.object, rasterLayer=predictors, 
    speciesData=block.data2$speciesData)

# apply in calibration of ensemble weights
# make sure that folds apply to subset of points
p.env <- block.data2$k$p 
a.env <- block.data2$k$a
k.env <- block.data2$k

ensemble.w2 <- ensemble.calibrate.weights(x=predictors, 
    p=p.env, a=a.env, k=k.env, 
    species.name="Bradypus",
    SINK=FALSE, PROBIT=TRUE,
    MAXENT=0, MAXNET=1, MAXLIKE=1, GBM=1, GBMSTEP=0, RF=0, CF=1,
    GLM=1, GLMSTEP=0, GAM=1, GAMSTEP=0, MGCV=0, MGCVFIX=0, 
    EARTH=0, RPART=0, NNET=1, FDA=0, SVM=0, SVME=0, GLMNET=0,
    BIOCLIM.O=1, BIOCLIM=1, DOMAIN=0, MAHAL=0, MAHAL01=0,
    ENSEMBLE.tune=TRUE,
    ENSEMBLE.best=0, ENSEMBLE.exponent=c(1, 2, 3),
    ENSEMBLE.min=0.7,
    factors="biome",
    Yweights="BIOMOD", 
    formulae.defaults=TRUE)

# confirm that correct folds were used
all.equal(ensemble.w2$groupp, block.data2$k$groupp)
all.equal(ensemble.w2$groupa, block.data2$k$groupa)

}

}