\name{ensemble.terra}
\alias{ensemble.terra}

\title{
Suitability mapping based on ensembles of modelling algorithms: consensus mapping via the terra package
}

\description{
The function \code{ensemble.terra} creates two consensus raster layers, one based on a (weighted) average of different suitability modelling algorithms, and a second one documenting the number of modelling algorithms that predict presence of the focal organisms. This function has the same behaviour as \code{\link[BiodiversityR]{ensemble.raster}}. 
}

\usage{
ensemble.terra(xn = NULL, 
    models.list = NULL, 
    input.weights = models.list$output.weights,
    thresholds = models.list$thresholds,
    RASTER.species.name = models.list$species.name, 
    RASTER.stack.name = "xnTitle", 
    RASTER.filetype = "GTiff", RASTER.datatype = "INT2S", RASTER.NAflag = -32767,
    RASTER.models.overwrite = TRUE,
    evaluate = FALSE, SINK = FALSE,
    p = models.list$p, a = models.list$a,
    pt = models.list$pt, at = models.list$at,
    CATCH.OFF = FALSE)
}

\arguments{
  \item{xn}{ SpatRaster object (\code{\link[terra]{rast}}) containing all layers that correspond to explanatory variables of an ensemble calibrated earlier with \code{\link{ensemble.calibrate.models}}. See also \code{\link[terra]{predict}}.   }
  \item{models.list}{ list with 'old' model objects such as \code{MAXENT} or \code{RF}. }
  \item{input.weights}{ array with numeric values for the different modelling algorithms; if \code{NULL} then values provided by parameters such as \code{MAXENT} and \code{GBM} will be used. As an alternative, the output from \code{ensemble.calibrate.weights} can be used. }
  \item{thresholds}{ array with the threshold values separating predicted presence for each of the algorithms. }
  \item{RASTER.species.name}{ First part of the names of the raster files that will be generated, expected to identify the modelled species (or organism). }
  \item{RASTER.stack.name}{ Last part of the names of the raster files that will be generated, expected to identify the predictor stack used. }
  \item{RASTER.filetype}{ Format of the raster files that will be generated. See  \code{\link[terra]{writeRaster}}.  }
  \item{RASTER.datatype}{ Format of the raster files that will be generated. See  \code{\link[terra]{writeRaster}}.  }
  \item{RASTER.NAflag}{ Value that is used to store missing data. See \code{\link[terra]{writeRaster}}.  }
  \item{RASTER.models.overwrite}{ Overwrite the raster files that correspond to each suitability modelling algorithm (if \code{TRUE}). (Overwriting actually implies that raster files are created or overwritten that start with "working_").  }
  \item{evaluate}{ if \code{TRUE}, then evaluate the created raster layers at locations \code{p}, \code{a}, \code{pt} and \code{at} (if provided). See also \code{\link[dismo]{evaluate}} }
  \item{SINK}{ Append the results to a text file in subfolder 'outputs' (if \code{TRUE}). The name of file is based on argument \code{RASTER.species.name}. In case the file already exists, then results are appended. See also \code{\link[base]{sink}}.  }
  \item{p}{ presence points used for calibrating the suitability models, typically available in 2-column (x, y) or (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{a}{ background points used for calibrating the suitability models, typically available in 2-column (x, y) or (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{pt}{ presence points used for evaluating the suitability models, typically available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}}  }
  \item{at}{ background points used for calibrating the suitability models, typicall available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{CATCH.OFF}{ Disable calls to function \code{\link[base]{tryCatch}}.  }
}

\details{
The basic function \code{ensemble.terra} fits individual suitability models for all models with positive input weights. In subfolder "models" of the working directory, suitability maps for the individual suitability modelling algorithms are stored. In subfolder "ensembles", a consensus suitability map based on a weighted average of individual suitability models is stored. In subfolder "ensembles/presence", a presence-absence (1-0) map will be provided. In subfolder "ensembles/count", a consensus suitability map based on the number of individual suitability models that predict presence of the focal organism is stored.

Note that values in suitability maps are integer values that were calculated by multiplying probabilities by 1000 (see also \code{\link[base]{trunc}}).
}

\value{
The basic function \code{ensemble.terra} mainly results in the creation of raster layers that correspond to fitted probabilities of presence of individual suitability models (in folder "models") and consensus models (in folder "ensembles"), and the number of suitability models that predict presence (in folder "ensembles"). Prediction of presence is based on a threshold usually defined by maximizing the sum of the true presence and true absence rates (see \code{threshold.method} and also \code{\link[dismo]{ModelEvaluation}}).
}

\references{
Kindt R. 2018. Ensemble species distribution modelling with transformed suitability values. Environmental Modelling & Software 100: 136-145. \doi{10.1016/j.envsoft.2017.11.009}

Buisson L, Thuiller W, Casajus N, Lek S and Grenouillet G. 2010. Uncertainty in ensemble forecasting of species distribution. Global Change Biology 16: 1145-1157
}

\author{ Roeland Kindt (World Agroforestry Centre)}

\seealso{
\code{\link{ensemble.raster}},
\code{\link{evaluation.strip.plot}}, \code{\link{ensemble.calibrate.models}}, \code{\link{ensemble.calibrate.weights}}, \code{\link{ensemble.batch}}
}

\examples{
\dontrun{
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
predictors@title <- "base"

# make a SpatRaster object
# Ideally this should not be created from files in the 'raster' grd format
# (so a better method would be to create instead from 'tif' files).

predictors.terra <- terra::rast(predictors)
# predictors@title <- "base"
crs(predictors.terra) <- c("+proj=longlat +datum=WGS84")
predictors.terra

# presence points
presence_file <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
pres <- read.table(presence_file, header=TRUE, sep=',')[,-1]

# choose background points
background <- dismo::randomPoints(predictors, n=1000, extf = 1.00)

# if desired, change working directory where subfolders of "models" and 
# "ensembles" will be created
# raster layers will be saved in subfolders of /models and /ensembles:
getwd()

# first calibrate the ensemble
# calibration is done in two steps
# in step 1, a k-fold procedure is used to determine the weights
# in step 2, models are calibrated for all presence and background locations

# Although a spatRaster (predictors.terra) object is used as input for 
# ensemble.calibrate.weights and ensemble.calibrate.models,
# internally the spatRaster will be converted to a rasterStack for these
# functions (among other things, to allow for dismo::prepareData)

# step 1: determine weights through 4-fold cross-validation
ensemble.calibrate.step1 <- ensemble.calibrate.weights(
    x=predictors.terra, p=pres, a=background, k=4, 
    SINK=TRUE, species.name="Bradypus",
    MAXENT=0, MAXNET=1, MAXLIKE=1, GBM=1, GBMSTEP=0, RF=1, CF=1,
    GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=1, MGCV=1, MGCVFIX=1, 
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, GLMNET=1,
    BIOCLIM.O=1, BIOCLIM=1, DOMAIN=1, MAHAL=0, MAHAL01=1,
    ENSEMBLE.tune=TRUE, PROBIT=TRUE,
    ENSEMBLE.best=0, ENSEMBLE.exponent=c(1, 2, 3),
    ENSEMBLE.min=c(0.65, 0.7),
    Yweights="BIOMOD",
    PLOTS=FALSE, formulae.defaults=TRUE)

# step 1 generated the weights for each algorithm
model.weights <- ensemble.calibrate.step1$output.weights
x.batch <- ensemble.calibrate.step1$x
p.batch <- ensemble.calibrate.step1$p
a.batch <- ensemble.calibrate.step1$a
MAXENT.a.batch <- ensemble.calibrate.step1$MAXENT.a
factors.batch <- ensemble.calibrate.step1$factors
dummy.vars.batch <- ensemble.calibrate.step1$dummy.vars

# step 2: calibrate models with all available presence locations
# weights determined in step 1 calculate ensemble in step 2
ensemble.calibrate.step2 <- ensemble.calibrate.models(
    x=x.batch, p=p.batch, a=a.batch, MAXENT.a=MAXENT.a.batch, 
    factors=factors.batch, dummy.vars=dummy.vars.batch, 
    SINK=TRUE, species.name="Bradypus",
    models.keep=TRUE,
    input.weights=model.weights,
    ENSEMBLE.tune=FALSE, PROBIT=TRUE,
    Yweights="BIOMOD",
    PLOTS=FALSE, formulae.defaults=TRUE)

# step 3: use previously calibrated models to create ensemble raster layers
# re-evaluate the created maps at presence and background locations
# (note that re-evaluation will be different due to truncation of raster layers
# as they wered saved as integer values ranged 0 to 1000)
ensemble.terra.results <- ensemble.terra(xn=predictors.terra, 
    models.list=ensemble.calibrate.step2$models, 
    input.weights=model.weights,
    SINK=TRUE, evaluate=TRUE,
    RASTER.species.name="Bradypus", RASTER.stack.name="base")

}
}