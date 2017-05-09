\name{ensemble.raster}
\alias{ensemble.raster}
\alias{ensemble.habitat.change}
\alias{ensemble.area}

\title{
Suitability mapping based on ensembles of modelling algorithms: consensus mapping
}

\description{
The basic function \code{ensemble.raster} creates two consensus raster layers, one based on a (weighted) average of different suitability modelling algorithms, and a second one documenting the number of modelling algorithms that predict presence of the focal organisms. Modelling algorithms include maximum entropy (MAXENT), boosted regression trees, random forests, generalized linear models (including stepwise selection of explanatory variables), generalized additive models (including stepwise selection of explanatory variables), multivariate adaptive regression splines, regression trees, artificial neural networks, flexible discriminant analysis, support vector machines, the BIOCLIM algorithm, the DOMAIN algorithm and the Mahalonobis algorithm. These sets of functions were developed in parallel with the \code{biomod2} package, especially for inclusion of the maximum entropy algorithm, but also to allow for a more direct integration with the BiodiversityR package, more direct handling of model formulae and greater focus on mapping. Researchers and students of species distribution are strongly encouraged to familiarize themselves with all the options of the \code{biomod2} and \code{dismo} packages. 
}

\usage{
ensemble.raster(xn = NULL, 
    models.list = NULL, 
    input.weights = models.list$output.weights,
    thresholds = models.list$thresholds,
    RASTER.species.name = models.list$species.name, 
    RASTER.stack.name = xn@title, 
    RASTER.format = "raster", RASTER.datatype = "INT2S", RASTER.NAflag = -32767,
    RASTER.models.overwrite = TRUE,
    KML.out = FALSE, KML.maxpixels = 100000, KML.blur = 10,
    evaluate = FALSE, SINK = FALSE,
    p = models.list$p, a = models.list$a,
    pt = models.list$pt, at = models.list$at)

ensemble.habitat.change(base.map = file.choose(), 
    other.maps = utils::choose.files(),
    change.folder = "ensembles/change",
    RASTER.format = "raster", RASTER.datatype = "INT1U", RASTER.NAflag = 255,
    KML.out = FALSE, KML.folder = "kml/change",
    KML.maxpixels = 100000, KML.blur = 10)

ensemble.area(x=NULL, km2=TRUE)

}

\arguments{
  \item{xn}{ RasterStack object (\code{\link[raster]{stack}}) containing all layers that correspond to explanatory variables of an ensemble calibrated earlier with \code{\link{ensemble.calibrate.models}}. See also \code{\link[raster]{predict}}.   }
  \item{models.list}{ list with 'old' model objects such as \code{MAXENT} or \code{RF}. }
  \item{input.weights}{ array with numeric values for the different modelling algorithms; if \code{NULL} then values provided by parameters such as \code{MAXENT} and \code{GBM} will be used. As an alternative, the output from \code{ensemble.calibrate.weights} can be used. }
  \item{thresholds}{ array with the threshold values separating predicted presence for each of the algorithms. }
  \item{RASTER.species.name}{ First part of the names of the raster files that will be generated, expected to identify the modelled species (or organism). }
  \item{RASTER.stack.name}{ Last part of the names of the raster files that will be generated, expected to identify the predictor stack used. }
  \item{RASTER.format}{ Format of the raster files that will be generated. See \code{\link[raster]{writeFormats}} and \code{\link[raster]{writeRaster}}.  }
  \item{RASTER.datatype}{ Format of the raster files that will be generated. See \code{\link[raster]{dataType}} and \code{\link[raster]{writeRaster}}.  }
  \item{RASTER.NAflag}{ Value that is used to store missing data. See \code{\link[raster]{writeRaster}}.  }
  \item{RASTER.models.overwrite}{ Overwrite the raster files that correspond to each suitability modelling algorithm (if \code{TRUE}). (Overwriting actually implies that raster files are created or overwritten that start with "working_").  }
  \item{KML.out}{ if \code{FALSE}, then no kml layers (layers that can be shown in Google Earth) are produced. If \code{TRUE}, then kml files will be saved in a subfolder 'kml'.   }
  \item{KML.maxpixels}{ Maximum number of pixels for the PNG image that will be displayed in Google Earth. See also \code{\link[raster]{KML}}.  }
  \item{KML.blur}{ Integer that results in increasing the size of the PNG image by \code{KML.blur^2}, which may help avoid blurring of isolated pixels. See also \code{\link[raster]{KML}}.  }
  \item{evaluate}{ if \code{TRUE}, then evaluate the created raster layers at locations \code{p}, \code{a}, \code{pt} and \code{at} (if provided). See also \code{\link[dismo]{evaluate}} }
  \item{SINK}{ Append the results to a text file in subfolder 'outputs' (if \code{TRUE}). The name of file is based on argument \code{RASTER.species.name}. In case the file already exists, then results are appended. See also \code{\link[base]{sink}}.  }
  \item{p}{ presence points used for calibrating the suitability models, typically available in 2-column (x, y) or (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{a}{ background points used for calibrating the suitability models, typically available in 2-column (x, y) or (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{pt}{ presence points used for evaluating the suitability models, typically available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}}  }
  \item{at}{ background points used for calibrating the suitability models, typicall available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }

  \item{base.map}{ filename with baseline map used to produce maps that show changes in suitable habitat }
  \item{other.maps}{ files with other maps used to produce maps that show changes in suitable habitat from a defined base.map  }
  \item{change.folder}{ folder where new maps documenting changes in suitable habitat will be stored. NOTE: please ensure that the base folder (eg: ../ensembles) exists already.  }
  \item{KML.folder}{ folder where new maps (in Google Earth format) documenting changes in suitable habitat will be stored. NOTE: please ensure that the base folder (eg: ../kml) exists already.   }

  \item{x}{ RasterLayer object (\code{\link[raster]{raster}}) in a longitude-latitude coordinate system  }
  \item{km2}{ Provide results in square km rather than square m. See also \code{\link[geosphere]{areaPolygon}}  }
}

\details{
The basic function \code{ensemble.raster} fits individual suitability models for all models with positive input weights. In subfolder "models" of the working directory, suitability maps for the individual suitability modelling algorithms are stored. In subfolder "ensembles", a consensus suitability map based on a weighted average of individual suitability models is stored. In subfolder "ensembles/presence", a presence-absence (1-0) map will be provided. In subfolder "ensembles/count", a consensus suitability map based on the number of individual suitability models that predict presence of the focal organism is stored.

Several of the features of \code{ensemble.raster} are also available from \code{\link{ensemble.calibrate.models}}. The main difference between the two functions is that \code{ensemble.raster} generates raster layers for individual suitability models, whereas the purpose of \code{\link{ensemble.calibrate.models}} is specifically to test different suitability modelling algorithms.

Note that values in suitability maps are integer values that were calculated by multiplying probabilities by 1000 (see also \code{\link[base]{trunc}}).

As the Mahalanobis function (\code{\link[dismo]{mahal}}) does not always provide values within a range of 0 - 1, the output values are rescaled by first subtracting the value of \code{1 - MAHAL.shape} from each prediction, followed by calculating the absolute value, followed by calculating the reciprocal value and finally multiplying this reciprocal value with \code{MAHAL.shape}. As this rescaling method does not estimate probabilities, inclusion in the calculation of a (weighted) average of ensemble probabilities may be problematic (the same applies to other distance-based methods). 

The \code{ensemble.habitat.change} function produces new raster layers that show changes in suitable and not suitable habitat between a base raster and a list of other rasters. The output uses the following coding: 0 = areas that remain unsuitable, 11 = areas that remain suitable, 10 = areas of lost habitat, 1 = areas of new habitat. (Codes are inspired on a binary classification of habitat suitability in base [1- or 0-] and other layer [-1 or -0], eg new habitat is coded 01 = 1).

With \code{KML.out = TRUE}, kml files are created in a subfolder named "KML". The colouring of the consensus suitability PNG is based on 20 intervals of size 50 between 0 and 1000. The colouring of the presence-absence PNG uses green for presence and red for absence. The colouring of the count suitability PNG uses black for zero (no models predict presence) and blue for the theoretical maximum number of models to predict presence (i.e. the count of all final weights), whereas intermediate numbers (1 to theoretical maximum - 1) are ranged from red to green. The colouring of the habitat change maps are: black (cells that are never suitable [value: 0]), green (cells that are always suitable [value: 11]), red (cells that are lost habitat [value: 10] and blue (cells that are new habitat [value: 1]).

The \code{ensemble.area} function calculates the area of different categories with \code{\link[geosphere]{areaPolygon}}

}

\value{
The basic function \code{ensemble.raster} mainly results in the creation of raster layers that correspond to fitted probabilities of presence of individual suitability models (in folder "models") and consensus models (in folder "ensembles"), and the number of suitability models that predict presence (in folder "ensembles"). Prediction of presence is based on a threshold usually defined by maximizing the sum of the true presence and true absence rates (see \code{threshold.method} and also \code{\link[dismo]{ModelEvaluation}}).

If desired by the user, the \code{ensemble.raster} function also saves details of fitted suitability models or data that can be plotted with the \code{\link{evaluation.strip.plot}} function.

}

\references{
Buisson L, Thuiller W, Casajus N, Lek S and Grenouillet G. 2010. Uncertainty in ensemble forecasting of species distribution. Global Change Biology 16: 1145-1157
}



\author{ Roeland Kindt (World Agroforestry Centre), Eike Luedeling (World Agroforestry Centre) and Evert Thomas (Bioversity International) }


\seealso{
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

# presence points
# presence points
presence_file <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
pres <- read.table(presence_file, header=TRUE, sep=',')[,-1]

# choose background points
background <- randomPoints(predictors, n=1000, extf = 1.00)

# if desired, change working directory where subfolders of "models" and 
# "ensembles" will be created
# raster layers will be saved in subfolders of /models and /ensembles:
getwd()

# first calibrate the ensemble
# calibration is done in two steps
# in step 1, a k-fold procedure is used to determine the weights
# in step 2, models are calibrated for all presence and background locations
# factor is not used as it is not certain whether correct levels will be used
# it may therefore be better to use dummy variables

# step 1: determine weights through 4-fold cross-validation
ensemble.calibrate.step1 <- ensemble.calibrate.weights(
    x=predictors, p=pres, a=background, k=4, 
    SINK=TRUE, species.name="Bradypus",
    MAXENT=1, MAXLIKE=0, GBM=1, GBMSTEP=0, RF=0, GLM=1, GLMSTEP=0, 
    GAM=1, GAMSTEP=0, MGCV=1, MGCVFIX=0,
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, GLMNET=1,
    BIOCLIM.O=0, BIOCLIM=1, DOMAIN=1, MAHAL=0, MAHAL01=1,
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
ensemble.raster.results <- ensemble.raster(xn=predictors, 
    models.list=ensemble.calibrate.step2$models, 
    input.weights=model.weights,
    SINK=TRUE, evaluate=TRUE,
    RASTER.species.name="Bradypus", RASTER.stack.name="base")

# use the base map to check for changes in suitable habitat
# this type of analysis is typically done with different predictor layers
# (for example, predictor layers representing different possible future climates)
# In this example, changes from a previous model (ensemble.raster.results)
# are contrasted with a newly calibrated model (ensemble.raster.results2)
# step 1: 4-fold cross-validation
ensemble.calibrate2.step1 <- ensemble.calibrate.weights(
    x=x.batch, p=p.batch, a=a.batch, MAXENT.a=MAXENT.a.batch, 
    factors=factors.batch, dummy.vars=dummy.vars.batch, 
    k=4, 
    SINK=TRUE, species.name="Bradypus",
    MAXENT=1, MAXLIKE=0, GBM=1, GBMSTEP=0, RF=0, GLM=1, GLMSTEP=0, 
    GAM=1, GAMSTEP=0, MGCV=1, MGCVFIX=0,
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, GLMNET=1,
    BIOCLIM.O=0, BIOCLIM=1, DOMAIN=1, MAHAL=0, MAHAL01=1,
    ENSEMBLE.tune=TRUE, PROBIT=TRUE,
    ENSEMBLE.best=0, ENSEMBLE.exponent=c(1, 2, 3),
    ENSEMBLE.min=c(0.65, 0.7),
    Yweights="BIOMOD",
    PLOTS=FALSE, formulae.defaults=TRUE)

model.weights2 <- ensemble.calibrate2.step1$output.weights

ensemble.calibrate2.step2 <- ensemble.calibrate.models(
    x=x.batch, p=p.batch, a=a.batch, MAXENT.a=MAXENT.a.batch, 
    factors=factors.batch, dummy.vars=dummy.vars.batch, 
    SINK=TRUE, species.name="Bradypus",
    models.keep=TRUE,
    input.weights=model.weights2,
    ENSEMBLE.tune=FALSE, PROBIT=TRUE,
    Yweights="BIOMOD",
    PLOTS=FALSE, formulae.defaults=TRUE)

ensemble.raster.results2 <- ensemble.raster(
    xn=predictors, 
    models.list=ensemble.calibrate2.step2$models, 
    input.weights=model.weights2,
    SINK=TRUE, evaluate=TRUE,
    RASTER.species.name="Bradypus", RASTER.stack.name="recalibrated")

base.file <- paste(getwd(), "/ensembles/presence/Bradypus_base.grd", sep="")
other.file <- paste(getwd(), "/ensembles/presence/Bradypus_recalibrated.grd", sep="")

changed.habitat <- ensemble.habitat.change(base.map=base.file, 
    other.maps=c(other.file),
    change.folder="ensembles/change")

change.file <- paste(getwd(), "/ensembles/change/Bradypus_recalibrated_presence.grd", sep="")

par.old <- graphics::par(no.readonly=T)
dev.new()
par(mfrow=c(2,2))
raster::plot(raster(base.file), breaks=c(-1, 0, 1), col=c("grey", "green"), 
    legend.shrink=0.8, main="base presence")
raster::plot(raster(other.file), breaks=c(-1, 0, 1), col=c("grey", "green"), 
    legend.shrink=0.8, main="other presence")
raster::plot(raster(change.file), breaks=c(-1, 0, 1, 10, 11), 
    col=c("grey", "blue", "red", "green"), 
    legend.shrink=0.8, main="habitat change", sub="11 remaining, 10 lost, 1 new")
graphics::par(par.old)

areas <- ensemble.area(raster(change.file))
areas
}
}