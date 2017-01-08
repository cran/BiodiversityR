\name{evaluation.strip.data}
\alias{evaluation.strip.data}
\alias{evaluation.strip.plot}

\title{Evaluation strips for ensemble suitability mapping}

\description{
These functions provide a dataframe which can subsequently be used to evaluate the relationship between environmental variables and the fitted probability of occurrence of individual or ensemble suitability modelling algorithms. The \code{biomod2} package provides an alternative implementation of this approach (\code{response.plot2}). 
}

\usage{
evaluation.strip.data(xn = NULL, ext = NULL,
    models.list = NULL, 
    input.weights = models.list$output.weights,
    steps=200
)

evaluation.strip.plot(data, TrainData=NULL,
    variable.focal = NULL, model.focal = NULL, 
    dev.new.width = 7, dev.new.height = 7, ...
)

}

\arguments{
  \item{xn}{ RasterStack object (\code{\link[raster]{stack}}) containing all layers that correspond to explanatory variables of an ensemble calibrated earlier with \code{\link{ensemble.calibrate.models}}. See also \code{\link[raster]{predict}}.   }
  \item{ext}{ an Extent object to limit the prediction to a sub-region of \code{xn} and the selection of background points to a sub-region of \code{x}, typically provided as c(lonmin, lonmax, latmin, latmax); see also \code{\link[raster]{predict}}, \code{\link[dismo]{randomPoints}} and \code{\link[raster]{extent}} }
  \item{models.list}{ list with 'old' model objects such as \code{MAXENT} or \code{RF}. }
  \item{input.weights}{ array with numeric values for the different modelling algorithms; if \code{NULL} then values provided by parameters such as \code{MAXENT} and \code{GBM} will be used. As an alternative, the output from \code{ensemble.calibrate.weights} can be used. }
  \item{steps}{ number of steps within the range of a continuous explanatory variable }

  \item{data}{ data set with ranges of environmental variables and fitted suitability models, typically returned by \code{evaluation.strip.data} }
  \item{TrainData}{ Data set representing the calibration data set. If provided, then a boxplot will be added for presence locations via \code{\link[graphics]{boxplot}} }
  \item{variable.focal}{ focal explanatory variable for plots with evaluation strips }
  \item{model.focal}{ focal model for plots with evaluation strips }
  \item{dev.new.width}{ Width for new graphics device (\code{\link[grDevices]{dev.new}}). If < 0, then no new graphics device is opened.  }
  \item{dev.new.height}{ Heigth for new graphics device (\code{\link[grDevices]{dev.new}}). If < 0, then no new graphics device is opened.  }
  \item{...}{ Other arguments passed to \code{\link[graphics]{plot}} }
}

\details{
These functions are mainly intended to be used internally by the ensemble.raster function.

\code{evaluation.strip.data} creates a data frame with variables (columns) corresponding to the environmental variables encountered in the RasterStack object (\code{x}) and the suitability modelling approaches that were defined. The variable of \code{focal.var} is an index of the variable for which values are ranged. The variable of \code{categorical} is an index for categorical (factor) variables.

A continuous (numeric) variable is ranged between its minimum and maximum values in the number of steps defined by argument \code{steps}. When a continuous variable is not the focal variable, then the average (\code{\link[base]{mean}}) is used.

A categorical (factor) variable is ranged for all the encountered levels (\code{\link[base]{levels}}) for this variable. When a categorical variable is not the focal variable, then the most frequent level is used.
}

\value{
function \code{evaluation.strip.data} creates a data frame, function code{evaluation.strip.data} allows for plotting. 
}

\references{
Elith J, Ferrier S, Huettmann F & Leathwick J. 2005. The evaluation strip: A new and robust method for plotting predicted responses from species distribution models. Ecological Modelling 186: 280-289 
}

\author{ Roeland Kindt (World Agroforestry Centre) }

\seealso{
\code{\link{ensemble.calibrate.models}} and \code{\link{ensemble.raster}}
}

\examples{
\dontrun{

# get predictor variables
library(dismo)
predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
    pattern='grd', full.names=TRUE)
predictors <- stack(predictor.files)
# subset based on Variance Inflation Factors
predictors <- subset(predictors, subset=c("bio5", "bio6", 
    "bio16", "bio17"))
predictors <- stack(predictors)
predictors
predictors@title <- "base"

# presence points
presence_file <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
pres <- read.table(presence_file, header=TRUE, sep=',')[,-1]

# the kfold function randomly assigns data to groups; 
# groups are used as calibration (1/5) and training (4/5) data
groupp <- kfold(pres, 5)
pres_train <- pres[groupp !=  1, ]
pres_test <- pres[groupp ==  1, ]

# choose background points
background <- randomPoints(predictors, n=1000, extf=1.00)
colnames(background)=c('lon', 'lat')
groupa <- kfold(background, 5)
backg_train <- background[groupa != 1, ]
backg_test <- background[groupa == 1, ]

# calibrate the models
# MAXLIKE not included as does not allow predictions for data.frames
# ENSEMBLE.min and ENSEMBLE.weight.min set very low to explore all
# algorithms.
# If focus is on actual ensemble, then set ENSEMBLE.min and 
# ENSEMBLE.weight.min to more usual values
ensemble.calibrate <- ensemble.calibrate.models(x=predictors, 
    p=pres_train, a=backg_train, 
    pt=pres_test, at=backg_test,
    ENSEMBLE.min=0.5, ENSEMBLE.weight.min = 0.001,
    MAXENT=1, MAXLIKE=0, GBM=1, GBMSTEP=1, RF=1, GLM=1, GLMSTEP=1, GAM=1, 
    GAMSTEP=1, MGCV=1, MGCVFIX=1, EARTH=1, RPART=1, NNET=1, FDA=1, 
    SVM=1, SVME=1, GLMNET=1,
    BIOCLIM.O=1, BIOCLIM=1, DOMAIN=1, MAHAL=0, MAHAL01=1,
    Yweights="BIOMOD", 
    PLOTS=FALSE, models.keep=TRUE)

# obtain data for plotting the evaluation strip
strip.data <- evaluation.strip.data(xn=predictors, steps=500,
    models.list=ensemble.calibrate$models)

# in case predictions for DOMAIN failed
# however, ENSEMBLE should also be recalculated
DOMAIN.model <- ensemble.calibrate$models$DOMAIN
strip.data$plot.data[, "DOMAIN"] <- dismo::predict(object=DOMAIN.model, 
    x=strip.data$plot.data)

# in case predictions for MAHAL01 failed
predict.MAHAL01 <- function(model, newdata, MAHAL.shape) {
    p <- dismo::predict(object=model, x=newdata)
    p <- p - 1 - MAHAL.shape
    p <- abs(p)
    p <- MAHAL.shape / p
    return(as.numeric(p))
}

MAHAL01.model <- ensemble.calibrate$models$MAHAL01
MAHAL.shape1 <- ensemble.calibrate$models$formulae$MAHAL.shape
strip.data$plot.data[, "MAHAL01"] <- predict.MAHAL01(model=MAHAL01.model, 
    newdata=strip.data$plot.data, MAHAL.shape=MAHAL.shape1)

# create graphs
evaluation.strip.plot(data=strip.data$plot.data, variable.focal="bio6",
    TrainData=strip.data$TrainData,
    type="o", col="red")
evaluation.strip.plot(data=strip.data$plot.data, model.focal="ENSEMBLE",
    TrainData=strip.data$TrainData,
    type="o", col="red")

}
}
