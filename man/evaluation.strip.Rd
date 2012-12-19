\name{evaluation.strip.data}
\alias{evaluation.strip.data}
\alias{evaluation.strip.plot}

\title{Evaluation strips for ensemble suitability mapping}

\description{
These functions provide a dataframe which can subsequently be used to evaluate the relationship between environmental variables and the fitted probability of occurrence of individual or ensemble suitability modelling algorithms. The \code{biomod2} package provides an alternative implementation of this approach (\code{response.plot2}). 
}

\usage{
evaluation.strip.data(x, ext = NULL, vars=names(x), factors = NULL, steps = 50, 
    modelnames = c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", 
    "EARTH", "RPART", "NNET", "FDA", "SVM", "BIOCLIM", "DOMAIN", "MAHAL"))

evaluation.strip.plot(data, modelnames = c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "EARTH", "RPART", "NNET", "FDA", "SVM",
    "BIOCLIM", "DOMAIN", "MAHAL"), variable = NULL, model = NULL, ...)
}

\arguments{
  \item{x}{ RasterStack object (\code{\link[raster]{stack}}) containing all layers that correspond to explanatory variables used for suitability modelling  }
  \item{ext}{ an Extent object to crop the area of interest to a sub-region of \code{x}, typically provided as c(lonmin, lonmax, latmin, latmax); see also \code{\link[raster]{crop}} and \code{\link[raster]{extent}} }
  \item{vars}{ Vector that indicates which variables should be included as columns in the data frame. Only variables that correspond to layers of the rasterStack will be included. }
  \item{factors}{ vector that indicates which variables are factors; see also \code{\link[dismo]{prepareData}} }
  \item{steps}{ number of steps within the range of a continuous explanatory variable }
  \item{modelnames}{ abbreviated names of the individual suitability models that are fitted. Should correspond to the names used in ensemble.grd. }
  \item{data}{ data set with ranges of environmental variables and fitted suitability models, typically returned by ensemble.grd }
  \item{variable}{ focal explanatory variable for plots with evaluation strips }
  \item{model}{ focal explanatory variable for plots with evaluation strips }
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
\code{\link{ensemble.raster}}
}

\examples{
# get predictor variables
library(dismo)
predictors <- stack(list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), pattern='grd', full.names=TRUE ))
names(predictors)

# create data set with ranges of explanatory variables
new.data <- evaluation.strip.data(predictors, factors="biome")
summary(new.data[new.data[,"focal.var"]==1,])
}
