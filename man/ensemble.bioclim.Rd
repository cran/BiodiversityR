\name{ensemble.bioclim}
\alias{ensemble.bioclim}
\alias{ensemble.bioclim.object}

\title{
Suitability mapping based on the BIOCLIM algorithm
}

\description{
Implementation of the BIOCLIM algorithm more similar to the original BIOCLIM algorithm and software than the implementation through \code{\link[dismo]{bioclim}}. Function \code{ensemble.bioclim} creates the suitability map. Function \code{ensemble.bioclim.object} provides the reference values used by the prediction function used by \code{\link[raster]{predict}} . 
}

\usage{
ensemble.bioclim(x = NULL, bioclim.object = NULL, 
    RASTER.object.name = bioclim.object$species.name, RASTER.stack.name = x@title,
    RASTER.format = "GTiff",
    CATCH.OFF = FALSE)

ensemble.bioclim.object(x = NULL, p = NULL, fraction = 0.9,
    quantiles = TRUE, 
    species.name = "Species001", 
    factors = NULL)
}

\arguments{
  \item{x}{ RasterStack object (\code{\link[raster]{stack}}) containing all environmental layers for which suitability should be calculated, or alternatively a data.frame containing the bioclimatic variables.  }
  \item{bioclim.object}{ Object listing optimal and absolute minima and maxima for bioclimatic variables, used by the prediction function that is used internally by \code{\link[raster]{predict}}. This object is created with \code{\link{ensemble.bioclim.object}}. }
  \item{RASTER.object.name}{ First part of the names of the raster file that will be generated, expected to identify the species or crop for which ranges were calculated }
  \item{RASTER.stack.name}{ Last part of the names of the raster file that will be generated, expected to identify the predictor stack used }
  \item{RASTER.format}{ Format of the raster files that will be generated. See \code{\link[raster]{writeFormats}} and \code{\link[raster]{writeRaster}}.  }
  \item{CATCH.OFF}{ Disable calls to function \code{\link[base]{tryCatch}}.  }

  \item{p}{ presence points used for calibrating the suitability models, typically available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}.  }
  \item{fraction}{ Fraction of range representing the optimal limits, default value of 0.9 as in the original BIOCLIM software.  }
  \item{quantiles}{ If \code{TRUE} then optimal limits are calculated as quantiles corresponding to \code{0.5-fraction/2} and \code{0.5+fraction/2} percentiles. If \code{FALSE} then optimal limits are calculated from the normal distribution with \code{mean - cutoff*sd} and \code{mean + cutoff*sd} with \code{cutoff} calculated as \code{qnorm(0.5+fraction/2)}. }
  \item{species.name}{ Name by which the model results will be saved.  }
  \item{factors}{ vector that indicates which variables are factors; these variables will be ignored by the BIOCLIM algorithm   } 
}

\details{
Function \code{ensemble.bioclim} maps suitability for a species based on optimal (percentiles, typically 5 and 95 percent) and absolute (minimum to maximum) limits for bioclimatic variables. If all values at a given location are within the optimal limits, suitability values are mapped as 1 (suitable). If not all values are within the optimal limits, but all values are within the absolute limits, suitability values are mapped as 0.5 (marginal). If not all values are within the absolute limits, suitability values are mapped as 0 (unsuitable).

Function \code{ensemble.bioclim.object} calculates the optimal and absolute limits. Optimal limits are calculated based on the parameter \code{fraction}, resulting in optimal limits that correspond to 0.5-fraction/2 and 0.5+fraction/2 (the default value of 0.9 therefore gives a lower limit of 0.05 and a upper limit of 0.95). Two methods are implemented to obtain optimal limits for the lower and upper limits. One method (\code{quantiles = FALSE}) uses mean, standard deviation and a cutoff parameter calculated with \code{\link[stats]{qnorm}}. The other method (\code{quantiles = TRUE}) calculates optimal limits via the \code{\link[stats]{quantile}} function. To handle possible asymmetrical distributions better, the second method is used as default.

When \code{x} is a RasterStack and point locations are provided, then optimal and absolute limits correspond to the bioclimatic values observed for the locations. When \code{x} is RasterStack and point locations are not provided, then optimal and absolute limits correspond to the bioclimatic values of the RasterStack. 

Applying to algorithm without providing point locations will provide results that are similar to the \code{\link[BiodiversityR]{ensemble.novel}} function, whereby areas plotted as not suitable will be the same areas that are novel.
}

\value{
Function \code{ensemble.bioclim.object} returns a list with following objects:
\item{lower.limits }{ vector with lower limits for each bioclimatic variable }
\item{upper.limits }{ vector with upper limits for each bioclimatic variable }
\item{minima }{ vector with minima for each bioclimatic variable }
\item{maxima }{ vector with maxima for each bioclimatic variable }
\item{means }{ vector with mean values for each bioclimatic variable }
\item{medians }{ vector with median values for each bioclimatic variable }
\item{sds }{ vector with standard deviation values for each bioclimatic variable }
\item{cutoff }{ cutoff value for the normal distribution }
\item{fraction }{ fraction of values within the optimal limits }
\item{species.name }{ name for the species }
}

\references{
Nix HA. 1986. A biogeographic analysis of Australian elapid snakes. In: Atlas of Elapid Snakes of Australia. (Ed.) R. Longmore, pp. 4-15. Australian Flora and Fauna Series Number 7. Australian Government Publishing Service: Canberra.

Booth TH, Nix HA, Busby JR and Hutchinson MF. 2014. BIOCLIM: the first species distribution modelling package, its early applications and relevance to most current MAXENT studies. Diversity and Distributions 20: 1-9
}

\author{ Roeland Kindt (World Agroforestry Centre) with inputs from Trevor Booth (CSIRO)}

\seealso{
\code{\link[dismo]{bioclim}}, \code{\link[BiodiversityR]{ensemble.bioclim.graph}} and \code{\link[BiodiversityR]{ensemble.novel}}
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
    "bio16", "bio17", "biome"))
predictors
predictors@title <- "base"

# presence points
presence_file <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
pres <- read.table(presence_file, header=TRUE, sep=',')[,-1]

background <- dismo::randomPoints(predictors, n=100)
colnames(background)=c('lon', 'lat')

pres.dataset <- data.frame(extract(predictors, y=pres))
names(pres.dataset) <- names(predictors)
pres.dataset$biome <- as.factor(pres.dataset$biome)

Bradypus.bioclim <- ensemble.bioclim.object(predictors, quantiles=T, 
    p=pres, factors="biome", species.name="Bradypus")
Bradypus.bioclim
# obtain the same results with a data.frame
Bradypus.bioclim2 <- ensemble.bioclim.object(pres.dataset, quantiles=T, 
    species.name="Bradypus")
Bradypus.bioclim2
# obtain results for entire rasterStack
Bradypus.bioclim3 <- ensemble.bioclim.object(predictors, p=NULL, quantiles=T, 
    factors="biome", species.name="America")
Bradypus.bioclim3

ensemble.bioclim(x=predictors, bioclim.object=Bradypus.bioclim)
ensemble.bioclim(x=predictors, bioclim.object=Bradypus.bioclim3)

par.old <- graphics::par(no.readonly=T)
graphics::par(mfrow=c(1,2))

rasterfull1 <- paste("ensembles//Bradypus_base_BIOCLIM_orig.tif", sep="")
raster::plot(raster(rasterfull1), breaks=c(-0.1, 0, 0.5, 1), 
    col=c("grey", "blue", "green"), main="original method")
rasterfull2 <- paste("ensembles//America_base_BIOCLIM_orig.tif", sep="")
raster::plot(raster(rasterfull2), breaks=c(-0.1, 0, 0.5, 1), 
    col=c("grey", "blue", "green"), main="America")

graphics::par(par.old)

# compare with implementation bioclim in dismo
bioclim.dismo <- bioclim(predictors, p=pres)
rasterfull2 <- paste("ensembles//Bradypus_base_BIOCLIM_dismo.tif", sep="")
raster::predict(object=predictors, model=bioclim.dismo, na.rm=TRUE, 
    filename=rasterfull2, progress='text', overwrite=TRUE)

par.old <- graphics::par(no.readonly=T)
graphics::par(mfrow=c(1,2))

raster::plot(raster(rasterfull1), breaks=c(-0.1, 0, 0.5, 1), 
    col=c("grey", "blue", "green"), main="original method")
raster::plot(raster(rasterfull2), main="dismo method")

graphics::par(par.old)

# use dummy variables to deal with factors
predictors <- stack(predictor.files)
biome.layer <- predictors[["biome"]]
biome.layer
ensemble.dummy.variables(xcat=biome.layer, most.frequent=0, freq.min=1,
    overwrite=TRUE)

predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
    pattern='grd', full.names=TRUE)
predictors <- stack(predictor.files)
predictors.dummy <- subset(predictors, subset=c("biome_1", "biome_2",  "biome_3",  
    "biome_4", "biome_5", "biome_7",  "biome_8",  "biome_9", "biome_10", 
    "biome_12", "biome_13", "biome_14"))
predictors.dummy
predictors.dummy@title <- "base_dummy"

Bradypus.dummy <- ensemble.bioclim.object(predictors.dummy, quantiles=T, 
    p=pres, species.name="Bradypus")
Bradypus.dummy
ensemble.bioclim(x=predictors.dummy, bioclim.object=Bradypus.dummy)

par.old <- graphics::par(no.readonly=T)
graphics::par(mfrow=c(1,2))

rasterfull3 <- paste("ensembles//Bradypus_base_dummy_BIOCLIM_orig.tif", sep="")
raster::plot(raster(rasterfull1), breaks=c(-0.1, 0, 0.5, 1), col=c("grey", "blue", "green"), 
    main="numeric predictors")
raster::plot(raster(rasterfull3), breaks=c(-0.1, 0, 0.5, 1), col=c("grey", "blue", "green"), 
    main="dummy predictors")

graphics::par(par.old)
}

}

