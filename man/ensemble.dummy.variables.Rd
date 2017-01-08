\name{ensemble.dummy.variables}
\alias{ensemble.dummy.variables}
\alias{ensemble.accepted.categories}
\alias{ensemble.simplified.categories}

\title{
Suitability mapping based on ensembles of modelling algorithms: handling of categorical data
}

\description{
The basic function \code{ensemble.dummy.variables} creates new raster layers representing dummy variables (coded 0 or 1) for all or the most frequent levels of a caterogical variable. Sometimes the creation of dummy variables is needed for proper handling of categorical data for some of the suitability modelling algorithms.
}

\usage{
ensemble.dummy.variables(xcat = NULL, 
    freq.min = 50, most.frequent = 5,
    new.levels = NULL, overwrite = TRUE, ...)

ensemble.accepted.categories(xcat = NULL, categories = NULL, 
    filename = NULL, overwrite = TRUE, ...)

ensemble.simplified.categories(xcat = NULL, p = NULL, 
    filename = NULL, overwrite = TRUE, ...)
}

\arguments{
  \item{xcat}{ RasterLayer object (\code{\link[raster]{raster}}) containing values for a categorical explanatory variable.  }
  \item{freq.min}{ Minimum frequency for a dummy raster layer to be created for the corresponding factor level. See also \code{\link[raster]{freq}}.  }
  \item{most.frequent}{ Number of dummy raster layers to be created (if larger than 0), corresponding to the same number of most frequent factor levels  See also \code{\link[raster]{freq}}.  }
  \item{new.levels}{ character vector giving factor levels that are not encountered in \code{xcat} and for which dummy layers should be created (this could be useful in dealing with novel conditions)  }
  \item{overwrite}{ overwrite an existing file name with the same name (if \code{TRUE}). See also \code{\link[raster]{writeRaster}}.  }
  \item{...}{ additional arguments for \code{\link[raster]{writeRaster}} or (for \code{ensemble.dummy.variables}, \code{\link[raster]{writeRaster}}).  }

  \item{categories}{ numeric vector providing the accepted levels of a categorical raster layer; expected to correspond to the levels encountered during calibration }
  \item{filename}{ name for the output file. See also \code{\link[raster]{writeRaster}}.  }

  \item{p}{ presence points that will be used for calibrating the suitability models, typically available in 2-column (x, y) or (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
}

\details{
The basic function \code{ensemble.dummy.variables} creates dummy variables from a \code{RasterLayer} object (see \code{\link[raster]{raster}}) that represents a categorical variable. With \code{freq.min} and \code{most.frequent} it is possible to limit the number of dummy variables that will be created. For example, \code{most.frequent = 5} results in five dummy variables to be created.  

Function \code{ensemble.accepted.categories} modifies the \code{RasterLayer} object (see \code{\link[raster]{raster}}) by replacing cell values for categories (levels) that are not accepted with missing values.

Function \code{ensemble.simplified.categories} modifies the \code{RasterLayer} object (see \code{\link[raster]{raster}}) by replacing cell values for categories (levels) where none of the presence points occur with the same level. This new level is coded by the maximum coding level for these 'outside categories'. 
}

\value{
The basic function \code{ensemble.dummy.variables} mainly results in the creation of raster layers that correspond to dummy variables.
}


\author{ Roeland Kindt (World Agroforestry Centre) and Evert Thomas (Bioversity International) }


\seealso{
\code{\link{ensemble.calibrate.models}}, \code{\link{ensemble.raster}}
}

\examples{
\dontrun{

# get predictor variables
library(dismo)
predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
    pattern='grd', full.names=TRUE)
predictors <- stack(predictor.files)
biome.layer <- predictors[["biome"]]
biome.layer

# create dummy layers for the 5 most frequent factor levels

ensemble.dummy.variables(xcat=biome.layer, most.frequent=5,
    overwrite=TRUE)

# check whether dummy variables were created
predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
    pattern='grd', full.names=TRUE)
predictors <- stack(predictor.files)
predictors
names(predictors)

# once dummy variables were created, avoid using the original categorical data layer
predictors <- subset(predictors, subset=c("bio5", "bio6", "bio16", "bio17", 
    "biome_1", "biome_2", "biome_7", "biome_8", "biome_13"))
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

# note that dummy variables with no variation are not used by DOMAIN
# note that dummy variables are not used by MAHAL and MAHAL01
# (neither are categorical variables)
ensemble.nofactors <- ensemble.calibrate.models(x=predictors, p=pres_train, a=backg_train, 
    pt=pres_test, at=backg_test,
    species.name="Bradypus",
    VIF=T,
    MAXENT=1, MAXLIKE=1, GBM=1, GBMSTEP=0, RF=1, GLM=1, GLMSTEP=0, GAM=1, 
    GAMSTEP=0, MGCV=1, MGCVFIX=0, EARTH=1, RPART=1, NNET=1, FDA=1, 
    SVM=1, SVME=1, BIOCLIM.O=1, BIOCLIM=1, DOMAIN=1, MAHAL=0, MAHAL01=1,
    Yweights="BIOMOD", 
    dummy.vars=c("biome_1", "biome_2", "biome_7", "biome_8", "biome_13"),
    PLOTS=FALSE, evaluations.keep=TRUE)
}

}