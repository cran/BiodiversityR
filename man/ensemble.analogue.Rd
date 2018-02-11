\name{ensemble.analogue}
\alias{ensemble.analogue}
\alias{ensemble.analogue.object}

\title{
Climate analogues from climatic distance raster layers.
}

\description{
Function \code{ensemble.analogue} creates the map with climatic distance and provides the locations of the climate analogues (defined as locations with smallest climatic distance to a reference climate). Function \code{ensemble.analogue.object} provides the reference values used by the prediction function used by \code{\link[raster]{predict}} . 
}

\usage{
ensemble.analogue(x = NULL, analogue.object = NULL, analogues = 1,
    RASTER.object.name = analogue.object$name, RASTER.stack.name = x@title,
    RASTER.format = "raster", RASTER.datatype = "INT2S", RASTER.NAflag = -32767,
    KML.out = T, KML.blur = 10, KML.maxpixels = 100000, 
    limits = c(1, 5, 20, 50), limit.colours = c('red', 'orange', 'blue', 'grey'),
    CATCH.OFF = FALSE)

ensemble.analogue.object(ref.location, future.stack, current.stack, name = "reference1", 
    method = "mahal", an = 10000, probs = c(0.025, 0.975), weights = NULL, z = 2)

}

\arguments{
  \item{x}{ RasterStack object (\code{\link[raster]{stack}}) containing all environmental layers (climatic variables) for which climatic distance should be calculated. }
  \item{analogue.object}{ Object listing reference values for the environmental layers and additional parameters (covariance matrix for \code{method = "mahal"} or normalization parameters for \code{method = "quantile"}) that are used by the prediction function that is used internally by \code{\link[raster]{predict}}. This object is created with \code{\link{ensemble.analogue.object}}. }
  \item{analogues}{ Number of analogue locations to be provided }
  \item{RASTER.object.name}{ First part of the names of the raster file that will be generated, expected to identify the area and time period for which ranges were calculated }
  \item{RASTER.stack.name}{ Last part of the names of the raster file that will be generated, expected to identify the predictor stack used }
  \item{RASTER.format}{ Format of the raster files that will be generated. See \code{\link[raster]{writeFormats}} and \code{\link[raster]{writeRaster}}.  }
  \item{RASTER.datatype}{ Format of the raster files that will be generated. See \code{\link[raster]{dataType}} and \code{\link[raster]{writeRaster}}.  }
  \item{RASTER.NAflag}{ Value that is used to store missing data. See \code{\link[raster]{writeRaster}}.  }
  \item{KML.out}{ If \code{TRUE}, then kml files will be saved in a subfolder 'kml/zones'.  }
  \item{KML.maxpixels}{ Maximum number of pixels for the PNG image that will be displayed in Google Earth. See also \code{\link[raster]{KML}}.  }
  \item{KML.blur}{ Integer that results in increasing the size of the PNG image by \code{KML.blur^2}, which may help avoid blurring of isolated pixels. See also \code{\link[raster]{KML}}.  }
  \item{limits}{ Limits indicating the accumulated number of closest analogue sites. These limits will correspond to different colours in the KML map. In the default setting, the closest analogue will be coloured red and the second to fifth closest analogues will be coloured orange.  }
  \item{limit.colours}{ Colours for the different limits based on number of analogues.  }
  \item{CATCH.OFF}{ Disable calls to function \code{\link[base]{tryCatch}}.  }

  \item{ref.location}{ Location of the reference location for which analogues are searched for and from which climatic distance will be calculated, typically available in 2-column (lon, lat) dataframe; see also \code{\link[raster]{extract}}. }
  \item{future.stack}{ RasterStack object (\code{\link[raster]{stack}}) containing the environmental layers (climatic variables) to obtain the conditions of the reference location. For climate change research, this RasterStack object corresponds to the future climatic conditions of the reference location. }
  \item{current.stack}{ RasterStack object (\code{\link[raster]{stack}}) containing all environmental layers (climatic variables) for which climatic distance should be calculated. For climate change research, this RasterStack object corresponds to the current climatic conditions and range where climate analogues are searched for. }
  \item{name}{ Name of the object, expect to expected to identify the area and time period for which ranges were calculated and where no novel conditions will be detected }
  \item{method}{ Method used to calculate climatic distance: \code{method = "mahal"} results in using the Mahalanobis distance (\code{\link[stats]{mahalanobis}}); \code{method = "quantile"} results in dividing the differences between reference climatic values and climatic values in the 'current' raster by a quantile range obtained from the 'current' raster; \code{method = "sd"} results in dividing the differences between reference climatic values and climatic values in the 'current' raster by standard deviations obtained from the 'current' raster; and \code{method = "none"} results in not dividing these differences. }
  \item{an}{ Number of randomly selected locations points to calculate the covariance matrix (\code{\link[stats]{cov}}) to be used with \code{\link[stats]{mahalanobis}}, therefore only used for \code{method = "mahal"}. See also \code{\link[dismo]{randomPoints}}.  }
  \item{probs}{ Numeric vector of probabilities [0,1] as used by \code{\link[stats]{quantile}}). Only used for \code{method = "quantile"}.  } 
  \item{weights}{ Numeric vector of weights by which each variable (difference) should be multiplied by (can be used to give equal weight to 12 monthly rainfall values and 24 minimum and maximum monthly temperature values). Not used for \code{method = "mahal"}.  } 
  \item{z}{ Parameter used as exponent for differences calculated between reference climatic variables and variables in the 'current' raster and reciprocal exponent for the sum of all differences. Default value of 2 corresponds to the Euclidean distance. Not used for \code{method = "mahal"}.  } 
}

\details{
Function \code{ensemble.analogues} maps the climatic distance from reference values determined by \code{ensemble.analogues.object} and provides the locations of the \code{analogues} closest analogues.

The \code{method = "mahal"} uses the Mahalanobis distance as environmental (climatic) distance: \code{\link[stats]{mahalanobis}}.

Other methods use a normalization method to handle scale differences between environmental (climatic) variables:

\eqn{ClimaticDistance = (\sum_i(weight_{i}*(|T_{i}-C_{i}| / norm_{i})^z))^(1/z)}{ClimaticDistance = (sum(weight[i]*(abs(T[i]-C[i])/norm[i])^z))^(1/z)}

where \eqn{T_{i}}{T[i]} are the target values for environmental (climatic) variable i,
\eqn{C_{i}}{C[i]} are the values in the current environmental layers where analogues are searched for,
\eqn{weight_{i}}{weight[i]} are the weights for environmental variable i, and
\eqn{norm_{i}}{norm[i]} are the normalization parameters for environmental variable i

}

\value{
Function \code{ensemble.analogue.object} returns a list with following objects:
\item{name }{ name for the reference location }
\item{ref.location }{ coordinates of the reference location }
\item{stack.name }{ name for time period for which values are extracted from the \code{future.stack} }
\item{method }{ method used for calculating climatic distance }
\item{target.values }{ target environmental values to select analogues for through minimum climatic distance }
\item{cov.mahal }{ covariance matrix }
\item{norm.values }{ parameters by which each difference between target and 'current' value will be divided }
\item{weight.values }{ weights by which each difference between target and 'current' value will be multiplied }
\item{z }{ parameter to be used as exponent for differences between target and 'current' values }
}

\author{ Roeland Kindt (World Agroforestry Centre) and Eike Luedeling (World Agroforestry Centre)}

\references{
Bos, Swen PM, et al. "Climate analogs for agricultural impact projection and adaptation-a reliability test." Frontiers in Environmental Science 3 (2015): 65.
Luedeling, Eike, and Henry Neufeldt. "Carbon sequestration potential of parkland agroforestry in the Sahel." Climatic Change 115.3-4 (2012): 443-461.
}

\seealso{
\code{\link{ensemble.novel}}
}

\examples{

\dontrun{
# get predictor variables
library(dismo)
predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
    pattern='grd', full.names=TRUE)
predictors <- stack(predictor.files)
predictors <- subset(predictors, subset=c("bio1", "bio5", "bio6", "bio7", "bio8", 
    "bio12", "bio16", "bio17"))
predictors
predictors@title <- "base"

# instead of searching for current analogue of future climate conditions,
# search for analogue in southern hemisphere
future.stack <- stack(crop(predictors, y=extent(-125, -32, 0, 40)))
future.stack@title <- "north"
current.stack <- stack(crop(predictors, y=extent(-125, -32, -56, 0)))
current.stack@title <- "south"

# reference location in Florida
# in this case future.stack and current.stack are both current
ref.loc <- data.frame(t(c(-80.19, 25.76)))
names(ref.loc) <- c("lon", "lat")

# climate analogue analysis based on the Mahalanobis distance
Florida.object.mahal <- ensemble.analogue.object(ref.location=ref.loc, 
    future.stack=future.stack, current.stack=current.stack, 
    name="FloridaMahal", method="mahal", an=10000)
Florida.object.mahal

Florida.analogue.mahal <- ensemble.analogue(x=current.stack, 
    analogue.object=Florida.object.mahal, analogues=50)
Florida.analogue.mahal

# climate analogue analysis based on the Euclidean distance and dividing each variable by the sd
Florida.object.sd <- ensemble.analogue.object(ref.location=ref.loc, 
    future.stack=future.stack, current.stack=current.stack, 
    name="FloridaSD", method="sd", z=2)
Florida.object.sd

Florida.analogue.sd <- ensemble.analogue(x=current.stack, 
    analogue.object=Florida.object.sd, analogues=50)
Florida.analogue.sd

# plot analogues on climatic distance maps
par(mfrow=c(1,2))
analogue.file <- paste(getwd(), "//ensembles//analogue//FloridaMahal_south_analogue.grd", sep="")
plot(raster(analogue.file), main="Mahalanobis climatic distance")
points(Florida.analogue.sd[3:50, "lat"] ~ Florida.analogue.sd[3:50, "lon"], 
    pch=1, col="red", cex=1)
points(Florida.analogue.mahal[3:50, "lat"] ~ Florida.analogue.mahal[3:50, "lon"], 
    pch=3, col="black", cex=1)
points(Florida.analogue.mahal[2, "lat"] ~ Florida.analogue.mahal[2, "lon"], 
    pch=22, col="blue", cex=2)
legend(x="topright", legend=c("closest", "Mahalanobis", "SD"), pch=c(22, 3 , 1), 
    col=c("blue" , "black", "red"))

analogue.file <- paste(getwd(), "//ensembles//analogue//FloridaSD_south_analogue.grd", sep="")
plot(raster(analogue.file), main="Climatic distance normalized by standard deviation")
points(Florida.analogue.mahal[3:50, "lat"] ~ Florida.analogue.mahal[3:50, "lon"], 
    pch=3, col="black", cex=1)
points(Florida.analogue.sd[3:50, "lat"] ~ Florida.analogue.sd[3:50, "lon"], 
    pch=1, col="red", cex=1)
points(Florida.analogue.sd[2, "lat"] ~ Florida.analogue.sd[2, "lon"], 
    pch=22, col="blue", cex=2)
legend(x="topright", legend=c("closest", "Mahalanobis", "SD"), pch=c(22, 3 , 1), 
    col=c("blue" , "black", "red"))
par(mfrow=c(1,1))
}

}

