\name{ensemble.bioclim.graph}
\alias{ensemble.bioclim.graph}
\alias{ensemble.bioclim.graph.data}

\title{Graphs of bioclimatic ranges of species and climates}

\description{
The main graph function makes graphs that show mean, median, minimum, maximum and lower and upper limits for species or climates. The \code{ensemble.bioclim.graph.data} function creates input data, using \code{\link[BiodiversityR]{ensemble.bioclim.object}} internally.
}

\usage{
ensemble.bioclim.graph(graph.data = NULL, focal.var = NULL, 
    species.climates.subset = NULL, cols = NULL,
    var.multiply = 1.0, ref.lines = TRUE)

ensemble.bioclim.graph.data(
    x=NULL, p=NULL, fraction = 0.9, 
    species.climate.name="Species001_base", factors = NULL)
}

\arguments{
  \item{graph.data}{ Input data with same variables as created by \code{\link[BiodiversityR]{ensemble.bioclim.graph}}  }
  \item{focal.var}{ Bioclimatic variable to be plotted in the graph }
  \item{species.climates.subset}{ Character vector with subset of names of species and climates to be plotted in the graph (if not provided, then all species and climates will be plotted).  }
  \item{cols}{ colours for the different species and climates }
  \item{var.multiply}{ multiplier for the values to be plotted; 0.1 should be used if the bioclimatic variable was multiplied by 10 in the raster layers as in WorldClim and AFRICLIM }
  \item{ref.lines}{ If \code{TRUE}, then horizontal reference lines will be added for the minimum and maximum values of the species or climate plotted on the extreme left in the graph }

  \item{x}{ RasterStack object (\code{\link[raster]{stack}}) containing all environmental layers for which statistics should be calculated; see also \code{\link[BiodiversityR]{ensemble.bioclim}}.  }
  \item{p}{ presence points used for calibrating the suitability models, typically available in 2-column (lon, lat) dataframe; see also \code{\link[BiodiversityR]{ensemble.bioclim}}.  }
  \item{fraction}{ Fraction of range representing the optimal limits, default value of 0.9 as in the original BIOCLIM software; see also \code{\link[BiodiversityR]{ensemble.bioclim}}.  }
  \item{species.climate.name}{ Name for the species or climate that will be used as label in the graph.  }
  \item{factors}{ vector that indicates which variables are factors; these variables will be ignored by the BIOCLIM algorithm; see also \code{\link[BiodiversityR]{ensemble.bioclim}}.   } 
}

\details{
The function creates a graph that shows mean, median, minimum, maximum and upper and lower limits for a range of species and climates. The graph can be useful in interpreting results of \code{\link[BiodiversityR]{ensemble.bioclim}} or \code{\link[BiodiversityR]{ensemble.novel}}.

In the graphs, means are indicated by an asterisk (\code{pch=8} and medians as larger circles (\code{pch=1})).
}

\value{
function \code{ensemble.bioclim.graph.data} creates a data frame, function \code{ensemble.bioclim.graph} allows for plotting. 
}

\author{ Roeland Kindt (World Agroforestry Centre) }

\seealso{
\code{\link[BiodiversityR]{ensemble.bioclim}}
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

# climates for north and south (use same process for future climates)
ext2 <- extent(-90, -32, 0, 23)
predictors2 <- crop(predictors, y=ext2)
predictors2 <- stack(predictors2)
predictors2@title <- "north"

ext3 <- extent(-90, -32, -33, 0)
predictors3 <- crop(predictors, y=ext3)
predictors3 <- stack(predictors3)
predictors3@title <- "south"

graph.data1 <- ensemble.bioclim.graph.data(predictors, p=pres, 
    factors="biome", species.climate.name="Bradypus")
graph.data2 <- ensemble.bioclim.graph.data(predictors, p=NULL, 
    factors="biome", species.climate.name="baseline")
graph.data3 <- ensemble.bioclim.graph.data(predictors2, p=NULL, 
    factors="biome", species.climate.name="north")
graph.data4 <- ensemble.bioclim.graph.data(predictors3, p=NULL, 
    factors="biome", species.climate.name="south")
graph.data.all <- rbind(graph.data1, graph.data2, graph.data3, graph.data4)

par.old <- graphics::par(no.readonly=T)
graphics::par(mfrow=c(2, 2))

ensemble.bioclim.graph(graph.data.all, focal.var="bio5", 
    var.multiply=0.1, cols=c("black", rep("blue", 3)))
ensemble.bioclim.graph(graph.data.all, focal.var="bio6", 
    var.multiply=0.1, cols=c("black", rep("blue", 3)))
ensemble.bioclim.graph(graph.data.all, focal.var="bio16", 
    var.multiply=1.0, cols=c("black", rep("blue", 3)))
ensemble.bioclim.graph(graph.data.all, focal.var="bio17", 
    var.multiply=1.0, cols=c("black", rep("blue", 3)))

graphics::par(par.old)

}
}
