\name{ensemble.PET.season}
\alias{ensemble.PET.season}

\title{
Calculate the balance between precipitation and potential evapotranspiration for the dry season with the largest balance (maximum climatological water deficit, accumulated aridity).
}

\description{
Internally, the function first determines different dry seasons, defined by consecutive months where precipitation is smaller than potential evapotranspiration. The function then returns the summation of monthly balances of precipitation minus potential evapotranspiration that is largest (most negative) of the different dry seasons. 
}

\usage{

ensemble.PET.season(PREC.stack = NULL, PET.stack = NULL, 
    filename = NULL, overwrite = TRUE,
    CATCH.OFF = FALSE, ...)

}

\arguments{
  \item{PREC.stack}{ stack object (\code{\link[raster]{stack}}) with monthly precipitation values.  }
  \item{PET.stack}{ stack object (\code{\link[raster]{stack}}) with monthly potential evapotranspiration values.  }
  \item{filename}{ Name for writing the resulting raster layer (as in \code{\link[raster]{writeRaster}}). }
  \item{overwrite}{ Replace a previous version of the same file.  }
  \item{CATCH.OFF}{ Disable calls to function \code{\link[base]{tryCatch}}.  }
  \item{...}{ Additional arguments for \code{\link[raster]{writeRaster}}.  }

}

\details{
Unlike the methodology described by Chave et al. (2014), the assumption is not made that there is a single drought season. Internally, the function first identifies dry months as months where the balance of precipitation minus potential evapotranspiration is negative. Then dry seasons are identified as consecutive dry months. For each dry season, the total sum of balances is calculated. The function finally identifies and returns the largest of these balances.

The algorithm of the function should obtain the same values of the Maximum Cumulative Water Deficit as from rules described by Aragao et al. 2007 (section 2.2), when using fixed monthly PET values of 100 mm instead of calculated monthly PET values (calculated, for example, from monthly mean temperatures and extraterrestrial solar radiation through the Hargreaves method).

Note that calculation may take a while for larger raster data sets.
}

\value{
The function returns and writes a raster layer
}

\references{
Chave J et al. 2014. Improved allometric models to estimate the aboveground biomass of tropical trees. Global Change Biology 20: 3177-3190.

Aragao LZ et al. 2007. Spatial patterns and fire response of recent Amazonian droughts. Geophysical Research Letters 34 L07701
}

\author{ Roeland Kindt (World Agroforestry Centre)}

\seealso{
\code{\link{ensemble.batch}}
}

\examples{

\dontrun{

## Not run: 

library(raster)
stack1 <- stack(monthly.prec.files)
stack2 <- stack(monthly.PET.files)
# note that the stacks should be of the same extend and resolution
ensemble.PET.season(PREC.stack=stack1, PET.stack=stack2,
    filename=paste(getwd(), '//Aridity.deficit.tif', sep=""))

}
}

