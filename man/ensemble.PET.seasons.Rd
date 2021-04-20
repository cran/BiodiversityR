\name{ensemble.PET.seasons}
\alias{ensemble.PET.seasons}
\alias{ensemble.prec.season}
\alias{ensemble.tmean.season}
\alias{ensemble.season.suitability}

\title{
Raster calculations of beginnings and lengths of growing seasons from the difference between precipitation (P) and potential evapotranspiration (PET), defining dry months with 2 * P < PET.
}

\description{
The main function of \code{ensemble.PET.seasons} calculates the number of growing seasons and their starts and lengths from the dry period criterion of 2 * P < PET (\url{http://www.fao.org/nr/climpag/cropfor/lgp_en.asp}). Functions \code{ensemble.PREC.season} and \code{ensemble.TMEAN.season} calculate the total precipitation and average temperature for provided starts and lengths of a growing season. Together with data on optimal and absolute precipitation and temperature limits for a selected crop (as available from FAO's ECOCROP database), these layers enable the calculation of crop suitability using methods detailed in Chapman et al. (2020).
}

\usage{

ensemble.PET.seasons(PREC.stack=NULL, PET.stack=NULL, 
    index=c("seasons", "start1", "length1", "start2", "length2", "start3", "length3"),
    filename=NULL, overwrite=TRUE,
    CATCH.OFF=FALSE, ...)

ensemble.prec.season(PREC.stack=NULL, 
    start.layer=NULL, length.layer=NULL,
    filename=NULL, overwrite=TRUE,
    CATCH.OFF=FALSE, ...)

ensemble.tmean.season(TMEAN.stack=NULL, 
    start.layer=NULL, length.layer=NULL,
    filename=NULL, overwrite=TRUE,
    CATCH.OFF=FALSE, ...)

ensemble.season.suitability(season.raster=NULL, 
    thresholds=NULL,
    filename=NULL, overwrite=TRUE,
    CATCH.OFF=FALSE, ...)
    
}

\arguments{
  \item{PREC.stack}{ stack object (\code{\link[raster]{stack}}) with monthly precipitation values.  }
  \item{PET.stack}{ stack object (\code{\link[raster]{stack}}) with monthly potential evapotranspiration values.  }
  \item{TMEAN.stack}{ stack object (\code{\link[raster]{stack}}) with monthly average temperature values.  }
  \item{index}{ selection of type of output - see details.  }
  \item{start.layer}{ raster layer with index of the month of the start of the growing season.  }
  \item{length.layer}{ raster layer with index of the length of the growing season.  }
  \item{season.raster}{ raster layer with seasonal precipitation or mean temperature.  }
  \item{thresholds}{ optimal and absolute thresholds of crop suitability, defined similarly as by ECOCROP.  }
  \item{filename}{ Name for writing the resulting raster layer (as in \code{\link[raster]{writeRaster}}). }
  \item{overwrite}{ Replace a previous version of the same file.  }
  \item{CATCH.OFF}{ Disable calls to function \code{\link[base]{tryCatch}}.  }
  \item{...}{ Additional arguments for \code{\link[raster]{writeRaster}}.  }

}

\details{
Function \code{ensemble.PET.seasons} calculates the number, starts and lengths of growing seasons after first internally determining dry periods from the criterion of 2 * P < PET. The function was developed with data sets with monthly precipitatin and PET values, but probably can also work with data sets of other temporal resolution. Where there are multiple gaps between dry seasons, different growing periods are identified.

The definition of dry periods is less strict than the definition of P < PET used in \code{ensemble.PET.season}, following the methodologies for this function. 

Argument \code{index} determines the contents of the output rasters:
- \code{seasons} selects the number of growing periods to be returned;
- \code{start1} selects the index of the start of the first or only growing period to be returned;
- \code{length1} selects the index of the end of the first or only growing period to be returned;
- \code{start2} selects the index of the start of the second growing period to be returned;
- \code{length2} selects the index of the end of the second growing period to be returned;
- \code{start3} selects the index of the start of the third growing period to be returned; and
- \code{length3} selects the index of the end of the third growing period to be returned.

The methodology of calculating crop suitability is directly based on Chapman et al. (2020), following their equations 2 (temperature suitability, based on the mean temperature of the growing season) and 3 (precipitation suitability, based on the total precipitation of the growing season). The methods of Chapman et al. (2020) are based on Ramirez-Villegas et al. (2013), including the calculation of crop suitability as the product of temperature suitability and crop suitability (their respective equations 1 and 3).

Crop thresholds are available from the FAO ECOCROP database, which are also available via function \code{\link[dismo]{getCrop}}.

Note that calculations can take a while for larger data sets.
}

\value{
The function returns and writes raster layers.
}

\references{

Ramirez-Villegas J, Jarvis A and Laderach P. 2013. Empirical approaches for assessing impacts of climate change on agriculture: The EcoCrop model and a case study with grain sorghum. Agricultural and Forest Meteorology \doi{10.1016/j.agrformet.2011.09.005}

Chapman et al. 2020. Impact of climate change on crop suitability in sub-Saharan Africa in parameterized and convection-permitting regional climate models. Environmental Research Letters \url{https://iopscience.iop.org/article/10.1088/1748-9326/ab9daf}

}

\author{ Roeland Kindt (World Agroforestry Centre)}

\seealso{
\code{\link{ensemble.PET.season}}
}

\examples{

\dontrun{

## Not run: 

library(raster)
P.stack <- stack(monthly.prec.files)
PE.stack <- stack(monthly.PET.files)

# Calculate average monthly values similarly as in 
TMIN.stack <- stack(monthly.tmin.files)
TMAX.stack <- stack(monthly.tmax.files)
T.stack <- stack(0.5*(TMIN.stack + TMAX.stack))

# step 1: determine number of seasons, start and length of season 1

seasons.raster <- ensemble.PET.seasons(PREC.stack=P.stack, PET.stack=PE.stack, 
    index="seasons", filename="seasons.tif", CATCH.OFF=TRUE)

start1.raster <- ensemble.PET.seasons(PREC.stack=P.stack, PET.stack=PE.stack, 
    index="start1", filename="start1.tif", CATCH.OFF=TRUE)

length1.raster <- ensemble.PET.seasons(PREC.stack=P.stack, PET.stack=PE.stack, 
    index="length1", filename="length1.tif", CATCH.OFF=TRUE)

start2.raster <- ensemble.PET.seasons(PREC.stack=P.stack, PET.stack=PE.stack, 
    index="start2", filename="start2.tif", CATCH.OFF=TRUE)

length2.raster <- ensemble.PET.seasons(PREC.stack=P.stack, PET.stack=PE.stack, 
    index="length2", filename="length2.tif", CATCH.OFF=TRUE)

# step 2: calculate total precipitation in first rainy season,
# then use this value to calculate precipitation suitability

prec.season <- ensemble.prec.season(PREC.stack=P.stack, 
    start.layer=start1.raster, length.layer=length1.raster, 
    filename="precSeason.tif", CATCH.OFF=FALSE)  

dismo::getCrop("Sorghum (med. altitude)")

prec.suit <- ensemble.season.suitability(season.raster=prec.season, 
    thresholds=c(300, 500, 1000, 3000), 
    filename="precSuitability.tif", CATCH.OFF=FALSE)  

# step 3: calculate average temperature in first rainy steason,
# then use this value to calculate temperature suitability

tmean.season <- ensemble.tmean.season(TMEAN.stack=T.stack, 
    start.layer=start1.raster, length.layer=length1.raster, 
    filename="tmeanSeason.tif", CATCH.OFF=FALSE) 

temp.suit <- ensemble.season.suitability(season.raster=tmean.season, 
    thresholds=c(10, 24, 35, 40), 
    filename="tempSuitability.tif", CATCH.OFF=FALSE) 
    
# step 4: seasonal crop suitability is product of precipitation suitability
# and temperature suitability
    
sorghum.suit <- prec.suit * temp.suit
plot(sorghum.suit)

}
}

