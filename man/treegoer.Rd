\name{treegoer.score}
\alias{treegoer.score}
\alias{treegoer.widen}
\alias{treegoer.filter}
\alias{treegoer.position}
\alias{treegoer.map}

\title{
Calculate climate scores with the Tree Globally Observed Environmental Ranges (TreeGOER) database.
}

\description{
Function \code{\link[BiodiversityR]{treegoer.score}} calculates a climate score via a similar algorithm that is used internally in the GlobalUsefulNativeTrees (GlobUNT) database (Kindt et al. 2023, \doi{10.1038/s41598-023-39552-1}). The function depends on \code{\link[BiodiversityR]{treegoer.filter}} and requires a data set (argument \code{treegoer.wide}) as created from the Tree Globally Observed Environmental Ranges (TreeGOER) database (Kindt 2023, \doi{10.1111/gcb.16914}) via \code{\link[BiodiversityR]{treegoer.widen}}.
}

\usage{

treegoer.score(site.data, 
  site.species = treegoer.wide$species,
  treegoer.wide, 
  filter.vars = c("bio05", "bio14", "climaticMoistureIndex"),
  upper.only.vars = NULL,
  lower.only.vars = NULL)

treegoer.filter(site.data,
  treegoer.wide, 
  filter.vars = c("bio05", "bio14", "climaticMoistureIndex"),
  upper.only.vars = NULL,
  lower.only.vars = NULL,
  limit.vars = c("Q05", "Q95")) 

treegoer.widen(treegoer,
  species = unique(treegoer$species)[1:100],
  filter.vars = c("bio05", "bio14", "climaticMoistureIndex"))
  
treegoer.position(site.data, 
  treegoer.wide, 
  focal.var  = "bio01")

treegoer.map(map.rast, 
  map.species=treegoer[1, "species"],
  treegoer, 
  filter.vars=c("bio05", "bio14", "climaticMoistureIndex"),
  upper.only.vars = NULL,
  lower.only.vars = NULL,
  verbose=FALSE)


}

\arguments{
  \item{site.data}{ Data set with 1 row, containing the environmental conditions at the planting site for the selected environmental variables of the TreeGOER database. This data set can be set by selecting a city from the CitiesGOER datase (\url{https://zenodo.org/records/10004594}) or a weather station  from the ClimateForecasts database (\url{https://zenodo.org/records/10726088}). }
  \item{site.species}{ Species for which the climate score will be calculated.  }
  \item{treegoer.wide}{ Data set created by \code{\link[BiodiversityR]{treegoer.widen}} from the TreeGOER database, or another data set with the same variables.  }
  \item{filter.vars}{ Environmental variables for which ranges (minimum, maximum and 0.05, 0.25, 0.75 and 0.95 quantile) are documented in the treegoer.wide data set.  }
  \item{limit.vars}{ Selection of the lower and upper limits for the environmental ranges, typically set as \code{c("MIN", "MAX")} (marginal bioclimatic species domain as in the BIOCLIM algorithm; see Booth 2018, \doi{10.1111/aec.12628}), \code{c("Q05", "Q95")} (core of the domain) or \code{c("QRT1", "QRT3")} (middle of the domain).  }
  \item{upper.only.vars}{ Selection of variables that will only be checked at the upper limits.  }
  \item{lower.only.vars}{ Selection of variables that will only be checked at the lower limits.  }
  \item{treegoer}{ Data set with environmental limits that was locally downloaded file (\code{TreeGOER_2023.txt}) that was downloaded from the archive (\url{https://zenodo.org/records/10008994}).  }
  \item{species}{ Selection of species to document in the wide format.  }
  \item{focal.var}{ Selection of variable to calculate the position of planting site in environmental space.  }
  \item{map.rast}{ SpatRaster object (\code{\link[terra]{rast}}) with layers showing the environmental conditions of planting sites.}
  \item{map.species}{ Species selected for mapping.}
  \item{verbose}{ Report progress on the creation of the suitability map (when \code{TRUE}).}
}

\details{

The calculation of the climate score uses an expanded version of the algorithms used by BIOCLIM (Booth 2018, \doi{10.1111/aec.12628}).

- A score of 3 indicates that for all selected variables, the planting site has environmental conditions that are within the middle (0.25 to 0.75 quantiles) of the species range.

- A score of 2 indicates that for all selected variables, the planting site has environmental conditions that are within the core (0.05 to 0.95 quantiles) of the species range. For some variables, the planting conditions are outside the middle of the species range.

- A score of 1 indicates that for all selected variables, the planting site has environmental conditions that are within the documented limits (minimum to maximum) of the species range. For some variables, the planting conditions are outside the core of the species range; the BIOCLIM algorithm defines this domain as the 'marginal domain'.

- A score of 0 indicates that for some of the selected variables, the planting site has environmental conditions that are outside the documented limits (< minimum or > maximum) of the species range. 

- A score of 0.5 indicates that for some of the selected variables, the planting site has environmental conditions that are outside the documented limits (< minimum or > maximum) of the species range. However, for none of the selected variables the planting has environmental conditions larger than those for variables that are checked only at the lower side, or smaller than those for variables that are checked only at the upper side.

- A score of -1 indicates that there was no information on the environmental ranges of the species. 

The calculation of the position of the planting site via \code{treegoer.position} is done as follows:

- For sites where the conditions of the planting location (PL) are above the median, the position is calculated as:

    (PL - MEDIAN) / (MAX - MEDIAN)

- For sites where the conditions of the planting location (PL) are below the median, the position is calculated as:

    (PL - MEDIAN) / (MEDIAN - MIN)

The sign (positive or negative) will therefore indicate the position (respectively, lower or higher) of the planting site with respect to the median of the species.

The magnitude of the metric will indicate the relative distance of the planting site with respect to the difference between median and extremes. Values that are lower than -1 will indicate novel conditions below the minimum. Values that are above +1 will indicate novel conditions above the maximum.

The same algorithms and similar scripts are used internally in the GlobalUsefulNativeTrees database (see Kindt et al. 2023). The internal scripts also resemble scripts provided here: \url{https://rpubs.com/Roeland-KINDT/1114902}.

Function \code{treegoer.map} creates a raster layer with the climate scores.

}

\value{
Function \code{treegoer.score} returns a data set that includes a climate score representing the position of the planting site within the environmental range of species documented by the Tree Globally Observed Environmental Ranges database. 
}

\author{ Roeland Kindt (World Agroforestry, CIFOR-ICRAF)}

\references{
Booth TH. 2018. Why understanding the pioneering and continuing contributions of BIOCLIM to species distribution modelling is important. Austral Ecology 43: 852-860.
\doi{10.1111/aec.12628}

Kindt R. 2023. TreeGOER: A database with globally observed environmental ranges for 48,129 tree species. Global Change Biology. 
\doi{10.1111/gcb.16914}

Kindt R., Graudal L, Lilleso J.P.-B. et al. 2023. GlobalUsefulNativeTrees, a database documenting 14,014 tree species, supports synergies between biodiversity recovery and local livelihoods in landscape restoration. Scientific Reports. 
\doi{10.1038/s41598-023-39552-1}

Kindt R. 2023. Using the Tree Globally Observed Environmental Ranges and CitiesGOER databases to Filter GlobalUsefulNativeTrees Species lists.
\url{https://rpubs.com/Roeland-KINDT/1114902}

Kindt R. 2023. CitiesGOER: Globally Observed Environmental Data for 52,602 Cities with a Population >= 5000 (version 2023.10).
\doi{10.5281/zenodo.10004594}

Kindt R. 2024. ClimateForecasts: Globally Observed Environmental Data for 15,504 Weather Station Locations (version 2024.03).
\doi{10.5281/zenodo.10776414}
}

\examples{

\dontrun{

# Example adapted from https://rpubs.com/Roeland-KINDT/1114902

# treegoer.file <- choose.files() 
# Provide the location where the TreeGOER file was downloaded locally
# (https://zenodo.org/records/10008994: TreeGOER_2023.txt)
treegoer <- fread(treegoer.file, sep="|", encoding="UTF-8")
nrow(treegoer)
length(unique(treegoer$species)) # 48129

# A data set of tree species
# Example has useful tree species filtered 
# for Kenya and human food from the GlobalUsefulNativeTrees database
# (https://worldagroforestry.org/output/globalusefulnativetrees)
# globunt.file <- choose.files()
globunt <- fread(globunt.file, sep="|", encoding="UTF-8")
nrow(globunt) # 461

# Environmental variables used for filtering or scoring species
focal.vars <- c("bio01", "bio12",
                "climaticMoistureIndex", "monthCountByTemp10", 
                "growingDegDays5",
                "bio05", "bio06", "bio16", "bio17",
                "MCWD")

# Use treegoer.widen()
treegoer.wide <- treegoer.widen(treegoer=treegoer,
                                species=globunt$Switchboard,
                                filter.vars=focal.vars)
names(treegoer.wide)

# Environmental conditions at the planting site
# Provide the locations where the CitiesGOER files were downloaded locally 
# (https://zenodo.org/records/10004594: CitiesGOER_baseline.xlsx). 
# Alternatively, the ClimateForecasts database can be used
# (https://zenodo.org/records/10726088: ClimateForecasts_baseline.xlsx)
# baseline.file <- choose.files()
site.baseline <- data.frame(read_excel(baseline.file,
                            sheet="Cities data",
                            skip=6))
# Set the planting location in Nairobi
site.planting <- site.baseline[site.baseline$Name == "Nairobi", ]
site.planting

# Calculate the climate scores
treegoer.scores <- treegoer.score(site.species=globunt$Switchboard,
                                  treegoer.wide=treegoer.wide,
                                  filter.vars=focal.vars, 
                                  site.data=site.planting)

# Calculate the climate score for a single environmental variable
treegoer.score <- treegoer.score(site.species=globunt$Switchboard,
                                  treegoer.wide=treegoer.wide,
                                  filter.vars="bio01", 
                                  site.data=site.planting)
                                  
                                  
library(dismo)
predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
    pattern='grd', full.names=TRUE)
predictors <- rast(predictor.files)
# subset based on Variance Inflation Factors
predictors <- subset(predictors, subset=c("bio5", "bio6", 
    "bio16", "bio17"))
predictors

names(predictors)[1:2] <- c("bio05", "bio06")
# WorldClim1 had temperature values multiplied by 10
if(t(minmax(predictors[["bio05"]]))[1] > 25) {
  predictors[["bio05"]] <- predictors[["bio05"]]/10
  predictors[["bio06"]] <- predictors[["bio06"]]/10
}

map.test2 <- treegoer.map(map.rast=predictors,
              map.species="Bertholletia excelsa",
             treegoer=treegoer,
             filter.vars=c("bio05", "bio06", "bio16", "bio17"),
             upper.only.vars=c("bio05", "bio06"),
             lower.only.vars=c("bio16", "bio17"))

map.test2
plot(map.test2)
                                  
}

}