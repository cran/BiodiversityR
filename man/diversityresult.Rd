\name{diversityresult}
\alias{diversityresult}
\alias{diversitycomp}
\alias{diversityvariables}

\title{Alternative Diversity Results}

\description{
Provides alternative methods of obtaining results on diversity statistics than provided directly by functions \code{\link[vegan]{diversity}}, \code{\link[vegan:diversity]{fisher.alpha}}, \code{\link[vegan]{specpool}} and \code{\link[vegan:diversity]{specnumber}} (all from \pkg{vegan}), although these same functions are called. Some other statistics are also calculated such as the reciprocal Berger-Parker diversity index and abundance (not a diversity statistic). The statistics can be calculated for the entire community, for each site separately, the mean of the sites can be calculated or a jackknife estimate can be calculated for the community.
}

\usage{
diversityresult(x, y = NULL, factor = NULL, level = NULL, 
    index=c("Shannon", "Simpson", "inverseSimpson", "Logalpha", "Berger", 
        "richness", "abundance", "Jevenness", "Eevenness", 
        "jack1", "jack2", "chao", "boot"),
    method=c("pooled", "each site", "mean", "sd", "max", "jackknife"), 
    sortit = FALSE, digits = 8)

diversityvariables(x, y, digits=8)

diversitycomp(x, y = NULL, 
    factor1 = NULL ,factor2 = NULL,
    index=c("Shannon", "Simpson", "inverseSimpson", "Logalpha", "Berger", 
        "richness", "abundance", "Jevenness", "Eevenness", 
        "jack1", "jack2", "chao", "boot"),
    method=c("pooled", "mean", "sd", "max", "jackknife"),
    sortit=FALSE, digits=8)
}

\arguments{
  \item{x}{ Community data frame with sites as rows, species as columns and species abundance as cell values. }
  \item{y}{ Environmental data frame. }
  \item{factor}{ Variable of the environmental data frame that defines subsets to calculate diversity statistics for. }
  \item{level}{ Level of the variable to create the subset to calculate diversity statistics. }
  \item{index}{ Type of diversity statistic with "richness" to calculate species richness, "abundance" to calculate abundance, "Shannon" to calculate the Shannon diversity index, "Simpson" to calculate 1-Simpson concentration index, "inverseSimpson" to calculate the reciprocal Simpson diversity index, "Logalpha" to calculate the log series alpha diversity index, "Berger" to calculate the reciprocal Berger-Parker diversity index, "Jevenness" to calculate one Shannon evenness index, "Eevenness" to calculate another Shannon evenness index, "jack1" to calculate the first-order jackknife gamma diversity estimator, "jack2" to calculate the second-order jackknife gamma diversity estimator, "chao" to calculate the Chao gamma diversity estimator and "boot" to calculate the bootstrap gamma diversity estimator. }
  \item{method}{ Method of calculating the diversity statistics: "pooled" calculates the diversity of the entire community (all sites pooled), "each site" calculates diversity for each site separetly, "mean" calculates the average diversity of the sites, "sd" calculates the standard deviation of the diversity of the sites, "max" calculates the maximum diversity of the sites, whereas "jackknife" calculates the jackknifed diversity for the entire data frame. }
  \item{sortit}{ Sort the sites by increasing values of the diversity statistic. }
  \item{digits}{ Number of digits in the results. }

  \item{factor1}{ Variable of the environmental data frame that defines subsets to calculate diversity statistics for. }
  \item{factor2}{ Optional second variable of the environmental data frame that defines subsets to calculate diversity statistics for in a crosstabulation with the other variable of the environmental data frame. }
}

\details{
These functions provide some alternative methods of obtaining results with diversity statistics, although functions \code{\link[vegan]{diversity}}, \code{\link[vegan:diversity]{fisher.alpha}}, \code{\link[vegan]{specpool}}, \code{\link[vegan:specpool]{estimateR}} and \code{\link[vegan:diversity]{specnumber}} (all from \pkg{vegan}) are called to calculate the various statistics.

Function \code{diversityvariables} adds variables to the environmental dataset (richness, Shannon, Simpson, inverseSimpson, Logalpha, Berger, Jevenness, Eevenness).

The reciprocal Berger-Parker diversity index is the reciprocal of the proportional abundance of the most dominant species.

J-evenness is calculated as: H / ln(S) where H is the Shannon diversity index and S the species richness.

E-evenness is calculated as: exp(H) / S where H is the Shannon diversity index and S the species richness.

The method of calculating the diversity statistics include following options: "all" calculates the diversity of the entire community (all sites pooled together), "s" calculates the diversity of each site separatedly, "mean" calculates the average diversity of the sites, whereas "Jackknife" calculates the jackknifed diversity for the entire data frame. Methods "s" and "mean" are not available for function diversitycomp. Gamma diversity estimators assume that the method is "all".

Functions \code{diversityresult} and \code{diversitycomp} allow to calculate diversity statistics for subsets of the community and environmental data sets. Function \code{diversityresult} calculates the diversity statistics for the specified level of a selected environmental variable. Function \code{diversitycomp} calculates the diversity statistics for all  levels of a selected environmental variable separatedly. When a second environmental variable is provided, function \code{diversitycomp} calculates diversity statistics as a crosstabulation of both variables.
}

\value{
The functions provide alternative methods of obtaining diversity results. For function \code{diversitycomp}, the number of sites is provided as "n".
}

\references{
Kindt, R. & Coe, R. (2005)
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{http://www.worldagroforestry.org/output/tree-diversity-analysis} 
}   

\author{Roeland Kindt (World Agroforestry Centre)}

\examples{

library(vegan)
data(dune.env)
data(dune)

diversityresult(dune, y=NULL, index="Shannon", method="each site", 
    sortit=TRUE, digits=5)
diversityresult(dune, y=dune.env, factor="Management", level="NM", 
    index="Shannon", method="each site", 
    sortit=TRUE, digits=5)
diversityresult(dune, y=NULL, index="Shannon", method="pooled", digits=5)
diversityresult(dune, y=dune.env, factor="Management", level="NM", 
    index="Shannon", method="pooled", digits=5)
diversityresult(dune, y=NULL, index="Shannon", method="mean", 
    digits=5)
diversityresult(dune, y=NULL, index="Shannon", method="sd", 
    digits=5)
diversityresult(dune, y=NULL, index="Shannon", method="jackknife", 
    digits=5)
diversityresult(dune, y=dune.env, factor="Management", level="NM", 
    index="Shannon", method="jackknife", digits=5)

diversitycomp(dune, y=dune.env, factor1="Moisture", index="Shannon",
    method="pooled", sortit=TRUE)
diversitycomp(dune, y=dune.env, factor1="Moisture", index="Shannon",
    method="mean", sortit=TRUE)
diversitycomp(dune, y=dune.env, factor1="Management", index="Shannon",
    method="jackknife", sortit=TRUE)

diversitycomp(dune, y=dune.env, factor1="Management", factor2="Moisture", 
    index="Shannon", method="pooled", digits=6)
diversitycomp(dune, y=dune.env, factor1="Management", factor2="Moisture", 
    index="Shannon", method="mean", digits=6)
}

\keyword{multivariate} 


