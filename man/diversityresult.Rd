\name{diversityresult}
\alias{diversityresult}
\alias{diversitycomp}

\title{Alternative Diversity Results}

\description{
Provides alternative methods of obtaining results on diversity statistics than provided directly by functions \code{\link[vegan]{diversity}}, \code{\link[vegan:diversity]{fisher.alpha}}, \code{\link[vegan]{specpool}}, \code{\link[vegan:specpool]{estimateR}} and \code{\link[vegan:diversity]{specnumber}} (all from \pkg{vegan}), although these same functions are called. Some other statistics are also calculated such as the reciprocal Berger-Parker diversity index and abundance (not a diversity statistic). The statistics can be calculated for the entire community, for each site separately, the mean of the sites can be calculated or a jackknife estimate can be calculated for the community.
}

\usage{
diversityresult(x,y="",factor,level,index="Shannon",method="all",sortit=F,
    digits=8)
diversitycomp(x,y="",factor1,factor2="",index="Shannon",
    method="all",sortit=F,...)
}

\arguments{
  \item{x}{ Community data frame with sites as rows, species as columns and species abundance as cell values. }
  \item{y}{ Environmental data frame. }
  \item{factor}{ Variable of the environmental data frame that defines subsets to calculate diversity statistics for. }
  \item{level}{ Level of the variable to create the subset to calculate diversity statistics. }
  \item{index}{ Type of diversity statistic with "richness" to calculate species richness, "abundance" to calculate abundance, "Shannon" to calculate the Shannon diversity index, "Simpson" to calculate 1-Simpson concentration index, "inverseSimpson" to calculate the reciprocal Simpson diversity index, "Logalpha" to calculate the log series alpha diversity index, "Berger" to calculate the reciprocal Berger-Parker diversity index, "Jevenness" to calculate one Shannon evenness index, "Eevenness" to calculate another Shannon evenness index, "Jack.1" to calculate the first-order jackknife gamma diversity estimator, "Jack.2" to calculate the first-order jackknife gamma diversity estimator, "Jack.2" to calculate the first-order jackknife gamma diversity estimator, "Chao" to calculate the Chao gamma diversity estimator, "Boot" to calculate the bootstrap gamma diversity estimator and "ACE" to calculate an abundance-based gamma diversity estimator. }
  \item{method}{ Method of calculating the diversity statistics: "all" calculates the diversity of the entire community (all sites pooled), "s" calculates the diversity of each site separatedly, "mean" calculates the average diversity of the sites, "sd" calculates the standard deviation of the diversity of the sites, whereas "Jackknife" calculates the jackknifed diversity for the entire data frame. Method "s" is not allowed for \code{diversitycomp}.  }
  \item{sortit}{ Sort the sites by increasing values of the diversity statistic. }
  \item{digits}{ Number of digits in the results. }
  \item{factor1}{ Variable of the environmental data frame that defines subsets to calculate diversity statistics for. }
  \item{factor2}{ Optional second variable of the environmental data frame that defines subsets to calculate diversity statistics for in a crosstabulation with the other variable of the environmental data frame. }
  \item{...}{ Other arguments passed to function \code{diversityresult}. }
}

\details{
These functions provide some alternative methods of obtaining results with diversity statistics, although functions \code{\link[vegan]{diversity}}, \code{\link[vegan:diversity]{fisher.alpha}}, \code{\link[vegan]{specpool}}, \code{\link[vegan:specpool]{estimateR}} and \code{\link[vegan:diversity]{specnumber}} (all from \pkg{vegan}) are called to calculate the various statistics.

The reciprocal Berger-Parker diversity index is the reciprocal of the proportional abundance of the most dominant species.

J-evenness is calculated as: H / ln(S) where H is the Shannon diversity index and S the species richness.

E-evenness is calculated as: exp(H) / S where H is the Shannon diversity index and S the species richness.

The method of calculating the diversity statistics include following options: "all" calculates the diversity of the entire community (all sites pooled together), "s" calculates the diversity of each site separatedly, "mean" calculates the average diversity of the sites, whereas "Jackknife" calculates the jackknifed diversity for the entire data frame. Methods "s" and "mean" are not available for function diversitycomp.

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
 
  \url{http://www.worldagroforestry.org/treesandmarkets/tree_diversity_analysis.asp} }   

\author{Roeland Kindt (World Agroforestry Centre)}

\examples{
library(vegan)
data(dune.env)
data(dune)
Diversity.1 <- diversityresult(dune, y=dune.env, factor='Management', 
    level='NM', index='Shannon' ,method='s', sortit=TRUE, digits=3)
Diversity.1
diversitycomp(dune, y=dune.env, factor1='Management', factor2="Moisture", 
    index='Shannon' ,method='all', sortit=TRUE, digits=3)
}

\keyword{multivariate} 


