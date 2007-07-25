\name{makecommunitydataset}
\alias{makecommunitydataset}

\title{Make a Community Dataset from a Stacked Dataset}

\description{
Makes a community data set from a stacked dataset (with separate variables for the site identities, the species identities and the abundance).
}

\usage{
makecommunitydataset(x,row,column,value,factor="",level="",drop=F)
}

\arguments{
  \item{x}{ Data frame. }
  \item{row}{ Name of the categorical variable for the rows of the crosstabulation (typically indicating sites) }
  \item{column}{ Name of the categorical variable for the columns of the crosstabulation (typically indicating species) }
  \item{value}{ Name of numerical variable for the cells of the crosstabulation (typically indicating abundance). The cells provide the sum of all values in the data frame. }
  \item{factor}{ Name of the variable to calculate a subset of the data frame. }
  \item{level}{ Value of the subset of the factor variable to calculate a subset of the data frame. }
  \item{drop}{ Drop rows without species (species with total abundance of zero are always dropped) }
}

\details{
This function calculates a cross-tabulation from a data frame, summing up all the values of the numerical variable identified as variable for the cell values. If factor="", then no subset is calculated from the data frame in the first step.
}

\value{
The function provides a community dataset from another data frame.
}

\references{
Kindt, R. & Coe, R. (2005) 
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{http://www.worldagroforestry.org/treesandmarkets/tree_diversity_analysis.asp}
}   

\author{Roeland Kindt (World Agroforestry Centre)}

\examples{
library(vegan)
data(dune.env)
makecommunitydataset(dune.env,row='Manure',column='Use',value='A1')
makecommunitydataset(dune.env,row='Manure',column='Use',value='A1',
    factor='Management',level='NM')
}

\keyword{multivariate} 


