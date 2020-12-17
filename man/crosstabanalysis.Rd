\name{crosstabanalysis}
\alias{crosstabanalysis}

\title{Presence-absence Analysis by Cross Tabulation}

\description{
This function makes a cross-tabulation of two variables after transforming the first variable to presence-absence and then returns results of \code{\link[stats]{chisq.test}}.
}

\usage{
crosstabanalysis(x,variable,factor)
}

\arguments{
  \item{x}{ Data set that contains the variables "variable" and "factor". }
  \item{variable}{ Variable to be transformed in presence-absence in the resulting cross-tabulation. }
  \item{factor}{ Variable to be used for the cross-tabulation together with the transformed variable. }
}

\value{
The function returns the results of \code{\link[stats]{chisq.test}} on a crosstabulation of two variables, after transforming the first variable to presence-absence first.
}

\references{
Kindt, R. & Coe, R. (2005)
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{https://www.worldagroforestry.org/output/tree-diversity-analysis} }   

\author{Roeland Kindt}

\examples{
library(vegan)
data(dune.env)
crosstabanalysis(dune.env,"Manure","Management")
}

\keyword{multivariate} 


