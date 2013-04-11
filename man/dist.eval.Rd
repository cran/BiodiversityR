\name{dist.eval}
\alias{dist.eval}

\title{Distance Matrix Evaluation}

\description{
Provides one test of a distance matrix, and continuous with \code{\link[vegan]{distconnected}} (\pkg{vegan}).
}

\usage{
dist.eval(x,dist)
}

\arguments{
  \item{x}{ Community data frame with sites as rows, species as columns and species abundance as cell values. }
  \item{dist}{ Method for calculating ecological distance with function \code{\link[vegan]{vegdist}}: partial match to "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "morisita", "horn" or "mountford". }
}

\details{
This functions provides two tests of a distance matrix. 

The first test checks whether any pair of sites that share some species have a larger distance than any other pair of sites that do not share any species. In case that cases are found, then a warning message is given. 

The second test is the one implemented by the \code{\link[vegan]{distconnected}} function (\pkg{vegan}). The  distconnected test is only calculated for distances that calculate a value of 1 if sites share no species (i.e. not manhattan or euclidean), using the threshold of 1 as an indication that the sites do not share any species. Interpretation of analysis of distance matrices that provided warnings should be cautious.
}

\value{
The function tests whether distance matrices have some desirable properties and provide warnings if this is not the case.
}

\references{
Kindt, R. & Coe, R. (2005)
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{http://www.worldagroforestry.org/resources/databases/tree-diversity-analysis}
 }   

\author{Roeland Kindt (World Agroforestry Centre)}

\examples{
library(vegan)
data(dune)
dist.eval(dune,"euclidean")
dist.eval(dune,"bray")
}

\keyword{multivariate} 


