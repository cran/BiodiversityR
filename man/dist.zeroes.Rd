\name{dist.zeroes}
\alias{dist.zeroes}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Distance Matrix Transformation}

\description{
Sample units without any species result in "NaN" values in the distance matrix for some of the methods of \code{\link[vegan]{vegdist}} (\pkg{vegan}). The function replaces "NA" by "0" if both sample units do not contain any species and "NA" by "1" if only one sample unit does not have any species.  
}

\usage{
dist.zeroes(comm,dist)
}

\arguments{
  \item{comm}{ Community data frame with sites as rows, species as columns and species abundance as cell values. }
  \item{dist}{ Distance matrix as calculated with function \code{\link[vegan]{vegdist}}. }
}

\details{
This functions changes a distance matrix by replacing "NaN" values by "0" if both sample units do not contain any species and by "1" if only one sample unit does not contain any species. 

Please note that there is a valid reason (deliberate removal of zero abundance values from calculations) that the original distance matrix contains "NaN", so you may not wish to do this transformation and remove sample units with zero abundances from further analysis.

}

\value{
The function provides a new distance matrix where "NaN" values have been replaced by "0" or "1".
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
matrix <- array(0,dim=c(5,3))
matrix[4,] <- c(1,2,3)
matrix[5,] <- c(1,0,0)
dist1 <- vegdist(matrix,method="kulc")
dist1
dist2 <- dist.zeroes(matrix,dist1)
dist2
}

\keyword{multivariate} 


