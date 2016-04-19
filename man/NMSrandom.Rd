\name{NMSrandom}
\alias{NMSrandom}

\title{Calculate the NMS Result with the Smallest Stress from Various Random Starts }

\description{
This function provides a simplified version of the method of calculating NMS results implemented by the function \code{\link[vegan]{metaMDS}} (\pkg{vegan}).
}

\usage{
NMSrandom(x,perm=100,k=2,stressresult=F,method="isoMDS")
}

\arguments{
  \item{x}{ Distance matrix. }
  \item{perm}{ Number of permutations to select the configuration with the lowest stress. }
  \item{k}{ Number of dimensions for the non metric scaling result; passed to \code{\link[MASS]{isoMDS}} or \code{\link[MASS]{sammon}}. }
  \item{stressresult}{ Provide the calculated stress for each permutation. }
  \item{method}{ Method for calculating the NMS:  \code{\link[MASS]{isoMDS}} or \code{\link[MASS]{sammon}}. }
}

\details{
This function is an easier method of calculating the best NMS configuration after various random starts than implemented in the \code{\link[vegan]{metaMDS}} function (\pkg{vegan}). The function uses a distance matrix (as calculated for example by function \code{\link[vegan]{vegdist}} from a community data set) and calculates random starting positions by function \code{\link[vegan:metaMDS]{initMDS}} (\pkg{vegan}) analogous to \code{\link[vegan]{metaMDS}}.
}

\value{
The function returns the NMS ordination result with the lowest stress (calculated by \code{\link[MASS]{isoMDS}} or \code{\link[MASS]{sammon}}.), or the stress of each NMS ordination.
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
library(MASS)
data(dune)
distmatrix <- vegdist(dune)
Ordination.model1 <- NMSrandom(distmatrix,perm=100,k=2)
Ordination.model1 <- add.spec.scores(Ordination.model1,dune, 
    method='wa.scores')
Ordination.model1
}

\keyword{multivariate} 


