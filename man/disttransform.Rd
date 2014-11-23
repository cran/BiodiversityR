\name{disttransform}
\alias{disttransform}

\title{Community Matrix Transformation}

\description{
Transforms a community matrix. Some transformation methods are described by distances for the original community matrix that result in the same distance matrix as calculated with the euclidean distance from the transformed community matrix.
In several cases (methods of "hellinger", "chord", "profiles" and "chi.square), the method makes use of function \code{\link[vegan]{decostand}}. In several other cases ("Braun.Blanquet", "Domin", "Hult", "Hill", "fix" and "coverscale.log"), the method makes use of function \code{\link[vegan]{coverscale}}.
}

\usage{
disttransform(x, method="hellinger")
}

\arguments{
  \item{x}{ Community data frame with sites as rows, species as columns and species abundance as cell values. }
  \item{method}{ Distance measure for the original community matrix that the euclidean distance will calculate for the transformed community matrix: partial match to "hellinger", "chord", "profiles", "chi.square", "log", "square", "pa", "Braun.Blanquet", "Domin", "Hult", "Hill", "fix" and "coverscale.log". }
}

\details{
This functions transforms a community matrix.

Some transformation methods ("hellinger", "chord", "profiles" and "chi.square") have the behaviour that the euclidean distance from the transformed matrix will equal a distance of choice for the original matrix. For example, using method "hellinger" and calculating the euclidean distance will result in the same distance matrix as by calculating the Hellinger distance from the original community matrix.

Transformation methods ("Braun.Blanquet", "Domin", "Hult", "Hill", "fix" and "coverscale.log") call function \code{\link[vegan]{coverscale}}.
}

\value{
The function returns a transformed community matrix.
}

\references{
Legendre, P. & Gallagher, E.D. (2001). Ecologically meaningful transformations 
  for ordination of species data. Oecologia 129: 271-280.

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
Community.1 <- disttransform(dune, method='hellinger')
Distmatrix.1 <- vegdist(Community.1,method='euclidean')
Distmatrix.1
}

\keyword{multivariate} 


