\name{add.spec.scores}
\alias{add.spec.scores}

\title{Add Species Scores to Unconstrained Ordination Results}

\description{
Calculates scores (coordinates) to plot species for PCoA or NMS results that do not naturally provide species scores. The function can also rescale PCA results to use the choice of rescaling used in \pkg{vegan} for the \code{\link[vegan:cca]{rda}} function (after calculating PCA results via PCoA with the euclidean distance first).
}

\usage{
add.spec.scores(ordi,comm,method="cor.scores",multi=1,Rscale=F,scaling="1")
}

\arguments{
  \item{ordi}{ Ordination result as calculated by \code{\link[stats]{cmdscale}},  \code{\link[MASS]{isoMDS}},  \code{\link[MASS]{sammon}},  \code{\link[vegan:metaMDS]{postMDS}},  \code{\link[vegan]{metaMDS}} or  \code{\link{NMSrandom}}. }
  \item{comm}{ Community data frame with sites as rows, species as columns and species abundance as cell values. }
  \item{method}{ Method for calculating species scores. Method "cor.scores" calculates the scores by the correlation between site scores and species vectors (via function \code{\link[stats]{cor}}), method "wa.scores" calculates the weighted average scores (via function \code{\link[vegan]{wascores}}) and method "pcoa.scores" calculates the scores by weighing the correlation between site scores and species vectors by variance explained by the ordination axes. }
  \item{multi}{ Multiplier for the species scores. }
  \item{Rscale}{ Use the same scaling method used by \pkg{vegan} for \code{\link[vegan:cca]{rda}}. }
  \item{scaling}{ Scaling method as used by \code{\link[vegan:cca]{rda}}. }
}

\value{
The function returns a new ordination result with new information on species scores. For PCoA results, the function calculates eigenvalues (not sums-of-squares as provided in results from function \code{\link[stats]{cmdscale}}), the percentage of explained variance per axis and the sum of all eigenvalues. PCA results (obtained by PCoA obtained by function \code{\link[stats]{cmdscale}} with the Euclidean distance) can be scaled as in function \code{\link[vegan:cca]{rda}}, or be left at the original scale.
}

\references{
Kindt, R. & Coe, R. (2005)
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{http://www.worldagroforestry.org/resources/databases/tree-diversity-analysis} }   

\author{Roeland Kindt}

\examples{
library(vegan)
data(dune)
distmatrix <- vegdist(dune, method='euc')
## Principal coordinates analysis with 19 axes to estimate total variance
Ordination.model1 <- cmdscale(distmatrix, k=19, eig=TRUE, add=FALSE)
Ordination.model1 <- add.spec.scores(Ordination.model1,dune, 
    method='pcoa.scores', Rscale=TRUE, scaling=1, multi=1)
Ordination.model1
## Compare Ordination.model1 with:
Ordination.model2 <- rda(dune)
summary(Ordination.model2, scaling=1)
}

\keyword{multivariate} 


