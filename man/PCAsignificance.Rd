\name{PCAsignificance}
\alias{PCAsignificance}
\alias{ordiequilibriumcircle}

\title{PCA Significance}

\description{
Calculates the number of significant axes from a Principal Components Analysis based on the broken-stick criterion, or adds an equilibrium circle to an ordination diagram.
}

\usage{
PCAsignificance(pca,axes=8)
ordiequilibriumcircle(pca,ordiplot,...)
}

\arguments{
  \item{pca}{ Principal Components Analysis result as calculated by \code{\link[vegan:cca]{rda}} (\pkg{vegan}). }
  \item{axes}{ Number of axes to calculate results for. }
  \item{ordiplot}{ Ordination plot created by \code{\link[vegan]{ordiplot}} (\pkg{vegan}) }
  \item{...}{ Other arguments passed to function \code{\link[graphics]{arrows}}. }
}

\details{
These functions provide two methods of providing some information on significance for a Principal Components Analysis (PCA).

Function \code{PCAsignificance} uses the broken-stick distribution to evaluate how many PCA axes are significant. This criterion is one of the most reliable to check how many axes are significant. PCA axes with larger percentages of (accumulated) variance than the broken-stick variances are significant (Legendre and Legendre, 1998).

Function \code{ordiequilibriumcircle} draws an equilibirum circle to a PCA ordination diagram. Only species vectors with heads outside of the equilibrium circle significantly contribute to the ordination diagram (Legendre and Legendre, 1998). Vectors are drawn for these species. The function considers the scaling methods used by \code{\link[vegan:cca]{rda}} for \code{scaling=1}. The method should only be used for \code{scaling=1} and PCA calculated by function \code{\link[vegan:cca]{rda}}. 
}

\value{
Function \code{PCAsignificance} returns a matrix with the variances that are explained by the PCA axes and by the broken-stick criterion.
 
Function \code{ordiequilibriumcircle} plots an equilibirum circle and returns a list with the radius and the scaling constant used by \code{\link[vegan:cca]{rda}}.
}

\references{
Legendre, P. & Legendre, L. (1998). Numerical Ecology. 2nd English Edition. Elsevier. 

Kindt, R. & Coe, R. (2005). 
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{http://www.worldagroforestry.org/treesandmarkets/tree_diversity_analysis.asp} }   

\author{Roeland Kindt (World Agroforestry Centre)}

\examples{
library(vegan)
data(dune)
Ordination.model1 <- rda(dune)
PCAsignificance(Ordination.model1)
plot1 <- ordiplot(Ordination.model1, choices=c(1,2), scaling=1)
ordiequilibriumcircle(Ordination.model1,plot1)
}

\keyword{multivariate} 


