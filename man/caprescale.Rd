\name{caprescale}
\alias{caprescale}

\title{Rescaling of Capscale Results to Reflect Total Sums of Squares Of Distance Matrix}

\description{
This is a simple function that rescales the ordination coordinates obtained from the distance-based redundancy analysis method implemented in \pkg{vegan} through \code{\link[vegan]{capscale}}. The rescaling of the ordination coordinates results in the distances between fitted site scores in ordination results (scaling=1 obtained via \code{\link[vegan]{ordiplot}} to be equal to the distances between sites on the axes corresponding to positive eigenvalues obtained from principal coordinates analysis (\code{\link[stats]{cmdscale}}).
}

\usage{
caprescale(x,verbose=FALSE)
}

\arguments{
  \item{x}{ Ordination result obtained with \code{\link[vegan]{capscale}}. }
  \item{verbose}{ Give some information on the pairwise distances among sites (TRUE) or not. }
}

\details{
The first step of distance-based redundancy analysis involves principal coordinates analysis whereby the distances among sites from a distance matrix are approximated by distances among sites in a multidimensional configuration (ordination). In case that the principal coordinates analysis does not result in negative eigenvalues, then the distances from the distance matrix are the same as the distances among the sites in the ordination. In case that the principal coordinates analysis results in negative eigenvalues, then the distances among the sites on all ordination axes are related to the sum of positive eigenvalues, a sum which is larger than the sum of squared distances of the distance matrix.
 
The distance-based redundancy analysis method implemented in \pkg{vegan} through \code{\link[vegan]{capscale}} uses a specific rescaling method for ordination results. Function \code{caprescale} modifies the results of  \code{\link[vegan]{capscale}} so that an ordination with scaling=1 (a distance biplot) obtained via\code{\link[vegan]{ordiplot}} preserves the distances reflected in the principal coordinates analysis implemented as the first step of the analysis. See Legendre and Legendre (1998) about the relationship between fitted site scores and eigenvalues.

}

\value{
The function modifies and returns an object obtained via  \code{\link[vegan]{capscale}}. 
}

\references{
Legendre, P. & Legendre, L. (1998). Numerical Ecology. Amsterdam: Elsevier. 853 pp.

Legendre, P. & Anderson, M.J. (1999). Distance-based redundancy analysis:  
testing multispecies responses in multifactorial ecological experiments. 
Ecological Monographs 69: 1-24.
}

\author{Roeland Kindt (World Agroforestry Centre)}

\examples{
library(vegan)
library(MASS)
data(dune)
data(dune.env)
Distmatrix.1 <- vegdist(dune,method='bray')
Ordination.model1 <- cmdscale(Distmatrix.1, k=19, eig=TRUE, add=FALSE)
# Sum of all eigenvalues
sum(Ordination.model1$eig)
# [1] 4.395807541512926
sum(Ordination.model1$eig[1:14])
# [1] 4.593946896588808
Distmatrix.2 <- as.matrix(vegdist(Ordination.model1$points[,1:14],method='euc'))
totalsumsquares1 <- sum(Distmatrix.2^2)/(2*20)
# Sum of distances among sites in principal coordinates analysis on axes
# corresponding to positive eigenvalues
totalsumsquares1
# [1] 4.593946896588808
Ordination.model2 <- capscale(dune ~ Management,dune.env,dist='bray', add=FALSE)
# Total sums of positive eigenvalues of the distance-based redundancy analysis
Ordination.model2$CA$tot.chi+Ordination.model2$CCA$tot.chi
# [1] 4.593946896588808
Ordination.model3 <- caprescale(Ordination.model2, verbose=TRUE)
sum1 <- summary(Ordination.model3,axes=17,scaling=1)$constraints
Distmatrix.3 <- as.matrix(vegdist(sum1 ,method='euc'))
totalsumsquares2 <- sum((Distmatrix.3)^2)/(2*20)/19
totalsumsquares2
# [1] 4.593946896588808

}

\keyword{multivariate} 


