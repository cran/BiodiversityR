\name{CAPdiscrim}
\alias{CAPdiscrim}

\title{Canonical Analysis of Principal Coordinates based on Discriminant Analysis}

\description{
This function provides a method for CAP as described by the authors of the ordination method. The CAP method implemented in \pkg{vegan} through \code{\link[vegan]{capscale}} conforms more to distance-based Redundancy Analysis (Legendre & Anderson, 1999) than to the original description for CAP (Anderson & Willis, 2003 ).
}

\usage{
CAPdiscrim(formula,data,dist="bray",axes=4,m=0,permutations=0)
}

\arguments{
  \item{formula}{ Formula with a community data frame (with sites as rows, species as columns and species abundance as cell values) or distance matrix on the left-hand side and a categorical variable on the right-hand side (only the first explanatory variable will be used). }
  \item{data}{ Environmental data set. }
  \item{dist}{ Method for calculating ecological distance with function \code{\link[vegan]{vegdist}}: partial match to "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "morisita", "horn" or "mountford". This argument is ignored in case that the left-hand side of the formula already is a distance matrix. }
  \item{axes}{ Number of PCoA axes (\code{\link[stats]{cmdscale}}) to provide in the result. }
  \item{m}{ Number of PCoA axes to be investigated by discriminant analysis (\code{\link[MASS]{lda}}). If m=0 then the number of axes that provides the best distinction between the groups is calculated (following the method of Anderson and Willis).  }
  \item{permutations}{ The number of permutations for significance testing.  }
}

\details{
This function provides a method of Constrained Analysis of Principal Coordinates (CAP) that conforms to the description of the method by the developers of the method, Anderson and Willis. The method investigates the results of a Principal Coordinates Analysis (function \code{\link[stats]{cmdscale}}) with linear discriminant analysis (\code{\link[MASS]{lda}}). Anderson and Willis advocate to use the number of principal coordinate axes that result in the best prediction of group identities of the sites.
  
For permutations > 0, the analysis is repeated by randomising the observations of the environmental data set. The significance is estimated by dividing the number of times the randomisation generated a larger percentage of correct predictions.

}

\value{
The function returns an object with information on CAP based on discriminant analysis. The object contains following elements:  

\item{PCoA}{ the positions of the sites as fitted by PCoA}
\item{m}{ the number of axes analysed by discriminant analysis}
\item{tot}{ the total variance (sum of all eigenvalues of PCoA)}
\item{varm}{ the variance of the m axes that were investigated}
\item{group}{ the original group of the sites}
\item{CV}{ the predicted group for the sites by discriminant analysis}
\item{percent}{ the percentage of correct predictions} 
\item{x}{ the positions of the sites provided by the discriminant analysis}
\item{F}{ the squares of the singulare values of the discriminant analysis}
\item{manova}{ the results for MANOVA with the same grouping variable}
\item{signi}{ the significance of the percentage of correct predictions}
\item{manova}{ a summary of the observed randomised prediction percentages}

The object can be plotted with \code{\link[vegan]{ordiplot}}, and species scores can be added by \code{\link{add.spec.scores}} .
}

\references{
Legendre, P. & Anderson, M.J. (1999). Distance-based redundancy analysis:  
testing multispecies responses in multifactorial ecological experiments. 
Ecological Monographs 69: 1-24.

Anderson, M.J. & Willis, T.J. (2003). Canonical analysis of principal coordinates: a useful method of constrained ordination for ecology. 
Ecology 84: 511-525.


Kindt, R. & Coe, R. (2005) 
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{http://www.worldagroforestry.org/treesandmarkets/tree_diversity_analysis.asp}
 }   

\author{Roeland Kindt (World Agroforestry Centre)}

\examples{
library(vegan)
library(MASS)
data(dune)
data(dune.env)
Ordination.model1 <- CAPdiscrim(dune~Management, data=dune.env,
    dist="bray",axes=2,m=0)
Ordination.model1
plot1 <- ordiplot(Ordination.model1)
ordisymbol(plot1,dune.env,"Management",legend=FALSE)
## CLICK IN THE GRAPH TO INDICATE THE POSITION FOR THE LEGEND
## IN CASE THAT THE OPTION WAS LEGEND=TRUE.
}

\keyword{multivariate} 


