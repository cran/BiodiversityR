\name{nested.anova.dbrda}
\alias{nested.anova.dbrda}
\alias{nested.npmanova}

\title{Nested Analysis of Variance via Distance-based Redundancy Analysis or Non-parametric Multivariate Analysis of Variance}

\description{
The functions provide nested analysis of variance for a two-level hierarchical model. The functions are implemented by estimating the correct F-ratio for the main and nested factors (assuming the nested factor is random) and using the recommended permutation procedures to test the significance of these F-ratios. F-ratios are estimated from variance estimates that are provided by distance-based redundancy analysis (\code{\link[vegan]{capscale}}) or non-parametric multivariate analysis of variance (\code{\link[vegan]{adonis}}).
}

\usage{
nested.anova.dbrda(formula, data, method="euc", add=FALSE,
    permutations=100, warnings=FALSE)
nested.npmanova(formula, data, method="euc", permutations=100, warnings=FALSE)
}

\arguments{

  \item{formula}{ Formula with a community data frame (with sites as rows, species as columns and species abundance as cell values) or (for \code{nested.anova.dbrda} only) distance matrix on the left-hand side and two categorical variables on the right-hand side (with the second variable assumed to be nested within the first). }
  \item{data}{ Environmental data set. }
  \item{method}{ Method for calculating ecological distance with function \code{\link[vegan]{vegdist}}: partial match to "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "morisita", "horn" or "mountford". This argument is ignored in case that the left-hand side of the formula already is a distance matrix. }
  \item{add}{ Should a constant be added to the off-diagonal elements of the distance-matrix (TRUE) or not. }
  \item{permutations}{ The number of permutations for significance testing.  }
  \item{warnings}{ Should warnings be suppressed (TRUE) or not. }
}

\details{
The functions provide two alternative procedures for multivariate analysis of variance on the basis of any distance measure. Function \code{nested.anova.dbrda} proceeds via \code{\link[vegan]{capscale}}, whereas \code{nested.npmanova} proceeds via \code{\link[vegan]{adonis}}. Both methods are complementary to each other as \code{nested.npmanova} always provides correct F-ratios and estimations of significance, whereas \code{nested.anova.dbrda} does not provide correct F-ratios and estimations of significance when negative eigenvalues are encountered or constants are added to the distance matrix, but always provides an ordination diagram. 
 
The F-ratio for the main factor is estimated as the mean square of the main factor divided by the mean square of the nested factor. The significance of the F-ratio of the main factor is tested by permuting entire blocks belonging to levels of the nested factor. The significance of the F-ratio of the nested factor is tested by permuting sample units within strata defined by levels of the main factor.

}


\value{
The functions provide an ANOVA table.
}

\references{
Legendre, P. & Anderson, M. J. (1999). Distance-based redundancy analysis: testing multispecies responses in multifactorial ecological experiments. Ecological Monographs 69, 1-24.

Anderson, M.J. (2001). A new method for non-parametric multivariate analysis of variance. Austral Ecology, 26: 32-46.

McArdle, B.H. and M.J. Anderson. (2001). Fitting multivariate models to community data: A comment on distance-based redundancy analysis. Ecology, 82: 290-297.
}

\author{Roeland Kindt (World Agroforestry Centre)}

\examples{
library(vegan)
data(warcom)
data(warenv)
# use larger number of permutations for real studies
nested.npmanova(warcom~rift.valley+popshort, data=warenv, method="jac", 
    permutations=10)
nested.anova.dbrda(warcom~rift.valley+popshort, data=warenv, method="jac", 
    permutations=10)
}


\keyword{multivariate} 


