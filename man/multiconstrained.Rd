\name{multiconstrained}
\alias{multiconstrained}

\title{Pairwise Comparisons for All Levels of a Categorical Variable by RDA, CCA or Capscale}

\description{
This function implements pairwise comparisons for categorical variable through \code{\link[vegan]{capscale}}, \code{\link[vegan]{cca}} or \code{\link[vegan:cca]{rda}} followed by \code{\link[vegan]{anova.cca}}. The function simply repeats constrained ordination analysis by selecting subsets of data that correspond to two factor levels. 
}

\usage{
multiconstrained(method="capscale", formula, data, distance = "bray"
    , comm = NULL, add = FALSE, multicomp="", contrast=0, ...)
}

\arguments{
  \item{method}{ Method for constrained ordination analysis; one of "rda", "cca"or "capscale". }
  \item{formula}{ Model formula as in \code{\link[vegan]{capscale}},  \code{\link[vegan]{cca}} or  \code{\link[vegan:cca]{rda}}. The LHS can be a community data matrix or a distance matrix for \code{\link[vegan]{capscale}}. }
  \item{data}{ Data frame containing the variables on the right hand side of the model formula as in \code{\link[vegan]{capscale}}, \code{\link[vegan]{cca}} or  \code{\link[vegan:cca]{rda}}. }
  \item{distance}{ Dissimilarity (or distance) index in vegdist used if the LHS of the formula is a data frame instead of dissimilarity matrix; used only with function \code{\link[vegan]{vegdist}} and partial match to "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "morisita", "horn" or "mountford". This argument is only used for \code{\link[vegan]{capscale}} in case that the LHS of the formula is a community matrix.  }
  \item{comm}{ Community data frame which will be used for finding species scores when the LHS of the formula was a dissimilarity matrix as only allowed for \code{\link[vegan]{capscale}}. This is not used if the LHS is a data frame.  }
  \item{add}{ Logical indicating if an additive constant should be computed, and added to the non-diagonal dissimilarities such that all eigenvalues are non-negative in underlying Principal Co-ordinates Analysis; only applicable in \code{\link[vegan]{capscale}}.  }
  \item{multicomp}{ Categorical variable used to construct the contrasts from. In case that this variable is missing, then the first explanatory variable of the formula will be used.  }
  \item{contrast}{ Return the ordination results for the particular contrast indicated by this number (e.g. with 5 levels, one can choose in between contrast 1-10). In case=0, then the first row of the \code{\link[vegan]{anova.cca}} results for all contrasts is provided.  }
  \item{...}{ Other parameters passed to \code{\link[vegan]{anova.cca}}.  }
}

\details{
This function provides a simple expansion of \code{\link[vegan]{capscale}}, \code{\link[vegan]{cca}} and \code{\link[vegan:cca]{rda}} by conducting the analysis for subsets of the community and environmental datasets that only contain two levels of a categoricl variable. 
 
When the choice is made to return results from all contrasts (contrast=0), then the first row of the \code{\link[vegan]{anova.cca}} tables for each contrast are provided. It is therefore possible to compare differences in results by modifying the "by" argument of this function (i.e. obtain the total of explained variance, the variance explained on the first axis or the variance explained by the variable alone). 
 
When the choice is made to return results from a particular contrast (contrast>0), then the ordination result is returned and two new datasets ("newcommunity" and "newenvdata") are created that only contain data for the two selected contrasts.
  
}

\value{
The function returns an ANOVA table that contains the first rows of the ANOVA tables obtained for all possible combinations of levels of the first variable. Alternatively, it returns an ordination result for the selected contrast and creates two new datasets ("newcommunity" and "newenvdata")
}

\references{
Legendre, P. & Anderson, M.J. (1999). Distance-based redundancy analysis:  
testing multispecies responses in multifactorial ecological experiments. 
Ecological Monographs 69: 1-24.
 
Anderson, M.J. & Willis, T.J. (2003). Canonical analysis of principal coordinates: a useful method of constrained ordination for ecology. 
Ecology 84: 511-525.
 }   

\author{Roeland Kindt (World Agroforestry Centre)}

\examples{
\dontrun{
library(vegan)
library(MASS)
data(dune)
data(dune.env)
multiconstrained(method="capscale", dune~Management, data=dune.env,
    distance="bray",add=TRUE)
multiconstrained(method="capscale", dune~Management+Condition(A1),
    data=dune.env, distance="bray", add=TRUE, contrast=3)
}
}

\keyword{multivariate} 


