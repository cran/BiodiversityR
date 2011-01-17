\name{balanced.specaccum}
\alias{balanced.specaccum}

\title{Balanced Species Accumulation Curves}

\description{
Provides species accumulation results calculated from balanced (equal subsample sizes) subsampling from each stratum. Sites can be accumulated in a randomized way, or alternatively sites belonging to the same stratum can be kept together Results are in the same format as \code{\link[vegan]{specaccum}} and can be plotted with \code{\link[vegan:specaccum]{plot.specaccum}} (\pkg{vegan}).
}

\usage{
balanced.specaccum(comm, permutations=100, strata=strata, grouped=TRUE, 
    reps=0, scale=NULL)
}

\arguments{
  \item{comm}{ Community data frame with sites as rows, species as columns and species abundance as cell values. }
  \item{permutations}{ Number of permutations to calculate the species accumulation curve. }
  \item{strata}{ Categorical variable used to specify strata. }
  \item{grouped}{ Should sites from the same stratum be kept together (TRUE) or not. }
  \item{reps}{ Number of subsamples to be taken from each stratum (see details). }
  \item{scale}{ Quantitative variable used to scale the sampling effort (see details). }
}

\details{
This function provides an alternative method of obtaining species accumulation results as provided by \code{\link[vegan]{specaccum}} and \code{\link{accumresult}}.
 
Balanced sampling is achieved by randomly selecting the same number of sites from each stratum. The number of sites selected from each stratum is determined by \code{reps}. Sites are selected from strata with sample sizes larger or equal than \code{reps}. In case that \code{reps} is smaller than 1 (default: 0), then the number of sites selected from each stratum is equal to the smallest sample size of all strata. Sites from the same stratum can be kept together (grouped=TRUE) or the order of sites can be randomized (grouped=FALSE).
 
The results can be scaled by the average accumulation of a quantitative variable (default is number of sites), as in \code{accumresult} (hint: add the abundance of each site to the environmental data frame to scale accumulation results by mean abundance). When sites are not selected from all strata, then the average is calculated only for the strata that provided sites.
 
}

\value{
The functions provide alternative methods of obtaining species accumulation curve results, although results are similar as obtained by functions \code{\link[vegan]{specaccum}} and \code{\link{accumresult}}.
}

\references{
Kindt, R. & Coe, R. (2005)
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{http://www.worldagroforestry.org/resources/databases/tree-diversity-analysis} }   

\author{Roeland Kindt (World Agroforestry Centre)}

\examples{
library(vegan)
data(dune.env)
data(dune)
Accum.1 <- balanced.specaccum(dune, strata=dune.env$Management, reps=3)
Accum.1
dune.env$site.totals <- apply(dune,1,sum)
Accum.2 <- balanced.specaccum(dune, strata=dune.env$Management, reps=3, scale=dune.env$site.totals)
Accum.2
}

\keyword{multivariate} 


