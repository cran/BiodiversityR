\name{accumresult}
\alias{accumresult}
\alias{accumplot}
\alias{accumcomp}

\title{Alternative Species Accumulation Curve Results}

\description{
Provides alternative methods of obtaining species accumulation results than provided by functions \code{\link[vegan]{specaccum}} and \code{\link[vegan:specaccum]{plot.specaccum}} (\pkg{vegan}).
}

\usage{
accumresult(x, y="", factor="", level, scale="", method="exact", permutations=100,
    conditioned=T, gamma="boot", ...)

accumplot(xr, addit=F, labels="", col=1, ci=2, pch=1, type="p", cex=1, 
    xlim=c(1, xmax), ylim=c(1, rich),
    xlab="sites", ylab="species richness", cex.lab=1, cex.axis=1, ...)

accumcomp(x, y="", factor, scale="", method="exact", permutations=100,
    conditioned=T, gamma="boot", plotit=T, labelit=T, legend=T, rainbow=T,
    xlim=c(1, max), ylim=c(0, rich),type="p", xlab="sites",
    ylab="species richness", cex.lab=1, cex.axis=1, ...)
}

\arguments{
  \item{x}{ Community data frame with sites as rows, species as columns and species abundance as cell values. }
  \item{y}{ Environmental data frame. }
  \item{factor}{ Variable of the environmental data frame that defines subsets to calculate species accumulation curves for. }
  \item{level}{ Level of the variable to create the subset to calculate species accumulation curves. }
  \item{scale}{ Continuous variable of the environmental data frame that defines the variable that scales the horizontal axis of the species accumulation curves. }
  \item{method}{ Method of calculating the species accumulation curve (as in function \code{\link[vegan]{specaccum}}). Method "collector" adds sites in the order they happen to be in the data, "random" adds sites in random order, "exact" finds the expected (mean) species richness, "coleman" finds the expected richness following Coleman et al. 1982, and "rarefaction" finds the mean when accumulating individuals instead of sites.  }
  \item{permutations}{ Number of permutations to calculate the species accumulation curve (as in function \code{\link[vegan]{specaccum}}). }
  \item{conditioned}{ Estimation of standard deviation is conditional on the empirical dataset for the exact SAC (as in function \code{\link[vegan]{specaccum}}). }
  \item{gamma}{ Method for estimating the total extrapolated number of species in the survey area (as in \code{\link[vegan]{specaccum}}). }

  \item{addit}{ Add species accumulation curve to an existing graph. }
  \item{xr}{ Result from \code{\link[vegan]{specaccum}} or \code{accumresult}. }
  \item{col}{ Colour for drawing lines of the species accumulation curve (as in function \code{\link[vegan:specaccum]{plot.specaccum}}). }
  \item{labels}{ Labels to plot at left and right of the species accumulation curves. }
  \item{ci}{ Multiplier used to get confidence intervals from standard deviatione (as in function \code{\link[vegan:specaccum]{plot.specaccum}}). }
  \item{pch}{ Symbol used for drawing the species accumulation curve (as in function \code{\link[graphics]{points}}). }
  \item{type}{ Type of plot (as in function \code{\link[graphics]{plot}}). }
  \item{cex}{ Character expansion factor (as in function \code{\link[graphics]{plot}}). }
  \item{xlim}{ Limits for the X = horizontal axis. }
  \item{ylim}{ Limits for the Y = vertical axis. }
  \item{xlab}{ Label for the X = horizontal axis (as in function \code{\link[graphics]{title}}). }
  \item{ylab}{ Label for the Y = vertical axis (as in function \code{\link[graphics]{title}}). }
  \item{cex.lab}{ The magnification to be used for X and Y labels relative to the current setting of \code{cex}. (as in function \code{\link[graphics]{par}}). }
  \item{cex.axis}{ The magnification to be used for axis annotation relative to the current setting of \code{cex} (as in function \code{\link[graphics]{par}}). }

  \item{plotit}{ Plot the results. }
  \item{labelit}{ Label the species accumulation curves with the levels of the categorical variable. }
  \item{legend}{ Add the legend (you need to click in the graph where the legend needs to be plotted). }
  \item{rainbow}{ Use rainbow colouring for the different curves. }
  \item{...}{ Other items passed to function \code{\link[vegan]{specaccum}} or \code{\link[vegan:specaccum]{plot.specaccum}}. }
}

\details{
These functions provide some alternative methods of obtaining species accumulation results, although function \code{\link[vegan]{specaccum}} is called by these functions to calculate the actual species accumulation curve.

Functions \code{accumresult} and \code{accumcomp} allow to calculate species accumulation curves for subsets of the community and environmental data sets. Function \code{accumresult} calculates the species accumulation curve for the specified level of a selected environmental variable. Method \code{accumcomp} calculates the species accumulation curve for all levels of a selected environmental variable separatedly. Both methods allow to scale the horizontal axis by multiples of the average of a selected continuous variable from the environmental dataset (hint: add the abundance of each site to the environmental data frame to scale accumulation results by mean abundance).

Functions \code{accumcomp} and \code{accumplot} provide alternative methods of plotting species accumulation curve results,  although function \code{\link[vegan:specaccum]{plot.specaccum}} is called by these functions. When you choose to add a legend, make sure that you click in the graph on the spot where you want to put the legend. 

}

\value{
The functions provide alternative methods of obtaining species accumulation curve results, although results are similar as obtained by functions \code{\link[vegan]{specaccum}} and \code{\link[vegan:specaccum]{plot.specaccum}}.
}

\references{
Kindt, R. & Coe, R. (2005)
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{http://www.worldagroforestry.org/output/tree-diversity-analysis} }   

\author{Roeland Kindt (World Agroforestry Centre)}

\examples{
library(vegan)
data(dune.env)
data(dune)
dune.env$site.totals <- apply(dune,1,sum)
Accum.1 <- accumresult(dune, y=dune.env, scale='site.totals', method='exact', conditioned=TRUE)
Accum.1
accumplot(Accum.1)
accumcomp(dune, y=dune.env, factor='Management', method='exact', legend=FALSE, conditioned=TRUE)
## CLICK IN THE GRAPH TO INDICATE WHERE THE LEGEND NEEDS TO BE PLACED FOR
## OPTION WHERE LEGEND=TRUE (DEFAULT).
}

\keyword{multivariate} 


