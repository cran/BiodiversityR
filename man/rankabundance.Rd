\name{rankabundance}
\alias{rankabundance}
\alias{rankabunplot}
\alias{rankabuncomp}

\title{Rank Abundance Curves}

\description{
Provides methods of calculating rank-abundance curves.
}

\usage{
rankabundance(x, y="", factor="", level, digits=1, t=qt(0.975, df=n-1))

rankabunplot(xr, addit=F, labels="", scale="abundance", scaledx=F, type="o",
    xlim=c(min(xpos), max(xpos)),
    ylim=c(0, max(x[,scale])),
    specnames=c(1:5), srt=0, ...)

rankabuncomp(x, y="", factor, scale="abundance", 
    scaledx=F, type="o", rainbow=T,
    legend=T, xlim=c(1, max1), ylim=c(0, max2), ...)
}

\arguments{
  \item{x}{ Community data frame with sites as rows, species as columns and species abundance as cell values. }
  \item{y}{ Environmental data frame. }
  \item{factor}{ Variable of the environmental data frame that defines subsets to calculate rank abundance curves for. }
  \item{level}{ Level of the variable to create the subset to calculate rank abundance curves. }
  \item{digits}{ Number of digits in the results. }
  \item{t}{ t-value to calculate confidence interval limits for the species proportion for cluster sampling (following Hayek and Buzas 1997).  }

  \item{xr}{ Result from \code{rankabundance}. }
  \item{addit}{ Add rank abundance curve to an existing graph. }
  \item{labels}{ Labels to plot at left of the rank abundance curves. }
  \item{scale}{ Method of scaling the vertical axis. Method "abundance" uses abundance, "proportion" uses proportional abundance (species abundance / total abundance), "logabun" calculates the logarithm of abundance using base 10 and "accumfreq" accumulates the proportional abundance. }
  \item{scaledx}{ Scale the horizontal axis to 100 percent of total number of species. }
  \item{type}{ Type of plot (as in function \code{\link[graphics]{plot}}) }
  \item{xlim}{ Limits for the horizontal axis. }
  \item{ylim}{ Limits for the vertical axis. }
  \item{specnames}{ Vector positions of species names to add to the rank-abundance curve. }
  \item{srt}{ The string rotation in degrees of the species names (as in \code{\link[graphics]{par}}). }

  \item{rainbow}{ Use rainbow colouring for the different curves. }
  \item{legend}{ Add the legend (you need to click in the graph where the legend needs to be plotted). }
  \item{...}{ Other arguments to be passed to functions \code{\link[graphics]{plot}} or \code{\link[graphics]{points}}. }
}

\details{
These functions provide methods of calculating and plotting rank-abundance curves.

The vertical axis can be scaled by various methods. Method "abundance" uses abundance, "proportion" uses proportional abundance (species abundance / total abundance), "logabun" calculates the logarithm of abundance using base 10 and "accumfreq" accumulates the proportional abundance.

The horizontal axis can be scaled by the total number of species, or by 100 percent of all species by option "scaledx".

The method of calculating the  confidence interval for species proportion is described in Hayek and Buzas (1997).

Functions \code{rankabundance} and \code{rankabuncomp} allow to calculate rank abundance curves for subsets of the community and environmental data sets. Function \code{rankabundance} calculates the rank abundance curve for the specified level of a selected environmental variable. Method \code{rankabuncomp} calculates the rank abundance curve for all levels of a selected environmental variable separatedly.

}

\value{
The functions provide information on rankabundance curves. Function \code{rankabundance} provides information on abundance, proportional abundance, logarithmic abundance and accumulated proportional abundance. The function also provides confidence interval limits for the proportion of each species (plower, pupper) and the proportion of species ranks (in percentage).
}

\references{
Hayek, L.-A. C. & Buzas, M.A. (1997). Surveying Natural Populations. Columbia University Press.

Kindt, R. & Coe, R. (2005)
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{http://www.worldagroforestry.org/output/tree-diversity-analysis} 
}   

\author{Roeland Kindt (World Agroforestry Centre)}

\examples{
library(vegan)
data(dune.env)
data(dune)
RankAbun.1 <- rankabundance(dune)
RankAbun.1
rankabunplot(RankAbun.1, scale='abundance', addit=FALSE, specnames=c(1,2,3))
rankabunplot(RankAbun.1, scale='logabun', addit=FALSE, specnames=c(1:30), 
    srt=45, ylim=c(1,100))
rankabuncomp(dune, y=dune.env, factor='Management', 
    scale='proportion', legend=FALSE)
## CLICK IN THE GRAPH TO INDICATE WHERE THE LEGEND NEEDS TO BE PLACED
## IF YOU OPT FOR LEGEND=TRUE.
}

\keyword{multivariate} 


