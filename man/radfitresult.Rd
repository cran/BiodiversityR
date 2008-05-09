\name{radfitresult}
\alias{radfitresult}

\title{Alternative Rank Abundance Fitting Results}

\description{
Provides alternative methods of obtaining rank abundance curves than provided by functions \code{\link[vegan]{radfit}}, \code{\link[vegan]{fisherfit}} and \code{\link[vegan:fisherfit]{prestonfit}} (\pkg{vegan}), although these same functions are called.
}

\usage{
radfitresult(x,y="",factor="",level,plotit=T)
}

\arguments{
  \item{x}{ Community data frame with sites as rows, species as columns and species abundance as cell values. }
  \item{y}{ Environmental data frame. }
  \item{factor}{ Variable of the environmental data frame that defines subsets to calculate fitted rank-abundance curves for. }
  \item{level}{ Level of the variable to create the subset to calculate fitted rank-abundance curves. }
  \item{plotit}{ Plot the results obtained by \code{\link[vegan:radfit]{plot.radfit}} . }
}

\details{
These functions provide some alternative methods of obtaining fitted rank-abundance curves, although functions \code{\link[vegan]{radfit}}, \code{\link[vegan]{fisherfit}} and \code{\link[vegan:fisherfit]{prestonfit}} (\pkg{vegan}) are called to calculate the actual results.
}

\value{
The function returns the results from three methods of fitting rank-abundance curves:
  \item{radfit}{ results of \code{\link[vegan]{radfit}}. }
  \item{fisherfit}{ results of \code{\link[vegan]{fisherfit}}. }
  \item{prestonfit}{ results of \code{\link[vegan]{prestonfit}}. }

Optionally, a plot is provided of the \code{\link[vegan]{radfit}} results by \code{\link[vegan:radfit]{plot.radfit}}.
}

\references{
Kindt, R. & Coe, R. (2005)
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{http://www.worldagroforestry.org/treesandmarkets/tree_diversity_analysis.asp} }   

\author{Roeland Kindt (World Agroforestry Centre)}

\examples{
library(vegan)
data(BCI)
BCIall <- t(as.matrix(colSums(BCI)))
radfitresult(BCIall)
}

\keyword{multivariate} 


