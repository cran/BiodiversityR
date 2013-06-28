\name{ordicoeno}
\alias{ordicoeno}

\title{Coenoclines for an Ordination Axis}

\description{
A graph is produced that summarizes (through GAM as implemented by \code{\link[mgcv]{gam}}) how the abundance of all species of the community data set change along an ordination axis (based on the position of sites along the axis and the information from the community data set).
}

\usage{
ordicoeno(x, ordiplot, axis = 1, legend = FALSE, cex = 0.8, ncol = 4, ...)
}

\arguments{
  \item{x}{ Community data frame with sites as rows, species as columns and species abundance as cell values. }
  \item{ordiplot}{ Ordination plot created by \code{\link[vegan]{ordiplot}}. }
  \item{axis}{ Axis of the ordination graph (1: horizontal, 2: vertical). }
  \item{legend}{ if \code{TRUE}, then add a legend to the plot. }
  \item{cex}{ the amount by which plotting text and symbols should be magnified relative to the default; see also \code{\link[graphics]{par}} }
  \item{ncol}{ the number of columns in which to set the legend items; see also \code{\link[graphics]{legend}}  }
  \item{...}{ Other arguments passed to functions \code{\link[graphics]{plot}} and \code{\link[graphics]{points}}. }
}

\details{
This functions investigates the relationship between the species vectors and the position of sites on an ordination axis. A GAM (\code{\link[mgcv]{gam}}) investigates the relationship by using the species abundances of each species as response variable, and the site position as the explanatory variable. The graph shows how the abundance of each species changes over the gradient of the ordination axis.
}

\value{
The function plots coenoclines and provides the expected degrees of freedom (complexity of the relationship) estimated for each species by GAM.
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
library(mgcv)
data(dune)
Ordination.model1 <- rda(dune)
plot1 <- ordiplot(Ordination.model1, choices=c(1,2), scaling=1)
ordicoeno(dune, ordiplot=plot1, legend=TRUE)
}

\keyword{multivariate} 


