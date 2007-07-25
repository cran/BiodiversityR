\name{residualssurface}
\alias{residualssurface}

\title{Show and Interpolate Two Dimensional Distribution of Residuals}

\description{
This function interpolates the spatial structure of residuals of a GLM through \code{\link[mgcv]{gam}} or \code{\link[spatial]{surf.ls}} and optionally provides a graph.
}

\usage{
residualssurface(model, data, x, y, gam = F, npol = 2, plotit = T, filled = F, bubble = F)
}

\arguments{
  \item{model}{ Result of GLM as calculated by \code{\link[stats]{glm}} or \code{\link[MASS]{glm.nb}}. }
  \item{data}{ Data set that contains the spatial coordinates of the sample units used for the original model (specified as "x" and "y"). }
  \item{x}{ Horizontal position of the sample units. }
  \item{y}{ Vertical position of the sample units. }
  \item{gam}{ Interpolate the spatial structure by \code{\link[mgcv]{gam}} (if "TRUE") or by \code{\link[spatial]{surf.ls}} (if "FALSE"). }
  \item{npol}{ Degree of polynomial surface as passed to \code{\link[spatial]{surf.ls}}. }
  \item{plotit}{ Plot the interpolated surface (through \code{\link[akima]{interp}} and the residuals) . }
  \item{filled}{ Fill the contours by \code{\link[graphics]{filled.contour}}. }
  \item{bubble}{ Provide a bubble graph of the residuals: circles indicate positive residuals, whereas squares indicate negative residuals. }
}

\details{
The function reports the results of a GAM or least-squares trend surface analysis of the spatial distribution of residuals of a model (through \code{\link[stats]{residuals}}). 

Optionally, a graph is produced that can contain the trend surface, filled contours and bubble graphs in addition to the spatial location of the sample units.

}


\value{
The function reports the results of a GAM or least-squares trend surface analysis of the spatial distribution of residuals. Optionally, a graph is provided.
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
library(mgcv)
library(akima)
data(faramea)
Count.model1 <- lm(Faramea.occidentalis ~ Precipitation,
    data=faramea, na.action=na.exclude)
surface.1 <- residualssurface(Count.model1, na.omit(faramea),
    'UTM.EW', 'UTM.NS', gam=TRUE, plotit=TRUE, bubble=TRUE)
}


\keyword{multivariate} 


