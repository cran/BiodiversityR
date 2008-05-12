\name{spatialsample}
\alias{spatialsample}

\title{Spatial Sampling within a Polygon}

\description{
Spatial sampling within a polygon provides several methods of selecting rectangular sample plots within a polygon. Using a GIS package may be preferred for actual survey design.
}

\usage{
spatialsample(x,method="random",n=5,xwidth=0.5,ywidth=0.5,xleft=0,
    ylower=0,xdist=0,ydist=0,plotit=T,plothull=F)
}

\arguments{
  \item{x}{ 2-column matrix with the coordinates of the vertices of the polygon. The first column contains the horizontal (x) position, the second column contains the vertical (y) position. }
  \item{method}{ Method of sampling, any of "random", "grid" or "random grid". }
  \item{n}{ Number of sample plots to be selected, or number of horizontal and vertical grid positions. }
  \item{xwidth}{ Horizontal width of the sample plots. }
  \item{ywidth}{ Vertical width of the sample plots. }
  \item{xleft}{ Horizontal starting position of the grid. }
  \item{ylower}{ Vertical starting position of the grid. }
  \item{xdist}{ Horizontal distance between grid locations. }
  \item{ydist}{ Vertical distance between grid locations. }
  \item{plotit}{ Plot the sample plots on the current graph. }
  \item{plothull}{ Plot a convex hull around the sample plots. }
}

\details{
Spatial sampling within a polygon provides several methods of selecting the position of sample plots.

Method "random" selects random positions of the sample plots using simple random sampling.

Method "grid" selects sample plots from a grid defined by "xleft", "ylower", "xdist" and "ydist". In case \code{xdist=0} or \code{ydist=0}, then the number of grid positions are defined by "n". In case "xleft" or "ylower" are below the minimum  position of any vertix of the polygon, then a random starting position is selected for the grid.

Method "random grid" selects sample plots at random from the sampling grid using the same methods of defining the grid as for method "grid".

}

\value{
The function returns a list of centres of rectangular sample plots.
}

\references{
Kindt, R. & Coe, R. (2005) 
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{http://www.worldagroforestry.org/treesandmarkets/tree_diversity_analysis.asp} }   

\author{Roeland Kindt (World Agroforestry Centre)}

\examples{
library(splancs)
area <- array(c(10,10,15,35,40,35,5,35,35,30,30,10), dim=c(6,2))
landuse1 <- array(c(10,10,15,15,30,35,35,30), dim=c(4,2))
landuse2 <- array(c(10,10,15,15,35,30,10,30,30,35,30,15), dim=c(6,2))
landuse3 <- array(c(10,10,30,35,40,35,5,10,15,30,30,10), dim=c(6,2))
plot(area[,1], area[,2], type="n", xlab="horizontal position", 
    ylab="vertical position", lwd=2, bty="l")
polygon(landuse1)
polygon(landuse2)
polygon(landuse3)
spatialsample(area, method="random", n=20, xwidth=1, ywidth=1, plotit=TRUE, 
    plothull=FALSE)
spatialsample(area, method="grid", xwidth=1, ywidth=1, plotit=TRUE, xleft=12, 
    ylower=7, xdist=4, ydist=4)
spatialsample(area, method="random grid", n=20, xwidth=1, ywidth=1, 
    plotit=TRUE, xleft=12, ylower=7, xdist=4, ydist=4)
}

\keyword{multivariate} 


