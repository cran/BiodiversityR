\name{ordisymbol}
\alias{ordisymbol}
\alias{ordibubble}
\alias{ordihull.centroids}
\alias{ordispider.centroids}
\alias{ordicluster2}
\alias{ordinearest}
\alias{ordivector}

\title{Add Other Graphical Items to Ordination Diagrams}

\description{
Functions to add some other graphical itmes to ordination diagrams than provided within \pkg{vegan} by \code{\link[vegan]{ordihull}}, \code{\link[vegan:ordihull]{ordispider}}, \code{\link[vegan:ordihull]{ordiarrows}}, \code{\link[vegan:ordihull]{ordisegments}}, \code{\link[vegan:ordihull]{ordigrid}}, \code{\link[vegan:ordihull]{ordiellipse}}, \code{\link[vegan:ordihull]{ordicluster}} and \code{\link[vegan:spantree]{lines.spantree}}.
}

\usage{
ordisymbol(ordiplot, y, factor, col=1, rainbow=T, legend=T,...)
ordibubble(ordiplot,var,...)
ordihull.centroids(ordmod, ordiplot, groups, draw = c("lines", "polygon"), 
    show.groups, return.outliers=F, return.index=F, ...)
ordispider.centroids(ordmod, ordiplot, groups, w = weights(ordiplot, "sites"), 
    show.groups, return.outliers=F, return.index=F,...) 
ordicluster2(ordiplot, cluster, mingroups=1, maxgroups=nrow(ordiplot$sites),...)
ordinearest(ordiplot, dist,...)
ordivector(ordiplot, spec, lty=2,...)
}

\arguments{
  \item{ordiplot}{ An ordination graph created by \code{\link[vegan]{ordiplot}} (\pkg{vegan}). }
  \item{y}{ Environmental data frame. }
  \item{factor}{ Variable of the environmental data frame that defines subsets to be given different symbols. }
  \item{var}{ Continous variable of the environmental dataset or species from the community dataset. }
  \item{col}{ Colour (as \code{\link[graphics]{points}}). }
  \item{rainbow}{ Use rainbow colours.  }
  \item{legend}{ Add the legend (you need to click in the graph where the legend needs to be plotted). }
  \item{ordmod}{ Constrained ordination model that provides constrained ordination axes (provided by \code{\link[vegan:cca]{rda}}, \code{\link[vegan]{cca}} or \code{\link[vegan]{capscale}}).}
  \item{groups}{ Factor giving the groups for which the graphical item is drawn (as in \code{\link[vegan]{ordihull}} and \code{\link[vegan:ordihull]{ordispider}}).}
  \item{draw}{ Use either lines or polygons (as in \code{\link[vegan]{ordihull}}).}
  \item{show.groups}{ Show only given groups. This can be a vector, or TRUE if you want to show items for which condition is TRUE (as in \code{\link[vegan]{ordihull}} and \code{\link[vegan:ordihull]{ordispider}}).}
  \item{w}{ Weights used to find the average within group (as in \code{\link[vegan:ordihull]{ordispider}}).}
  \item{return.outliers}{ Return the sites that are not closest to their own centroid.}
  \item{return.index}{ Return results of a test whether sites are closest to their own centroid (TRUE) or not.}
  \item{cluster}{ Cluster object. }
  \item{mingroups}{ Minimum of clusters to be plotted. }
  \item{maxgroups}{ Maximum of clusters to be plotted.. }
  \item{dist}{ Distance matrix. }
  \item{spec}{ Species name from the community dataset. }
  \item{lty}{ Line type as specified for \code{\link[graphics]{par}}. }
  \item{...}{ Other arguments passed to functions \code{\link[graphics]{points}}, \code{\link[graphics]{symbols}}, \code{\link[vegan]{ordihull}} or \code{\link[graphics]{arrows}}. }
}

\details{
Function \code{ordisymbol} plots different levels of the specified variable in different symbols and different colours  (if \code{rainbow} option was selected). When you opt to add the legend, then you need to click in the ordination graph to indicate the position of the legend.
 
Function \code{ordibubble} draws bubble diagrams indicating the value of the specified continuous variable. Circles indicate positive values, squares indicate negative values.
 
Functions \code{ordihull.centroid} and \code{ordispider.centroid} have the same behaviour as their \pkg{vegan} equivalents (\code{\link[vegan]{ordihull}} and \code{\link[vegan:ordihull]{ordispider}}), except that hulls and spiderplots are only drawn for those sites that are closest to their own group centroid. The function thus allows to identify multivariate outliers. Optionally details about the outliers are provided. The function only returns results for ordination methods that provide constrained ordination axes.
 
Function \code{ordicluster2} provides an alternative method of overlaying information from hierarchical clustering on an ordination diagram than provided by function \code{\link[vegan:ordihull]{ordicluster}}. The method draws convex hulls around sites that are grouped into the same cluster. You can select the minimum and maximum number of clusters that are plotted (i.e. the range of clustering steps to be shown).
 
Function \code{ordinearest} draws a vector from each site to the site that is nearest to it as determined from a distance matrix. When you combine the method with \code{\link[vegan:spantree]{lines.spantree}} using the same distance measure, then you can evaluate in part how the minimum spanning tree was constructed. 
 
Function \code{ordivector} draws a vector for the specified species on the ordination diagramme and draws perpendicular lines from each site to a line that connects the origin and the head of species vector. This method helps in the biplot interpretation of a species vector as described by Jongman, ter Braak and van Tongeren (1995).  
}

\value{
These functions add graphical items to an existing ordination diagram.
}

\references{
Jongman, R.H.G, ter Braak, C.J.F & van Tongeren, O.F.R. (1987). 
    Data Analysis in Community and Landscape Ecology. Pudog, Wageningen.

Kindt, R. & Coe, R. (2005). 
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{http://www.worldagroforestry.org/treesandmarkets/tree_diversity_analysis.asp} }   

\author{Roeland Kindt (World Agroforestry Centre) and Jari Oksanen (\code{ordinearest})}

\examples{
library(vegan)
data(dune)
data(dune.env)
Ordination.model1 <- rda(dune)
plot1 <- ordiplot(Ordination.model1, choices=c(1,2), scaling=2)
ordisymbol(plot1,dune.env,"Management",legend=FALSE)
## Make sure you click in the graph to show where the legend should be placed
## In case that you choose option of legend=TRUE
plot2 <- ordiplot(Ordination.model1, choices=c(1,2), scaling=1)
distmatrix <- vegdist(dune, method='bray')
cluster <- hclust(distmatrix, method='single')
ordicluster2(plot2, cluster)
ordinearest(plot2,distmatrix,col=2)
ordivector(plot2,"Agrsto",lty=2)
}

\keyword{multivariate} 


