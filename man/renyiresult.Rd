\name{renyiresult}
\alias{renyiresult}
\alias{renyiplot}
\alias{renyiaccumresult}
\alias{renyicomp}

\title{Alternative Renyi Diversity Results}

\description{
Provides some alternative methods of obtaining results on Renyi diversity profile values than provided by \code{\link[vegan]{renyi}} (\pkg{vegan}).
}

\usage{
renyiresult(x,y="",factor,level,method="all",
    scales=c(0,0.25,0.5,1,2,4,8,Inf),evenness=F,...)
renyiplot(xr,addit=F,pch=1,ylim=c(0,m),labelit=T,legend=T,col=1,cex=1,
    rainbow=T,evenness=F,...)
renyiaccumresult(x,y="",factor,level,
    scales=c(0,0.25,0.5,1,2,4,8,Inf),permutations=100,...)
renyicomp(x,y,factor,sites=Inf,
    scales=c(0,0.25,0.5,1,2,4,8,Inf),permutations=100,plotit=T,...)
}

\arguments{
  \item{x}{ Community data frame with sites as rows, species as columns and species abundance as cell values. }
  \item{y}{ Environmental data frame. }
  \item{factor}{ Variable of the environmental data frame that defines subsets to calculate diversity profiles for. }
  \item{level}{ Level of the variable to create the subset to calculate diversity profiles. }
  \item{method}{ Method of calculating the diversity profiles: "all" calculates the diversity of the entire community (all sites pooled together), "s" calculates the diversity of each site separatedly. }
  \item{scales}{ Scale parameter values as in function \code{\link[vegan]{renyi}} (\pkg{vegan}). }
  \item{evenness}{ Calculate or plot the evenness profile. }
  \item{xr}{ Result from \code{\link[vegan]{renyi}} or \code{renyiresult}. }
  \item{addit}{ Add diversity profile to an existing graph. }
  \item{pch}{ Symbol used for drawing the diversity profiles (as in function \code{\link[graphics]{points}}). }
  \item{ylim}{ Limits of the vertical axis. }
  \item{labelit}{ Provide site labels (site names) at beginning and end of the diversity profiles. }
  \item{legend}{ Add the legend (you need to click in the graph where the legend needs to be plotted). }
  \item{col}{ Colour for the diversity profile (as in function \code{\link[graphics]{points}}). }
  \item{cex}{ Character expansion factor (as in function \code{\link[graphics]{points}}).  }
  \item{rainbow}{ Use rainbow colours for the diversity profiles. }
  \item{sites}{ Maximum number of sites to provide profile values. }
  \item{permutations}{ Number of permutations for the Monte-Carlo simulations for accumulated renyi diversity profiles (estimated by \code{\link[vegan:renyi]{renyiaccum}}). }
  \item{plotit}{ Plot the results (you need to click in the graph where the legend should be plotted). }
  \item{...}{ Other arguments to be passed to functions \code{\link[vegan]{renyi}} or \code{\link[graphics]{plot}}. }
}

\details{
These functions provide some alternative methods of obtaining results with diversity profiles, although function \code{\link[vegan]{renyi}} is always used to calculate the diversity profiles.

The method of calculating the diversity profiles: "all" calculates the diversity profile of the entire community (all sites pooled together), whereas "s" calculates the diversity profile of each site separatedly. The evenness profile is calculated by subtracting the profile value at scale 0 from all the profile values.

Functions \code{renyiresult}, \code{renyiaccumresult} and \code{renyicomp} allow to calculate diversity profiles for subsets of the community and environmental data sets. functions \code{renyiresult} and \code{renyiaccumresult} calculate the diversity profiles for the specified level of a selected environmental variable. Method \code{renyicomp} calculates the diversity profile for all levels of a selected environmental variable separatedly.

Functions \code{renyicomp} and \code{renyiaccumresult} calculate accumulation curves for the Renyi diversity profile by randomised pooling of sites and calculating diversity profiles for the pooled sites as implemented in \code{\link[vegan:renyi]{renyiaccum}}. The method is similar to the random method of species accumulation (\code{\link[vegan]{specaccum}}). If the number of "sites" is not changed from the default, it is replaced by the sample size of the level with the fewest number of sites.

}

\value{
The functions provide alternative methods of obtaining Renyi diversity profiles.
}

\references{
Kindt R., Degrande A., Turyomurugyendo L., Mbosso C., Van Damme P., Simons A.J. (2001). Comparing species richness and evenness contributions to on-farm tree diversity for data sets with varying sample sizes from Kenya, Uganda, Cameroon and Nigeria with randomised diversity profiles. Paper presented at IUFRO conference on forest biometry, modeling and information science, 26-29 June, University of Greenwich, UK

Kindt R. (2002). Methodology for tree species diversification planning for African agroecosystems. Thesis submitted in fulfilment of the requirement of the degree of doctor (PhD) in applied biological sciences. Faculty of agricultural and applied biological sciences, Ghent University, Ghent (Belgium), 332+xi pp.

Kindt R., Van Damme P. & Simons A.J. (2006). Tree diversity in western Kenya: using diversity profiles to characterise richness and evenness. Biodiversity and Conservation 15: 1253-1270.

Kindt, R. & Coe, R. (2005) 
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{http://www.worldagroforestry.org/treesandmarkets/tree_diversity_analysis.asp} 

}   

\author{Roeland Kindt (World Agroforestry Centre)}

\examples{
library(vegan)
data(dune.env)
data(dune)
Renyi.1 <- renyiresult(dune, y=dune.env, factor='Management', level='NM', 
    method='s')
Renyi.1
renyiplot(Renyi.1, evenness=FALSE, addit=FALSE, pch=1,col='1', cex=1, 
    legend=FALSE)
## CLICK IN THE GRAPH TO INDICATE WHERE THE LEGEND NEEDS TO BE PLACED
## IN CASE THAT YOU OPT FOR LEGEND=TRUE
}

\keyword{multivariate} 



