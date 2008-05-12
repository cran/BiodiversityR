\name{BCI.env}
\alias{BCI.env}
\docType{data}
\title{Barro Colorado Island Quadrat Descriptions}
\description{
Environmental characteristics and UTM coordinates of a 50 ha sample plot (consisting of 50 1-ha quadrats) from Barro Colorado Island of Panama. Dataset \code{\link[vegan]{BCI}} provides the tree species composition (trees with diameter at breast height equal or larger than 10 cm) of the same plots.
}
\usage{data(BCI.env)}
\format{
  A data frame with 50 observations on the following 6 variables.
  \describe{
    \item{\code{UTM.EW}}{a numeric vector}
    \item{\code{UTM.NS}}{a numeric vector}
    \item{\code{Precipitation}}{a numeric vector}
    \item{\code{Elevation}}{a numeric vector}
    \item{\code{Age.cat}}{a factor with levels \code{c1} \code{c2} \code{c3}}
    \item{\code{Geology}}{a factor with levels \code{pT} \code{Tb} \code{Tbo} \code{Tc} \code{Tcm} \code{Tct} \code{Tgo} \code{Tl} \code{Tlc}}
  }
}

\source{
  \url{http://www.sciencemag.org/cgi/content/full/295/5555/666/DC1}
}
\references{
Pyke CR, Condit R, Aguilar S and Lao S. (2001). Floristic composition across a climatic gradient in a neotropical lowland forest. Journal of Vegetation Science 12: 553-566.
 
Condit, R, Pitman, N, Leigh, E.G., Chave, J., Terborgh, J., Foster,  R.B., Nunez, P., Aguilar, S., Valencia, R., Villa, G., Muller-Landau, H.C., Losos, E. & Hubbell, S.P. (2002). Beta-diversity in tropical forest trees. \emph{Science} 295: 666-669.
 
Kindt, R. & Coe, R. (2005)
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{http://www.worldagroforestry.org/treesandmarkets/tree_diversity_analysis.asp} 

}
\examples{
data(BCI.env)
}
\keyword{datasets}

