\name{faramea}
\alias{faramea}
\docType{data}
\title{Faramea occidentalis abundance in Panama}
\description{
This dataset describes the abundance (number of trees with diameter at breast height equal or larger than 10 cm) of the tree species Faramea occidentalis as observed in a 1-ha quadrat survey from the Barro Colorada Island of Panama. For each quadrat, some environmental characteristics are also provided.
}
\usage{data(faramea)}
\format{
  A data frame with 45 observations on the following 8 variables.
  \describe{
    \item{\code{UTM.EW}}{a numeric vector}
    \item{\code{UTM.NS}}{a numeric vector}
    \item{\code{Precipitation}}{a numeric vector}
    \item{\code{Elevation}}{a numeric vector}
    \item{\code{Age}}{a numeric vector}
    \item{\code{Age.cat}}{a factor with levels \code{c1} \code{c2} \code{c3}}
    \item{\code{Geology}}{a factor with levels \code{pT} \code{Tb} \code{Tbo} \code{Tc} \code{Tcm} \code{Tgo} \code{Tl}}
    \item{\code{Faramea.occidentalis}}{a numeric vector}
  }
}
\details{
Although the original survey documented tree species composition of all 1-ha subplots of larger (over 1 ha) sample plot, only the first (and sometimes the last) quadrats of the larger plots were included. This selection was made to avoid that larger sample plots dominated the analysis. This selection of sites is therefore different from the selection of the 50 1-ha quadrats of the largest sample plot of the same survey (\code{\link[vegan]{BCI}} and \code{\link{BCI.env}})
 
This dataset is the main dataset used for the examples provided in chapters 6 and 7 of the Tree Diversity Analysis manual (Kindt & Coe, 2005). 
}
\source{
\url{http://www.sciencemag.org/cgi/content/full/295/5555/666/DC1}
}

\references{
Pyke CR, Condit R, Aguilar S and Lao S. (2001). Floristic composition across a climatic gradient in a neotropical lowland forest. Journal of Vegetation Science 12: 553-566.
 
Condit, R, Pitman, N, Leigh, E.G., Chave, J., Terborgh, J., Foster, R.B., Nunez, P., Aguilar, S., Valencia, R., Villa, G., Muller-Landau, H.C., Losos, E. & Hubbell, S.P. (2002). Beta-diversity in tropical forest trees. \emph{Science} 295: 666-669.
 
Kindt, R. & Coe, R. (2005)
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{http://www.worldagroforestry.org/treesandmarkets/tree diversity analysis.asp} 

}
\examples{
data(faramea)
}
\keyword{datasets}

