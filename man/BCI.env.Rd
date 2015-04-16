\name{BCI.env}
\alias{BCI.env}
\docType{data}
\title{Barro Colorado Island Quadrat Descriptions}
\description{
Topography-derived variables and UTM coordinates and UTM coordinates of a 50 ha sample plot (consisting of 50 1-ha quadrats) from Barro Colorado Island of Panama. Dataset \code{\link[vegan]{BCI}} provides the tree species composition (trees with diameter at breast height equal or larger than 10 cm) of the same plots.
}
\usage{data(BCI.env)}
\format{
  A data frame with 50 observations on the following 6 variables.
  \describe{

    \item{\code{UTM.EW}}{UTM easting}
    \item{\code{UTM.NS}}{UTM northing}
    \item{\code{elevation}}{mean of the elevation values of the four cell corners}
    \item{\code{convex}}{mean elevation of the target cell minus the mean elevation of the eight surrounding cells}
    \item{\code{slope}}{mean angular deviation from horizontal of each of the four triangular planes formed by connecting three of its corners}
    \item{\code{aspectEW}}{the sine of aspect}
    \item{\code{aspectNS}}{the cosine of aspect}
  }
}

\references{
Pyke C.R., Condit R., Aguilar S. and Lao S. (2001). Floristic composition across a climatic gradient in a neotropical lowland forest. Journal of Vegetation Science 12: 553-566.

Condit R., Pitman N., Leigh E.G., Chave J., Terborgh J., Foster R.B., Nunez P., Aguilar S., Valencia R., Villa G., Muller-Landau H.C., Losos E. and Hubbell, S.P. (2002). Beta-diversity in tropical forest trees. Science 295: 666-669.

De Caceres M., P. Legendre, R. Valencia, M. Cao, L.-W. Chang, G. Chuyong, R. Condit, Z. Hao, C.-F. Hsieh, S. Hubbell, D. Kenfack, K. Ma, X. Mi, N. Supardi Noor, A. R. Kassim, H. Ren, S.-H. Su, I-F. Sun, D. Thomas, W. Ye and F. He. (2012). The variation of tree beta diversity across a global network of forest plots. Global Ecology and Biogeography 21: 1191-1202
}

\examples{
data(BCI.env)
}

\keyword{datasets}

