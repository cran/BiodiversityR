\name{transfgradient}
\alias{transfgradient}
\docType{data}
\title{Gradient for Hypothetical Example of Turover of Species Composition}
\description{
This dataset documents the site sequence of 19 sites on a gradient determined from unimodal species distributions. The dataset is accompanied by \code{\link{transfspecies}} that documents the species composition of the sites. This is a hypothetical example that allows to investigate how well ecological distance measures or ordination methods recover the expected best sequence of sites.
}
\usage{data(transfgradient)}
\format{
  A data frame with 19 observations on the following variable.
  \describe{
    \item{\code{gradient}}{a numeric vector}
  }
}
\source{
Legendre, P. & Gallagher, E.D. (2001) Ecologically meaningful transformations for ordination of species data. Oecologia 129: 271-280. 
}
\references{
Figure 3a.
}
\examples{
data(transfspecies)
data(transfgradient)
plot(transfspecies[,1]~transfgradient[,1],xlab="gradient",
    ylab="species abundance",type="n",ylim=c(0.5,8.5))
for (i in 1:9) {points(transfgradient[,1],transfspecies[,i],type="o",pch=i)}
}
\keyword{datasets}
