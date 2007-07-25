\name{transfspecies}
\alias{transfspecies}
\docType{data}
\title{Hypothetical Example of Turover of Species Composition}
\description{
This dataset documents the species composition of 19 sites that follow a specific sequence of sites as determined from unimodal species distributions. The dataset is accompanied by \code{\link{transfgradient}} that documents the gradient in species turnover. This is a hypothetical example that allows to investigate how well ecological distance measures or ordination methods recover the expected best sequence of sites.
}
\usage{data(transfspecies)}
\format{
  A data frame with 19 observations on the following 9 variables.
  \describe{
    \item{\code{species1}}{a numeric vector}
    \item{\code{species2}}{a numeric vector}
    \item{\code{species3}}{a numeric vector}
    \item{\code{species4}}{a numeric vector}
    \item{\code{species5}}{a numeric vector}
    \item{\code{species6}}{a numeric vector}
    \item{\code{species7}}{a numeric vector}
    \item{\code{species8}}{a numeric vector}
    \item{\code{species9}}{a numeric vector}
  }
}
\details{
The example in the Tree Diversity Analysis manual only looks at the ecological distance from the first site. Hence, only the first 10 sites that share some species with this site should be selected.
 
This dataset enables investigations of how well ecological distance measures and ordination diagrams reconstruct the gradient (sequence of sites). The gradient expresses how the sites would be arranged based on their species composition.
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
