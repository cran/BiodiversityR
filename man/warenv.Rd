\name{warenv}
\alias{warenv}
\docType{data}
\title{Warburgia ugandensis Population Structure}
\description{
This data set contains population and regional locations for 100 individuals of the Warburgia ugandensis tree species (a medicinal tree species native to Eastern Africa). This data set is associated with \code{warcom} that contains scores for 185 AFLP loci.
}
\usage{data(warenv)}
\format{
  A data frame with 100 observations on the following 4 variables.
  \describe{
    \item{\code{population}}{a factor with levels \code{Kibale} \code{Kitale} \code{Laikipia} \code{Lushoto} \code{Mara}}
    \item{\code{popshort}}{a factor with levels \code{KKIT} \code{KLAI} \code{KMAR} \code{TLUS} \code{UKIB}}
    \item{\code{country}}{a factor with levels \code{Kenya} \code{Tanzania} \code{Uganda}}
    \item{\code{rift.valley}}{a factor with levels \code{east} \code{west}}
  }
}

\source{
Muchugi, A.N. (2007) Population genetics and taxonomy of important medicinal tree species of the genus Warburgia. PhD Thesis. Kenyatta University, Kenya.
}

\examples{
data(warenv)
}
\keyword{datasets}
