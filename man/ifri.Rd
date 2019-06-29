\name{ifri}
\alias{ifri}
\docType{data}
\title{Example data from the International Forestry Resources and Institutions (IFRI) research network}

\description{
This data set contains information on the number of stems (individuals) and basal areas for 34 vegetation plots inventoried in February 1997 in Lothlorien forest, 37 vegetation plots inventoried in February 1996 in May Creek Forest and 36 vegetation plots inventoried in May 1995 in Yellowwood State Forest. All three sites are in Indiana, USA. Data were gathered through IFRI inventory protocols to record any tree, palm and woody climber with diameter at breast height greater than or equal to 10 cm in 10-m radius circular plots; only tree species data were kept in the example data sets (IFRI research instruments and IFRI manual section P: Forest Plot Form, section D1: Tree, Palm and Woody Climber Information). 
}

\usage{data(ifri)}

\format{
  A data frame with 486 observations on the following 5 variables.
  \describe{
    \item{\code{forest}}{a factor with 3 levels: "LOT" (Lothlorien forest), "MCF" (May Creek Forest) and "YSF" (Yellowwood State Forest) }
    \item{\code{plotID}}{a factor with 107 levels providing an identification code for a 314.16 square metres (10 m radius) vegetation plot }
    \item{\code{species}}{a factor with 50 levels providing an 8 character code for a tree species }
    \item{\code{count}}{a numeric vector providing the number of stems (individuals) for each species in each vegetation plot }
    \item{\code{basal}}{a numeric vector providing the basal area (calculated from the diameter at breast height) in square cm for each species in each vegetation plot }
  }
}

\source{
IFRI (2014) Data from the International Forestry Resources and Institutions (IFRI) research network. \url{http://ifri.forgov.org/}
}

\examples{
data(ifri)
}

\keyword{datasets}
