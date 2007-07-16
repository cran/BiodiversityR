\name{nnetrandom}
\alias{nnetrandom}

\title{Calculate the NNET Result with the Smallest Value from Various Random Starts}

\description{
This function provides the best solution from various calls to the \code{\link[nnet]{nnet}} feed-forward artificial neural networks function (\pkg{nnet}).
}

\usage{
nnetrandom(formula,data,tries=10,leave.one.out=F,...)
}

\arguments{
  \item{formula}{ Formula as passed to \code{\link[nnet]{nnet}}. }
  \item{data}{ Data as passed to \code{\link[nnet]{nnet}}. }
  \item{tries}{ Number of calls to \code{\link[nnet]{nnet}} to obtain the best solution. }
  \item{leave.one.out}{ Calculate leave-one-out predictions. }
  \item{...}{ Other arguments passed to \code{\link[nnet]{nnet}}. }
}

\details{
This function makes various calls to \code{\link[nnet]{nnet}}. If desired by the user, leave-one-out statistics are provided that report the prediction if one particular sample unit was not used for iterating the networks.
}

\value{
The function returns the same components as \code{\link[nnet]{nnet}}, but adds the following components:
  
\item{range}{ Summary of the observed "values". }
\item{tries}{ Number of different attempts to iterate an ANN. }
\item{CV}{ Predicted class when not using the respective sample unit for iterating ANN. }
\item{succesful}{ Test whether leave-one-out statistics provided the same class as the original class. }
}

\author{Roeland Kindt (World Agroforestry Centre)}

\examples{
data(faramea)
faramea <- na.omit(faramea)
faramea$presence <- as.numeric(faramea$Faramea.occidentalis > 0)
attach(faramea)
library(nnet)
result <- nnetrandom(presence ~ Elevation, data=faramea, size=2, 
    skip=FALSE, entropy=TRUE, trace=FALSE, maxit=1000, tries=1000, 
    leave.one.out=FALSE)
summary(result)
result$fitted.values
result$value
result2 <- nnetrandom(presence ~ Elevation, data=faramea, size=2, 
    skip=FALSE, entropy=TRUE, trace=FALSE, maxit=1000, tries=50, 
    leave.one.out=TRUE)
result2$range
result2$CV
result2$successful
}

\keyword{multivariate} 


