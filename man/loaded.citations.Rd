\name{loaded.citations}
\alias{loaded.citations}

\title{Give Citation Information for all Loaded Packages}

\description{
This function provides citation information for all loaded packages.
}

\usage{
loaded.citations()
}

\details{
The function checks for the loaded packages via \code{\link[base:zpackages]{.packages}}. Citation information is provided for the base package and for all the non-standard packages via \code{\link[utils]{citation}}.
}

\value{
The function provides a list of all loaded packages and the relevant citation information.
}

\author{Roeland Kindt (World Agroforestry Centre) }

\keyword{multivariate} 


