\name{dist.eval}
\alias{dist.eval}
\alias{prepare.bioenv}


\title{Distance Matrix Evaluation}

\description{
Function \code{dist.eval} provides one test of a distance matrix, and then continues with \code{\link[vegan]{distconnected}} (\pkg{vegan}). Function \code{prepare.bioenv} converts selected variables to numeric variables and then excludes all categorical variables in preparation of applying \code{\link[vegan]{bioenv}} (\pkg{vegan}).
}

\usage{
dist.eval(x, dist)
prepare.bioenv(env, as.numeric = c())
}

\arguments{
  \item{x}{ Community data frame with sites as rows, species as columns and species abundance as cell values. }
  \item{env}{ Environmental data frame with sites as rows and variables as columns. }
  \item{dist}{ Method for calculating ecological distance with function \code{\link[vegan]{vegdist}}: partial match to "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "morisita", "horn" or "mountford". }
  \item{as.numeric}{ Vector with names of variables in the environmental data set to be converted to numeric variables. }

}

\details{
Function \code{dist.eval} provides two tests of a distance matrix: 

(i) The first test checks whether any pair of sites that share some species have a larger distance than any other pair of sites that do not share any species. In case that cases are found, then a warning message is given. 

(ii) The second test is the one implemented by the \code{\link[vegan]{distconnected}} function (\pkg{vegan}). The  distconnected test is only calculated for distances that calculate a value of 1 if sites share no species (i.e. not manhattan or euclidean), using the threshold of 1 as an indication that the sites do not share any species. Interpretation of analysis of distance matrices that provided these warnings should be cautious.

Function \code{prepare.bioenv} provides some simple methods of dealing with categorical variables prior to applying \code{\link[vegan]{bioenv}}.

}

\value{
The function tests whether distance matrices have some desirable properties and provide warnings if this is not the case.
}

\references{
Kindt, R. & Coe, R. (2005)
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{http://www.worldagroforestry.org/output/tree-diversity-analysis}
 }   

\author{Roeland Kindt (World Agroforestry Centre)}

\examples{
library(vegan)
data(dune)
dist.eval(dune,"euclidean")
dist.eval(dune,"bray")

\dontrun{
data(dune.env)
dune.env2 <- dune.env[,c('A1', 'Moisture', 'Manure')]
dune.env2$Moisture <- as.numeric(dune.env2$Moisture)
dune.env2$Manure <- as.numeric(dune.env2$Manure)
sol <- bioenv(dune ~ A1 + Moisture + Manure, dune.env2)
sol
summary(sol)
dune.env3 <- prepare.bioenv(dune.env, as.numeric=c('Moisture', 'Manure'))
bioenv(dune, dune.env3)
}


}

\keyword{multivariate} 


