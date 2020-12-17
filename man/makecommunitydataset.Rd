\name{makecommunitydataset}
\alias{makecommunitydataset}
\alias{stackcommunitydataset}

\title{Make a Community Dataset from a Stacked Dataset}

\description{
Makes a community data set from a stacked dataset (with separate variables for the site identities, the species identities and the abundance).
}

\usage{
makecommunitydataset(x, row, column, value, factor="", level="", drop=F)
stackcommunitydataset(comm, remove.zeroes=FALSE, order.sites=FALSE, order.species=FALSE)
}

\arguments{
  \item{x}{ Data frame. }
  \item{row}{ Name of the categorical variable for the rows of the crosstabulation (typically indicating sites) }
  \item{column}{ Name of the categorical variable for the columns of the crosstabulation (typically indicating species) }
  \item{value}{ Name of numerical variable for the cells of the crosstabulation (typically indicating abundance). The cells provide the sum of all values in the data frame. }
  \item{factor}{ Name of the variable to calculate a subset of the data frame. }
  \item{level}{ Value of the subset of the factor variable to calculate a subset of the data frame. }
  \item{drop}{ Drop rows without species (species with total abundance of zero are always dropped) }

  \item{comm}{ Community data set }
  \item{remove.zeroes}{ Should rows with zero abundance be removed? }
  \item{order.sites}{ Should sites be ordered alphabetically? }
  \item{order.species}{ Should species be ordered alphabetically? }

}

\details{
\code{makecommunitydataset} calculates a cross-tabulation from a data frame, summing up all the values of the numerical variable identified as variable for the cell values. If factor="", then no subset is calculated from the data frame in the first step.

\code{stackcommunitydataset} reverses the actions of \code{makecommunitydataset} and recreates the data in stacked format. 
}

\value{
The function provides a community dataset from another data frame.
}

\references{
Kindt, R. & Coe, R. (2005) 
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{https://www.worldagroforestry.org/output/tree-diversity-analysis}
}   

\author{Roeland Kindt (World Agroforestry Centre)}

\examples{

\dontrun{
dune.file <- normalizePath(paste(system.file(package="BiodiversityR"), 
    '/etc/dunestacked.csv', sep=''))
dune.stacked <- read.csv(dune.file)

# dune.stacked has different variables for sites, species and abundance
head(dune.stacked)
dune.comm2 <- makecommunitydataset(dune.stacked, row='sites', column='species', 
    value='abundance')

# recreate the original stack
dune.stacked2 <- stackcommunitydataset(dune.comm2, remove.zeroes=T)

}

}

\keyword{multivariate} 


