\name{import.from.Excel}
\alias{import.from.Excel}
\alias{import.from.Access}

\title{Import Community or Environmental Data from Excel or Access}

\description{
This function provides a method for importing community or environmental data from Excel (through functions \code{\link[RODBC:odbcConnect]{odbcConnectExcel}} or \code{\link[RODBC:odbcConnect]{odbcConnectAccess}}). In case of stacked data, a community dataset is formed by function \code{\link{makecommunitydataset}}.
}

\usage{
import.from.Excel(file = file.choose(), sheet = "community", sitenames = "sites", 
    column = "species", value = "abundance", factor = "", level = "", cepnames = FALSE)
import.from.Access(file = file.choose(), table = "community", sitenames = "sites", 
    column = "species", value = "abundance", factor = "", level = "", cepnames = FALSE)
}

\arguments{
  \item{file}{ Location of the Excel file. }
  \item{sheet}{ Name of the sheet of the Excel file to import: one of "community", "environmental" or "stacked". }
  \item{table}{ Name of the table of the Access file to import: one of "community", "environmental" or "stacked". }
  \item{sitenames}{ Name of categorical variable that provides the names for the sites. }
  \item{column}{ Name of the categorical variable for the columns of the crosstabulation (typically indicating species); passed to \code{\link{makecommunitydataset}}. }
  \item{value}{ Name of numerical variable for the cells of the crosstabulation (typically indicating abundance). The cells provide the sum of all values in the data frame; passed to \code{\link{makecommunitydataset}}. }
  \item{factor}{ Name of the variable to calculate a subset of the data frame; passed to \code{\link{makecommunitydataset}}. }
  \item{level}{ Value of the subset of the factor variable to calculate a subset of the data frame; passed to \code{\link{makecommunitydataset}}. }
  \item{cepnames}{ Should the names of columns be abbreviated via \code{\link[vegan]{make.cepnames}} (TRUE) or not (FALSE). }
}

\details{
This function provides a method of method of importing community or environmental datasets through functions \code{\link[RODBC:odbcConnect]{odbcConnectExcel}} and \code{\link[RODBC]{sqlFetch}}.

For stacked datasets, a community data set is created with function \code{\link{makecommunitydataset}}. For community data with more species than the limited number of columns in Excel, this may be the only option of importing a community dataset.

An additional advantage of the function is that the community and environmental data can be stored in the same file. 

You may want to check compatibility of the community and environmental datasets with functions \code{\link{check.datasets}} and modify the community dataset through \code{\link{same.sites}}.

}

\value{
The function returns a dataset that can be used as community or environmental data for analysis in BiodiversityR. 
}

\references{
Kindt, R. & Coe, R. (2005)
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{http://www.worldagroforestry.org/treesandmarkets/tree_diversity_analysis.asp}
 }   

\author{Roeland Kindt (World Agroforestry Centre)}

\examples{
## Example is not run as requires RODBC
## Check for file "dune.xls" in the etc folder of BiodiversityR
##
#### Default 'ordiplot3d'
library(RODBC)
\dontrun{
dunenew <- import.from.Excel(sheet="community", sitenames="sites")
dune.envnew <- import.from.Excel(sheet="environmental", sitenames="sites")
}
}

\keyword{multivariate} 


