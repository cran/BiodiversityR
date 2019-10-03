\name{importancevalue}
\alias{importancevalue}
\alias{importancevalue.comp}

\title{Importance Value}

\description{
Calculates the importance values of tree species based on frequency (calculated from number of plots), density (calculated from number of individuals) and dominance (calculated from basal area). See details.
}

\usage{
importancevalue(x, site="plotID", species="species", 
    count="count", basal="basal", 
    factor="forest", level="") 

importancevalue.comp(x, site="plotID", species="species", 
    count="count", basal="basal",
    factor="forest")

}

\arguments{
  \item{x}{ data frame with information on plot identities, species identities, number of individuals and basal areas }
  \item{site}{ factor variable providing the identities of survey plots }
  \item{species}{ factor variable providing the identities of tree species }
  \item{count}{ number of individuals for each tree species in each survey plot }
  \item{basal}{ basal area for each tree species in each survey plot }
  \item{factor}{ factor variable used to define subsets (typically different forest reserves) }
  \item{level}{ level of the factor variable used to create a subset from the original data }
}

\details{
The importance value is calculated as the sum from (i) the relative frequency; (ii) the relative density; and (iii) the relative dominance. The importance value ranges between 0 and 300.

Frequency is calculated as the number of plots where a species is observed divided by the total number of survey plots. Relative frequency is calculated by dividing the frequency by the sum of the frequencies of all species, multiplied by 100 (to obtain a percentage). 

Density is calculated as the total number of individuals of a species. Relative density is calculated by dividing the density by the sum of the densities of all species, multiplied by 100 (to obtain a percentage). 

Dominance is calculated as the total basal area of a species. Relative dominance is calculated by dividing the dominance by the sum of the dominance of all species, multiplied by 100 (to obtain a percentage). 

Functions \code{importancevalue.comp} applies function \code{importancevalue} to all available levels of a factor variable.
}

\value{
Provides information on the importance value for all tree species
}

\references{

Curtis, J.T. & McIntosh, R. P. (1951) An Upland Forest Continuum in the Prairie-Forest Border Region of Wisconsin. Ecology 32: 476-496.

Kent, M. (2011) Vegetation Description and Data Analysis: A Practical Approach. Second edition. 428 pages.

}   

\author{Roeland Kindt (World Agroforestry Centre), Peter Newton (University of Michigan)}

\seealso{
\code{\link{ifri}}
}

\examples{
data(ifri)
importancevalue(ifri, site='plotID', species='species', count='count', 
    basal='basal', factor='forest', level='YSF')
importancevalue.comp(ifri, site='plotID', species='species', count='count', 
    basal='basal', factor='forest')

# When all survey plots are the same size, importance value
# is not affected. Counts and basal areas now calculated per square metre
ifri$count <- ifri$count/314.16
ifri$basal <- ifri$basal/314.16

importancevalue(ifri, site='plotID', species='species', count='count', 
    basal='basal', factor='forest', level='YSF')
importancevalue.comp(ifri, site='plotID', species='species', count='count', 
    basal='basal', factor='forest')

# Calculate diversity profiles from importance values
imp <- importancevalue.comp(ifri, site='plotID', species='species', 
    count='count', basal='basal', factor='forest')
vals <- imp[["values"]]
for (i in 1:length(vals)) {
    imp.i <- data.frame(imp[[vals[i]]])
    name.i <- paste(vals[[i]], ".Renyi", sep="")
    imp[[name.i]] <- renyi(imp.i$importance.value)
}

# LOT more diverse
imp$LOT.Renyi - imp$MCF.Renyi
imp$LOT.Renyi - imp$YSF.Renyi

# YSF and MCF different richness and evenness
imp$YSF.Renyi - imp$MCF.Renyi


}

\keyword{multivariate} 


