\name{BiodiversityRGUI}
\alias{BiodiversityRGUI}

\title{GUI for Biodiversity, Community Ecology and Ensemble Suitability Analysis }

\description{
This function provides a GUI (Graphical User Interface) for some of the functions 
of \pkg{vegan}, some other packages and some new functions to run biodiversity 
analysis, including species accumulation curves, diversity indices, Renyi 
profiles, rank-abundance curves, GLMs for analysis of species abundance and
presence-absence, distance matrices, Mantel tests, cluster and ordination 
analysis (including constrained ordination methods such as RDA, CCA, db-RDA and
CAP). In 2012 methods for ensemble suitability 
The function depends and builds on \pkg{Rcmdr}, performing all analyses
on the community and environmental datasets that the user selects. A thorough
description of the package and the biodiversity and ecological methods that it
accomodates (including examples) is provided in the freely available 
Tree Diversity Analysis manual (Kindt and Coe, 2005) that is accessible via the 
help menu.
}

\usage{
BiodiversityRGUI(changeLog = FALSE, backward.compatibility.messages = FALSE)
}

\arguments{
  \item{changeLog}{ Show the changeLog file }
  \item{backward.compatibility.messages}{ Some notes on backward compatiblity }
}

\details{
The function launches the R-Commander GUI with an extra menu for common statistical methods for biodiversity and community ecology analysis as described in the Tree Diversity Analysis manual of Roeland Kindt and Richard Coe (available via \url{https://www.worldagroforestry.org/output/tree-diversity-analysis}]) and expanded systematically with new functions that became available from the vegan community ecology package.
 
Since 2012, functions for ensemble suitability modelling were included in BiodiversityR. In 2016, a GUI was created for ensemble suitabilty modelling.

The R-Commander is launched by changing the location of the Rcmdr "etc" folder to the "etc" folder of BiodiversityR. As the files of the "etc" folder of BiodiversityR are copied from the Rcmdr, it is possible that newest versions of the R-Commander will not be launched properly. In such situations, it is possible that copying all files from the Rcmdr "etc" folder again and adding the BiodiversityR menu options to the Rcmdr-menus.txt is all that is needed to launch the R-Commander again. However, please alert Roeland Kindt about the issue.
 
BiodiversityR uses two data sets for biodiversity and community ecology analysis: the community dataset (or community matrix or species matrix) and the environmental dataset (or environmental matrix). The environmental dataset is the same dataset that is used as the "active dataset" of The R-Commander. (Note that you could sometimes use the same dataset as both the community and environmental dataset. For example, you could use the community dataset as environmental dataset as well to add information about specific species to ordination diagrams. As another example, you could use the environmental dataset as community dataset if you first calculated species richness of each site, saved this information in the environmental dataset, and then use species richness as response variable in a regression analysis.) Some options of analysis of ecological distance allow the community matrix to be a distance matrix (the community data set will be interpreted as distance matrix via \code{\link[stats:dist]{as.dist}} prior to further analysis). 

For ensemble suitability modelling, different data sets should be created and declared such as the calibration stack, the presence data set and the absence data set. The ensemble suitability modelling menu gives some guidelines on getting started with ensemble suitability modelling.
}

\value{
Besides launching the graphical user interface, the function gives some notes on backward compatibility.
}

\references{
Kindt, R. & Coe, R. (2005)
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{https://www.worldagroforestry.org/output/tree-diversity-analysis} 
}   

\author{Roeland Kindt (with some help from Jari Oksanen) }

\keyword{multivariate} 



