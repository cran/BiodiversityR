\name{BiodiversityR-package}
\alias{BiodiversityR-package}
\alias{BiodiversityR}
\docType{package}
\title{
GUI for biodiversity and community ecology analysis
}
\description{
This package provides a GUI (Graphical User Interface, via the R-Commander) and some utility functions (often based on the vegan package) for statistical analysis of biodiversity and ecological communities, including species accumulation curves, diversity indices, Renyi profiles, GLMs for analysis of species abundance and presence-absence, distance matrices, Mantel tests, and cluster, constrained and unconstrained ordination analysis. A book on biodiversity and community ecology analysis is available for free download from the website.
}
\details{
\tabular{ll}{
Package: \tab BiodiversityR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2007-4-2\cr
License: \tab GPL v2\cr
Depends: \tab R (>= 2.4.0), vegan (>= 1.9-31)\cr
Suggests: \tab MASS, mgcv, cluster, Rcmdr (>= 1.3-4), car, RODBC, rpart, effects, multcomp, ellipse, maptree\cr
Packaged: \tab Wed Jul 18 18:21:26 2007; RKINDT\cr
Built: \tab R 2.4.0; ; 2007-07-18 18:21:28; windows\cr
}

We warmly thank all that provided inputs that lead to improvement of the Tree Diversity Analysis manual that describes common methods for biodiversity and community ecology analysis and its accompanying software. We especially appreciate the comments received during training sessions with draft versions of this manual and the accompanying software in Kenya, Uganda and Mali. We are equally grateful to the thoughtful reviews by Dr Simoneta Negrete-Yankelevich (Instituto de Ecologia, Mexico) and Dr Robert Burn (Reading University, UK) of the draft version of this manual, and to Hillary Kipruto for help in editing of this manual. We also want to specifically thank Mikkel Grum, Jane Poole and Paulo van Breugel for helping in testing the packaged version of the software. We also want to give special thanks for all the support that was given by Jan Beniest, Tony Simons and Kris Vanhoutte in realizing the book and software. 
 
We highly appreciate the support of the Programme for Cooperation with International Institutes (SII), Education and Development Division of the Netherlands Ministry of Foreign Affairs, and VVOB (The Flemish Association for Development Cooperation and Technical Assistance, Flanders, Belgium) for funding the development for this manual. We also thank VVOB for seconding Roeland Kindt to the World Agroforestry Centre (ICRAF). The tree diversity analysis manual was inspired by research, development and extension activities that were initiated by ICRAF on tree and landscape diversification. We want to acknowledge the various donor agencies that have funded these activities, especially VVOB, DFID, USAID and EU.
 
We are grateful for the developers of the R Software for providing a free and powerful statistical package that allowed development of BiodiversityR. We also want to give special thanks to Jari Oksanen for developing the vegan package and John Fox for developing the Rcmdr package, which are key packages that are used by BiodiversityR.


Index:
\preformatted{
BCI.env                 Barro Colorado Island Quadrat Descriptions
BiodiversityRGUI        GUI for Biodiversity Analysis and Ordination
CAPdiscrim              Canonical Analysis of Principal Coordinates
                        based on Discriminant Analysis
NMSrandom               Calculate the NMS Result with the Smallest
                        Stress from Various Random Starts
PCAsignificance         PCA Significance
accumresult             Alternative Species Accumulation Curve Results
add.spec.scores         Add Species Scores to Unconstrained Ordination
                        Results
crosstabanalysis        Presence-absence Analysis by Cross Tabulation
deviancepercentage      Calculate Percentage and Significance of
                        Deviance Explained by a GLM
dist.eval               Distance Matrix Evaluation
dist.zeroes             Distance Matrix Transformation
distdisplayed           Compare Distance Displayed in Ordination
                        Diagram with Distances of Distance Matrix
disttransform           Community Matrix Transformation
diversityresult         Alternative Diversity Results
faramea                 Faramea occidentalis abundance in Panama
import.from.Excel       Import Community or Environmental Data from
                        Excel or Access
loaded.citations        Give Citation Information for all Loaded
                        Packages
makecommunitydataset    Make a Community Dataset from a Stacked Dataset
multiconstrained        Pairwise Comparisons for All Levels of a
                        Categorical Variable by RDA, CCA or Capscale
nnetrandom              Calculate the NNET Result with the Smallest
                        Value from Various Random Starts
ordicoeno               Coenoclines for an Ordination Axis
ordisymbol              Add Other Graphical Items to Ordination
                        Diagrams
radfitresult            Alternative Rank Abundance Fitting Results
rankabundance           Rank Abundance Curves
removeNAcomm            Synchronize Community and Environmental
                        Datasets
renyiresult             Alternative Renyi Diversity Results
residualssurface        Show and Interpolate Two Dimensional
                        Distribution of Residuals
spatialsample           Spatial Sampling within a Polygon
transfgradient          Gradient for Hypothetical Example of Turover of
                        Species Composition
transfspecies           Hypothetical Example of Turover of Species
                        Composition
}
}
\author{
Roeland Kindt

Maintainer: Roeland Kindt (World Agroforestry Centre [ICRAF])
<R.Kindt@CGIAR.org>
 
}
\references{
Kindt, R. & Coe, R. (2005) Tree diversity analysis: A manual and software for common statistical methods for ecological and biodiversity studies. 
 
\url{http://www.worldagroforestry.org/treesandmarkets/tree_diversity_analysis.asp} 
 
We suggest to use this citation for this software as well (together with citations of all other packages that were used)  
}
\keyword{ package }



