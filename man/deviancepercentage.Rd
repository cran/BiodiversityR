\name{deviancepercentage}
\alias{deviancepercentage}

\title{Calculate Percentage and Significance of Deviance Explained by a GLM}

\description{
This function calculates the percentage of deviance explained by a GLM model and calculates the significance of the model.
}

\usage{
deviancepercentage(x,data,test="F",digits=2)
}

\arguments{
  \item{x}{ Result of GLM as calculated by \code{\link[stats]{glm}} or \code{\link[MASS]{glm.nb}}. }
  \item{data}{ Data set to be used for the null model (preferably the same data set used by the 'full' model). }
  \item{test}{ Test statistic to be used for the comparison between the null model and the 'full' model as estimated by \code{\link[stats]{anova.glm}} or \code{\link[MASS]{anova.negbin}}: partial match of one of "Chisq", "F" or "Cp". }
  \item{digits}{ Number of digits in the calculation of the percentage. }
}

\details{
The function calculates the percentage of explained deviance and the significance of the 'full' model by contrasting it with the null model.

For the null model, the data is subjected to \code{\link[stats:na.fail]{na.omit}}. You should check whether the same data are used for the null and 'full' models.

}


\value{
The function calculates the percentage of explained deviance and the significance of the 'full' model by contrasting it with the null model by ANOVA. The results of the ANOVA are also provided.
}

\references{
Kindt, R. & Coe, R. (2005)
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{http://www.worldagroforestry.org/resources/databases/tree-diversity-analysis} }   

\author{Roeland Kindt}

\examples{
library(vegan)
data(dune)
data(dune.env)
dune.env$Agrsto<- dune$Agrsto
Count.model1 <- glm(Agrsto ~ Management + A1, family=quasipoisson(link=log), 
    data=dune.env, na.action=na.omit)
summary(Count.model1)
deviancepercentage(Count.model1,dune.env)
}

\keyword{multivariate} 


