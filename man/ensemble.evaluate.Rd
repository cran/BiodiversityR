\name{ensemble.evaluate}
\alias{ensemble.evaluate}
\alias{ensemble.SEDI}

\title{
Model evaluation including True Skill Statistic (TSS), AUCdiff and Symmetric Extremal Dependence Index (SEDI).
}

\description{
The main function of \code{ensemble.evaluate} calculates various model evaluation statistics. Function \code{ENSEMBLE.SEDI} calculates the Symmetric Extremal Dependence Index (SEDI) from the True Positive Rate (TPR = Sensitivity = Hit Rate) and the False Positive Rate (FPR = False Alarm Rate = 1 - Specificity).
}

\usage{

ensemble.evaluate(eval, fixed.threshold = NULL, eval.train = NULL) 

ensemble.SEDI(TPR, FPR, small = 1e-9) 

}

\arguments{
  \item{eval}{ ModelEvaluation object (\code{\link[dismo]{evaluate}}), ideally obtained via model testing data that were not used for calibrating the model. }
  \item{fixed.threshold}{ Absence-presence threshold to create the confusion matrix. See also (\code{\link[dismo]{threshold}} and \code{\link{ensemble.threshold}}). }
  \item{eval.train}{ ModelEvaluation object (\code{\link[dismo]{evaluate}}), ideally obtained via model calibration data that were used for calibrating the model. }

  \item{TPR}{ True Presence Rate, equal to correctly predicted presence observations divided by total number of presence observations. Also known as Sensitivity or Hit Rate. }
  \item{FPR}{ False Presence Rate, equal to wrongly predicted absence observations divided by total number of absence observations. Also known as False Alarm Rate. }
  \item{small}{ small amount that replaces zeroes in calculations. }
}

\details{
Details for the True Skill Statistic (TSS = TPR + TNR - 1 = TPR - FPR), Symmetric Extremal Dependence Index (SEDI), False Negative Rate (omission or miss rate) and AUCdiff (AUCtrain - AUCtest) are available from Ferro and Stephenson (2011), Wunderlich et al. (2019) and Castellanos et al. (2019). 

Values for TSS and SEDI are given for the fixed absence-presence threshold, as well as their maximal values across the entire range of candidate threshold values calculate by \code{\link[dismo]{evaluate}}.

In case that \code{fixed.threshold} is not provided, it is calculated from the calibration ModelEvaluation as the threshold that maximizes the sum of TPR (sensitivity) and TNR (specificity) (and thus also maximizes the TSS for the calibration). 
}

\value{
A numeric vector with following values.

AUC: Area Under The Curve for the testing ModelEvaluation
TSS: maximum value of the True Skill Statistic over range of threshold values
SEDI: maximum value of the Symmetric Extremal Dependence Index over range of threshold values
TSS.fixed: True Skill Statistic at the fixed threshold value
SEDI.fixed: SEDI at the fixed threshold value
FNR.fixed: False Negative Rate (= omission rate) at the fixed threshold value
MCR.fixed: Missclassification Rate at the fixed threshold value
AUCdiff: Difference between AUC for calibration and the testing data
}

\references{

Ferro CA, Stephenson DB. 2011. Extremal Dependence Indices: Improved Verification Measures for Deterministic Forecasts of Rare Binary Events. Wea. Forecasting 26: 699-713. 

Wunderlich RF, Lin Y-P, Anthony J, Petway JR. 2019. Two alternative evaluation metrics to replace the true skill statistic in the assessment of species distribution models. Nature Conservation 35: 97-116. \url{https://doi.org/10.3897/natureconservation.35.33918}

Castellanos AA, Huntley JW, Voelker G, Lawing AM. 2019. Environmental filtering improves ecological niche models across multiple scales. \url{https://doi.org/10.1111/2041-210X.13142}

Kindt R. 2018. Ensemble species distribution modelling with transformed suitability values. Environmental Modelling & Software 100: 136-145. \url{https://doi.org/10.1016/j.envsoft.2017.11.009}
}

\author{ Roeland Kindt (World Agroforestry Centre)}

\seealso{
\code{\link{ensemble.batch}}
}

\examples{

## check examples from Ferro and Stephenson (2011)
## see their Tables 2 - 5

TPR.Table2 <- 55/100
FPR.Table2 <- 45/900
ensemble.SEDI(TPR=TPR.Table2, FPR=FPR.Table2)

TPR.Table4 <- 195/300
FPR.Table4 <- 105/700
ensemble.SEDI(TPR=TPR.Table4, FPR=FPR.Table4)

\dontrun{
## Not run: 
# get predictor variables
library(dismo)
predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
    pattern='grd', full.names=TRUE)
predictors <- stack(predictor.files)
# subset based on Variance Inflation Factors
predictors <- subset(predictors, subset=c("bio5", "bio6", 
    "bio16", "bio17", "biome"))
predictors
predictors@title <- "predictors"

# presence points
presence_file <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
pres <- read.table(presence_file, header=TRUE, sep=',')[,-1]

# the kfold function randomly assigns data to groups; 
# groups are used as calibration (1/4) and training (3/4) data
groupp <- kfold(pres, 4)
pres_train <- pres[groupp !=  1, ]
pres_test <- pres[groupp ==  1, ]

# choose background points
background <- randomPoints(predictors, n=1000, extf=1.00)
colnames(background)=c('lon', 'lat')
groupa <- kfold(background, 4)
backg_train <- background[groupa != 1, ]
backg_test <- background[groupa == 1, ]

# formulae for random forest and generalized linear model
# compare with: ensemble.formulae(predictors, factors=c("biome"))

rfformula <- as.formula(pb ~ bio5+bio6+bio16+bio17)

glmformula <- as.formula(pb ~ bio5 + I(bio5^2) + I(bio5^3) + 
    bio6 + I(bio6^2) + I(bio6^3) + bio16 + I(bio16^2) + I(bio16^3) + 
    bio17 + I(bio17^2) + I(bio17^3) )

# fit four ensemble models (RF, GLM, BIOCLIM, DOMAIN)
# factors removed for BIOCLIM, DOMAIN, MAHAL
ensemble.nofactors <- ensemble.calibrate.models(x=predictors, p=pres_train, a=backg_train, 
    pt=pres_test, at=backg_test,
    species.name="Bradypus",
    ENSEMBLE.tune=TRUE,
    ENSEMBLE.min = 0.65,
    MAXENT=0, MAXNET=0, MAXLIKE=0, GBM=0, GBMSTEP=0, RF=1, CF=0, 
    GLM=1, GLMSTEP=0, GAM=0, GAMSTEP=0, MGCV=0, MGCVFIX=0, 
    EARTH=0, RPART=0, NNET=0, FDA=0, SVM=0, SVME=0, GLMNET=0,
    BIOCLIM.O=0, BIOCLIM=1, DOMAIN=1, MAHAL=0, MAHAL01=0,
    Yweights="BIOMOD",
    factors="biome",
    evaluations.keep=TRUE, models.keep=FALSE,
    RF.formula=rfformula,
    GLM.formula=glmformula)

# with option evaluations.keep, all model evaluations are saved in the ensemble object
attributes(ensemble.nofactors$evaluations)

# Get evaluation statistics for the ENSEMBLE model
eval.ENSEMBLE <- ensemble.nofactors$evaluations$ENSEMBLE.T
eval.calibrate.ENSEMBLE <- ensemble.nofactors$evaluations$ENSEMBLE.C
ensemble.evaluate(eval=eval.ENSEMBLE, eval.train=eval.calibrate.ENSEMBLE)

# TSS is maximum where specificity + sensitivity is maximum
threshold.specsens <- threshold(eval.ENSEMBLE, stat="spec_sens")
ensemble.evaluate(eval=eval.ENSEMBLE, fixed.threshold=threshold.specsens,
    eval.train=eval.calibrate.ENSEMBLE)

# usual practice to calculate threshold from calibration data
ensemble.evaluate(eval=eval.ENSEMBLE, eval.train=eval.calibrate.ENSEMBLE)

}


}

