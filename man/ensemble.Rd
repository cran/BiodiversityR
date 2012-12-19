\name{ensemble.test}
\alias{ensemble.test}
\alias{ensemble.test.splits}
\alias{ensemble.test.gbm}
\alias{ensemble.test.nnet}
\alias{ensemble.formulae}

\title{
Suitability mapping based on ensembles of modelling algorithms: comparison of different algorithms
}

\description{
The basic function \code{ensemble.test} allows to evaluate different algorithms for (species) suitability modelling, including maximum entropy (MAXENT), boosted regression trees, random forests, generalized linear models (including stepwise selection of explanatory variables), generalized additive models (including stepwise selection of explanatory variables), multivariate adaptive regression splines, regression trees, artificial neural networks, flexible discriminant analysis, support vector machines, the BIOCLIM algorithm, the DOMAIN algorithm and the Mahalanobis algorithm. These sets of functions were developed in parallel with the \code{biomod2} package, especially for inclusion of the maximum entropy algorithm, but also to allow for a more direct integration with the BiodiversityR package, more direct handling of model formulae and greater focus on mapping. Researchers and students of species distribution are strongly encouraged to familiarize themselves with all the options of the BIOMOD and dismo packages.
}


\usage{
ensemble.test(x, p, a = NULL, an = 1000, ext = NULL, k = 0, pt = NULL, at = NULL, PLOTS = TRUE, evaluations.keep = FALSE, 
    MAXENT = 1, GBM = 1, GBMSTEP = 1, RF = 1, GLM = 1, GLMSTEP = 1, 
    GAM = 1, GAMSTEP = 1, MGCV = 1, EARTH = 1, RPART = 1, NNET = 1, FDA = 1, SVM=1, BIOCLIM = 1, DOMAIN = 1, MAHAL=1, 
    Yweights = "BIOMOD", factors = NULL, formulae.defaults=TRUE,
    GBM.formula = NULL, GBM.n.trees = 3000,    
    GBMSTEP.gbm.x = 2:(1 + nlayers(x)), GBMSTEP.tree.complexity = 5, GBMSTEP.learning.rate = 0.01, GBMSTEP.bag.fraction = 0.5, 
    RF.formula = NULL, RF.ntree = 750, RF.mtry = log(nlayers(x)), 
    GLM.formula = NULL, GLM.family = binomial(link = "logit"), 
    GLMSTEP.steps = 1000, STEP.formula = NULL, GLMSTEP.scope = NULL, GLMSTEP.k = 2, 
    GAM.formula = NULL, GAM.family = binomial(link = "logit"), 
    GAMSTEP.steps = 1000, GAMSTEP.scope = NULL, 
    MGCV.formula = NULL, 
    EARTH.formula = NULL, EARTH.glm = list(family = binomial(link = "logit"), maxit=1000), 
    RPART.formula = NULL, RPART.xval = 50, 
    NNET.formula = NULL, NNET.size = 8, NNET.decay = 0.01, 
    FDA.formula = NULL, SVM.formula=NULL)

ensemble.test.splits(x, p, a = NULL, an = 1000, ext = NULL, k = 5, digits = 2, PLOTS = FALSE, SCRIPT = TRUE, 
    MAXENT = 1, GBM = 1, GBMSTEP = 1, RF = 1, GLM = 1, GLMSTEP = 1, GAM = 1, 
    GAMSTEP = 1, MGCV = 1, EARTH = 1, RPART = 1, NNET = 1, FDA = 1, SVM=1, BIOCLIM = 1, DOMAIN = 1, MAHAL=1, 
    Yweights = "BIOMOD", factors = NULL, formulae.defaults=TRUE,
    GBM.formula = NULL, GBM.n.trees = 3000,    
    GBMSTEP.gbm.x = 2:(1 + nlayers(x)), GBMSTEP.tree.complexity = 5, GBMSTEP.learning.rate = 0.01, GBMSTEP.bag.fraction = 0.5, 
    RF.formula = NULL, RF.ntree = 750, RF.mtry = log(nlayers(x)), 
    GLM.formula = NULL, GLM.family = binomial(link = "logit"), 
    GLMSTEP.steps = 1000, STEP.formula = NULL, GLMSTEP.scope = NULL, GLMSTEP.k = 2, 
    GAM.formula = NULL, GAM.family = binomial(link = "logit"), 
    GAMSTEP.steps = 1000, GAMSTEP.scope = NULL, 
    MGCV.formula = NULL, 
    EARTH.formula = NULL, EARTH.glm=list(family=binomial(link="logit"), maxit=1000), 
    RPART.formula = NULL, RPART.xval = 50, 
    NNET.formula = NULL, NNET.size = 8, NNET.decay = 0.01, 
    FDA.formula = NULL, SVM.formula=NULL)

ensemble.test.gbm(x, p, a = NULL, an = 1000, ext = NULL, k = 5, digits = 2, PLOTS = FALSE, Yweights = "BIOMOD", factors=NULL, 
    GBMSTEP.gbm.x = 2:(1 + nlayers(x)), complexity = c(3:6), learning = c(0.005, 0.002, 0.001), GBMSTEP.bag.fraction = 0.5)

ensemble.test.nnet(x, p, a = NULL, an = 1000, ext = NULL, k = 5, digits = 2, PLOTS = FALSE, Yweights = "BIOMOD", factors=NULL, 
    formulae.defaults=TRUE, NNET.formula = NULL, 
    size = c(2, 4, 6, 8), decay = c(0.1, 0.05, 0.01, 0.001))

ensemble.formulae(x, factors = NULL, dummy.vars = NULL)
}

\arguments{
  \item{x}{ RasterStack object (\code{\link[raster]{stack}}) containing all layers that correspond to explanatory variables  }
  \item{p}{ presence points used for calibrating the suitability models, typically available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{a}{ background points used for calibrating the suitability models, typically available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{an}{ number of background points for calibration to be selected with \code{\link[dismo]{randomPoints}} in case argument \code{a} is missing}
  \item{ext}{ an Extent object to limit the selection of background points to a sub-region of \code{x}, typically provided as c(lonmin, lonmax, latmin, latmax); see also \code{\link[dismo]{randomPoints}} and \code{\link[raster]{extent}} }
  \item{k}{ If larger than 1, the number of groups to split between calibration (k-1) and evaluation (1) data sets (for example, \code{k=5} results in 4/5 of presence and background points to be used for calibrating the models, and 1/5 of presence and background points to be used for evaluating the models). See also \code{\link[dismo]{kfold}}.  }
  \item{pt}{ presence points used for evaluating the suitability models, available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{at}{ background points used for evaluating the suitability models, available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{PLOTS}{ Plot the results of evaluation for the various suitability models (if \code{TRUE}). See also \code{\link[dismo]{evaluate}}.  }
  \item{evaluations.keep}{ Keep the results of evaluations (if \code{TRUE}). See also \code{\link[dismo]{evaluate}}. }
  \item{MAXENT}{ number: if larger than 0, then a maximum entropy model (\code{\link[dismo]{maxent}}) will be fitted among ensemble }
  \item{GBM}{ number: if larger than 0, then a boosted regression trees model (\code{\link[gbm]{gbm}}) will be fitted among ensemble }
  \item{GBMSTEP}{ number: if larger than 0, then a stepwise boosted regression trees model (\code{\link[dismo]{gbm.step}}) will be fitted among ensemble }
  \item{RF}{ number: if larger than 0, then a random forest model (\code{\link[randomForest]{randomForest}}) will be fitted among ensemble }
  \item{GLM}{ number: if larger than 0, then a generalized linear model (\code{\link[stats]{glm}}) will be fitted among ensemble }
  \item{GLMSTEP}{ number: if larger than 0, then a stepwise generalized linear model (\code{\link[MASS]{stepAIC}}) will be fitted among ensemble }
  \item{GAM}{ number: if larger than 0, then a generalized additive model (\code{\link[gam]{gam}}) will be fitted among ensemble }
  \item{GAMSTEP}{ number: if larger than 0, then a stepwise generalized additive model (\code{\link[gam]{step.gam}}) will be fitted among ensemble }
  \item{MGCV}{ number: if larger than 0, then a generalized additive model (\code{\link[mgcv]{gam}}) will be fitted among ensemble }
  \item{EARTH}{ number: if larger than 0, then a multivariate adaptive regression spline model (\code{\link[earth]{earth}}) will be fitted among ensemble }
  \item{RPART}{ number: if larger than 0, then a recursive partioning and regression tree model (\code{\link[rpart]{rpart}}) will be fitted among ensemble }
  \item{NNET}{ number: if larger than 0, then an artificial neural network model (\code{\link[nnet]{nnet}}) will be fitted among ensemble }
  \item{FDA}{ number: if larger than 0, then a flexible discriminant analysis model (\code{\link[mda]{fda}}) will be fitted among ensemble }
  \item{SVM}{ number: if larger than 0, then a support vector machine model (\code{\link[kernlab]{ksvm}}) will be fitted among ensemble }
  \item{BIOCLIM}{ number: if larger than 0, then the BIOCLIM algorithm (\code{\link[dismo]{bioclim}}) will be fitted among ensemble }
  \item{DOMAIN}{ number: if larger than 0, then the DOMAIN algorithm (\code{\link[dismo]{domain}}) will be fitted among ensemble }
  \item{MAHAL}{ number: if larger than 0, then the Mahalanobis algorithm (\code{\link[dismo]{mahal}}) will be fitted among ensemble }
  \item{Yweights}{ chooses how cases of presence and background (absence) are weighted; \code{"BIOMOD"} results in equal weighting of all presence and all background cases, \code{"equal"} results in equal weighting of all cases. The user can supply a vector of weights similar to the number of cases in the calibration data set. }
  \item{factors}{ vector that indicates which variables are factors; see also \code{\link[dismo]{prepareData}}  }
  \item{formulae.defaults}{ Suggest formulae for most of the models (if \code{TRUE}). See also \code{\link{ensemble.formulae}}.   }
  \item{GBM.formula}{ formula for the boosted regression trees algorithm; see also \code{\link[gbm]{gbm}}  }
  \item{GBM.n.trees}{ total number of trees to fit for the boosted regression trees model; see also \code{\link[gbm]{gbm}}  }
  \item{GBMSTEP.gbm.x}{ indices of column numbers with explanatory variables for stepwise boosted regression trees; see also \code{\link[dismo]{gbm.step}}  }
  \item{GBMSTEP.tree.complexity}{ complexity of individual trees for stepwise boosted regression trees; see also \code{\link[dismo]{gbm.step}}  }
  \item{GBMSTEP.learning.rate}{ weight applied to individual trees for stepwise boosted regression trees; see also \code{\link[dismo]{gbm.step}} }
  \item{GBMSTEP.bag.fraction}{ proportion of observations used in selecting variables for stepwise boosted regression trees; see also \code{\link[dismo]{gbm.step}}  }
  \item{RF.formula}{ formula for random forest algorithm; see also \code{\link[randomForest]{randomForest}}  }
  \item{RF.ntree}{ number of trees to grow for random forest algorithm; see also \code{\link[randomForest]{randomForest}}  }
  \item{RF.mtry}{ number of variables randomly sampled as candidates at each split for random forest algorithm; see also \code{\link[randomForest]{randomForest}}  }
  \item{GLM.formula}{ formula for the generalized linear model; see also \code{\link[stats]{glm}}  }
  \item{GLM.family}{ description of the error distribution and link function for the generalized linear model; see also \code{\link[stats]{glm}}  }
  \item{GLMSTEP.steps}{ maximum number of steps to be considered for stepwise generalized linear model; see also \code{\link[MASS]{stepAIC}}  }
  \item{STEP.formula}{ formula for the "starting model" to be considered for stepwise generalized linear model; see also \code{\link[MASS]{stepAIC}}  }
  \item{GLMSTEP.scope}{ range of models examined in the stepwise search; see also \code{\link[MASS]{stepAIC}}  }
  \item{GLMSTEP.k}{ multiple of the number of degrees of freedom used for the penalty (only k = 2 gives the genuine AIC); see also \code{\link[MASS]{stepAIC}}  }
  \item{GAM.formula}{ formula for the generalized additive model; see also \code{\link[gam]{gam}}  }
  \item{GAM.family}{ description of the error distribution and link function for the generalized additive model; see also \code{\link[gam]{gam}}  }
  \item{GAMSTEP.steps}{ maximum number of steps to be considered in the stepwise generalized additive model; see also \code{\link[gam]{step.gam}}  }
  \item{GAMSTEP.scope}{ range of models examined in the step-wise search n the stepwise generalized additive model; see also \code{\link[gam]{step.gam}}  }
  \item{MGCV.formula}{ formula for the generalized additive model; see also \code{\link[mgcv]{gam}}  }
  \item{EARTH.formula}{ formula for the multivariate adaptive regression spline model; see also \code{\link[earth]{earth}}  }
  \item{EARTH.glm}{ list of arguments to pass on to \code{\link[stats]{glm}}; see also \code{\link[earth]{earth}}  }
  \item{RPART.formula}{ formula for the recursive partioning and regression tree model; see also \code{\link[rpart]{rpart}}  }
  \item{RPART.xval}{ number of cross-validations for the recursive partioning and regression tree model; see also \code{\link[rpart]{rpart.control}}  }
  \item{NNET.formula}{ formula for the artificial neural network model; see also \code{\link[nnet]{nnet}}  }
  \item{NNET.size}{ number of units in the hidden layer for the artificial neural network model; see also \code{\link[nnet]{nnet}}  }
  \item{NNET.decay}{ parameter of weight decay for the artificial neural network model; see also \code{\link[nnet]{nnet}}  }
  \item{FDA.formula}{ formula for the flexible discriminant analysis model; see also \code{\link[mda]{fda}}  }
  \item{SVM.formula}{ formula for the support vector machine model; see also \code{\link[kernlab]{ksvm}}  }
  \item{digits}{ number of decimal place after multiplying results by 100 (percentages); see also \code{\link[base]{round}}  }
  \item{SCRIPT}{ provide suggested parameters (weights) for the \code{\link{ensemble.raster}} function (if \code{TRUE})  }
  \item{complexity}{ vector with values of complexity of individual trees (\code{tree.complexity}) for boosted regression trees; see also \code{\link[dismo]{gbm.step}}  }
  \item{learning}{ vector with values of weights applied to individual trees (\code{learning.rate}) for boosted regression trees; see also \code{\link[dismo]{gbm.step}}  }
  \item{size}{ vector with values of number of units in the hidden layer for the artificial neural network model; see also \code{\link[nnet]{nnet}}  }
  \item{decay}{ vector with values of weight decay for the artificial neural network model; see also \code{\link[nnet]{nnet}}  }
  \item{dummy.vars}{ vector that indicates which variables are dummy variables }
}

\details{
The basic function \code{ensemble.test} first calibrates individual suitability models based on presence locations \code{p} and background locations \code{a}, then evaluates these suitability models based on presence locations \code{pt} and background locations \code{at}. While calibrating and testing individual models, results obtained via the \code{\link[dismo]{evaluate}} function are shown in the GUI, and possibly plotted (\code{PLOTS}) or saved (\code{evaluations.keep}).

Function \code{ensemble.test.splits} splits the presence and background locations in a user-defined (\code{k}) number of subsets, then sequentially calibrates individual suitability models with (\code{k-1}) combined subsets and evaluates with the remaining one subset, whereby each subset is used once for evaluation in the user-defined number (\code{k}) of runs. For example, \code{k = 5} results in splitting the locations in 5 subsets, then using one of these subsets in turn for evaluations (see also \code{\link[dismo]{kfold}}).

Function \code{ensemble.test.gbm} allows to test various combinations of parameters \code{tree.complexity} and \code{learning.rate} for the \code{\link[dismo]{gbm.step}} model.

Function \code{ensemble.test.nnet} allows to test various combinations of parameters \code{size} and \code{decay} for the \code{\link[nnet]{nnet}} model.

Function \code{ensemble.formulae} provides suggestions for formulae that can be used for \code{ensemble.test} and \code{\link{ensemble.raster}}.
}

\value{
Function \code{ensemble.test} (potentially) returns a list with results from evaluations (via \code{\link[dismo]{evaluate}}) of calibration and test runs of individual suitability models.

Function \code{ensemble.test.splits} returns a matrix with, for each individual suitability model, the AUC of each run and the average AUC over the runs. Models are sorted by the average AUC. The average AUC for each model can be used as input weights for the \code{\link{ensemble.raster}} function (these weights are provided when the option of \code{SCRIPT} is selected.

Functions \code{ensemble.test.gbm} and \code{ensemble.test.nnet} return a matrix with, for each combination of model parameters, the AUC of each run and the average AUC. Models are sorted by the average AUC.
}

\references{
Buisson L, Thuiller W, Casajus N, Lek S and Grenouillet G. 2010. Uncertainty in ensemble forecasting of species distribution. Global Change Biology 16: 1145-1157
}



\author{ Roeland Kindt (World Agroforestry Centre) }

\seealso{
\code{\link{ensemble.raster}}
}

\examples{
\dontrun{
# based on examples in the dismo package
# simplified example of ensemble modelling with 4 modeling algorithms

# get predictor variables
library(dismo)
predictors <- stack(list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), pattern='grd', full.names=TRUE ))
# predictors

# presence points
presence <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
pres <- read.table(presence, header=TRUE, sep=',')[,-1]

# the kfold function randomly assigns data to groups; subdivide in calibration and training data
groupp <- kfold(pres, 5)
pres_train <- pres[groupp != 1, ]
pres_test <- pres[groupp == 1, ]

# choose background points
ext = extent(-90, -32, -33, 23)
background <- randomPoints(predictors, n=1000, ext=ext, extf = 1.00)
colnames(background) = c('lon', 'lat')
groupa <- kfold(background, 5)
backg_train <- background[groupa != 1, ]
backg_test <- background[groupa == 1, ]

# formulae for random forest and generalized linear model
# compare with: ensemble.formulae(predictors, factors=c("biome"))

rfformula <- as.formula(pb ~ bio1+bio5+bio6+bio7+bio8+bio12+bio16+bio17)

glmformula <- as.formula(pb ~ bio1 + I(bio1^2) + I(bio1^3) +
    bio5 + I(bio5^2) + I(bio5^3) + bio6 + I(bio6^2) + I(bio6^3) + 
    bio7 + I(bio7^2) + I(bio7^3) + bio8 + I(bio8^2) + I(bio8^3) +
    bio12 + I(bio12^2) + I(bio12^3) + bio16 + I(bio16^2) + I(bio16^3) + 
    bio17 + I(bio17^2) + I(bio17^3) )

# fit four ensemble models (RF, GLM, BIOCLIM, DOMAIN)
ensemble.nofactors <- ensemble.test(x=predictors, p=pres_train, a=backg_train, pt=pres_test, at=backg_test,
    MAXENT=0, GBM=0, GBMSTEP=0, RF=1, GLM=1, GLMSTEP=0, GAM=0, GAMSTEP=0, MGCV=0, 
    EARTH=0, RPART=0, NNET=0, FDA=0, SVM=0, BIOCLIM=1, DOMAIN=1, MAHAL=0,
    Yweights="BIOMOD", factors=NULL,
    PLOTS=FALSE, evaluations.keep=TRUE,
    RF.formula=rfformula,
    GLM.formula=glmformula)

# test performance of different suitability models; data are split in 4 subsets, each used once for evaluation
# GAMSTEP option is not available at the moment (problems with assignments)
ensemble.nofactors2 <- ensemble.test.splits(x=predictors, p=pres, a=background, k=4, 
    MAXENT=1, GBM=1, GBMSTEP=1, RF=1, GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=0, MGCV=1, 
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, BIOCLIM=1, DOMAIN=1, MAHAL=1,
    Yweights="BIOMOD", factors=NULL,
    PLOTS=FALSE, formulae.defaults=TRUE,
    GBMSTEP.learning.rate=0.002)
ensemble.nofactors2

}
}
