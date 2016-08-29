\name{ensemble.test}
\alias{ensemble.test}
\alias{ensemble.test.splits}
\alias{ensemble.test.gbm}
\alias{ensemble.test.nnet}
\alias{ensemble.drop1}
\alias{ensemble.formulae}
\alias{ensemble.weights}
\alias{ensemble.strategy}
\alias{ensemble.threshold}

\title{
Suitability mapping based on ensembles of modelling algorithms: comparison of different algorithms and calibration
}

\description{
The basic function \code{ensemble.test} allows to evaluate different algorithms for (species) suitability modelling, including maximum entropy (MAXENT), boosted regression trees, random forests, generalized linear models (including stepwise selection of explanatory variables), generalized additive models (including stepwise selection of explanatory variables), multivariate adaptive regression splines, regression trees, artificial neural networks, flexible discriminant analysis, support vector machines, the BIOCLIM algorithm, the DOMAIN algorithm and the Mahalanobis algorithm. These sets of functions were developed in parallel with the \code{biomod2} package, especially for inclusion of the maximum entropy algorithm, but also to allow for a more direct integration with the BiodiversityR package, more direct handling of model formulae and greater focus on mapping. Researchers and students of species distribution are strongly encouraged to familiarize themselves with all the options of the BIOMOD and dismo packages.
}


\usage{
ensemble.test(x = NULL, p = NULL, a = NULL, an = 1000, excludep = FALSE, ext = NULL, 
    k = 0, pt = NULL, at = NULL, CIRCLES.at = FALSE, CIRCLES.d = 100000, 
    TrainData = NULL, TestData = NULL, 
    VIF = FALSE, COR = FALSE,
    SINK = FALSE, PLOTS = TRUE, 
    threshold.method = "spec_sens", threshold.sensitivity = 0.9, 
    threshold.PresenceAbsence = FALSE,
    evaluations.keep = FALSE, 
    models.list = NULL, models.keep = FALSE, 
    models.save = FALSE, species.name = "Species001",
    AUC.weights = TRUE, ENSEMBLE.tune = FALSE, 
    ENSEMBLE.best = 0, ENSEMBLE.min = 0.7, ENSEMBLE.exponent = 1, 
    ENSEMBLE.weight.min = 0.05,
    input.weights = NULL, 
    MAXENT = 1, GBM = 1, GBMSTEP = 1, RF = 1, GLM = 1, GLMSTEP = 1, 
    GAM = 1, GAMSTEP = 1, MGCV = 1, MGCVFIX = 0, EARTH = 1, 
    RPART = 1, NNET = 1, FDA = 1, SVM = 1 , SVME = 1,
    BIOCLIM = 1, DOMAIN = 1, MAHAL = 1, 
    GEODIST = 0, 
    PROBIT = FALSE,
    Yweights = "BIOMOD", 
    layer.drops = NULL, factors = NULL, dummy.vars = NULL,
    formulae.defaults = TRUE, maxit = 100,
    MAXENT.a = NULL, MAXENT.an = 10000, MAXENT.BackData = NULL, 
    MAXENT.path=paste(getwd(), "/models/maxent_", species.name,  sep=""), 
    GBM.formula = NULL, GBM.n.trees = 2001, 
    GBMSTEP.gbm.x = 2:(ncol(TrainData.vars)+1), GBMSTEP.tree.complexity = 5, 
    GBMSTEP.learning.rate = 0.005, GBMSTEP.bag.fraction = 0.5, 
    GBMSTEP.step.size = 100, 
    RF.formula = NULL, RF.ntree = 751, 
    RF.mtry = floor(sqrt(ncol(TrainData.vars))), 
    GLM.formula = NULL, GLM.family = binomial(link = "logit"), 
    GLMSTEP.steps = 1000, STEP.formula = NULL, GLMSTEP.scope = NULL, 
    GLMSTEP.k = 2, 
    GAM.formula = NULL, GAM.family = binomial(link = "logit"), 
    GAMSTEP.steps = 1000, GAMSTEP.scope = NULL, GAMSTEP.pos = 1, 
    MGCV.formula = NULL, MGCV.select = FALSE, 
    MGCVFIX.formula = NULL, 
    EARTH.formula = NULL, 
    EARTH.glm = list(family = binomial(link = "logit"), maxit = maxit), 
    RPART.formula = NULL, RPART.xval = 50, 
    NNET.formula = NULL, NNET.size = 8, NNET.decay = 0.01, 
    FDA.formula = NULL, 
    SVM.formula = NULL, 
    SVME.formula = NULL, 
    MAHAL.shape = 1, 
    RASTER.format = "raster")

ensemble.test.splits(x = NULL, p = NULL, a = NULL, an = 1000, 
    CIRCLES.at = FALSE, CIRCLES.d = 100000, 
    excludep = FALSE, ext = NULL, 
    k = 4, 
    TrainData = NULL,
    VIF = FALSE, COR = FALSE,
    SINK = FALSE, PLOTS = FALSE, 
    data.keep = FALSE,
    species.name = "Species001",
    threshold.method = "spec_sens", threshold.sensitivity = 0.9, 
    threshold.PresenceAbsence = FALSE,
    AUC.weights = TRUE, ENSEMBLE.tune = FALSE, 
    ENSEMBLE.best = 0, ENSEMBLE.min = 0.7, ENSEMBLE.exponent = 1,
    ENSEMBLE.weight.min = 0.05,
    input.weights = NULL,
    MAXENT = 1, GBM = 1, GBMSTEP = 1, RF = 1, GLM = 1, GLMSTEP = 1, 
    GAM = 1, GAMSTEP = 1, MGCV = 1, MGCVFIX = 0, EARTH = 1, 
    RPART = 1, NNET = 1, FDA = 1, SVM = 1 , SVME = 1, 
    BIOCLIM = 1, DOMAIN = 1, MAHAL = 1,
    PROBIT = FALSE,
    Yweights = "BIOMOD", 
    layer.drops = NULL, factors = NULL, dummy.vars = NULL,
    formulae.defaults = TRUE, maxit = 100,
    MAXENT.a = NULL, MAXENT.an = 10000, MAXENT.BackData = NULL, 
    MAXENT.path = paste(getwd(), "/models/maxent_", species.name,  sep=""), 
    GBM.formula = NULL, GBM.n.trees = 2001, 
    GBMSTEP.gbm.x = 2:(ncol(TrainData1)), GBMSTEP.tree.complexity = 5, 
    GBMSTEP.learning.rate = 0.005, 
    GBMSTEP.bag.fraction = 0.5, GBMSTEP.step.size = 100, 
    RF.formula = NULL, RF.ntree = 751, RF.mtry = floor(sqrt(ncol(TrainData1)-1)), 
    GLM.formula = NULL, GLM.family = binomial(link = "logit"), 
    GLMSTEP.steps = 1000, STEP.formula = NULL, GLMSTEP.scope = NULL, GLMSTEP.k = 2, 
    GAM.formula = NULL, GAM.family = binomial(link = "logit"), 
    GAMSTEP.steps = 1000, GAMSTEP.scope = NULL, GAMSTEP.pos = 1, 
    MGCV.formula = NULL, MGCV.select = FALSE, 
    MGCVFIX.formula = NULL, 
    EARTH.formula = NULL, 
    EARTH.glm = list(family = binomial(link = "logit"), maxit = maxit), 
    RPART.formula = NULL, RPART.xval = 50, 
    NNET.formula = NULL, NNET.size = 8, NNET.decay = 0.01, 
    FDA.formula = NULL, 
    SVM.formula = NULL, 
    SVME.formula = NULL, 
    MAHAL.shape = 1)

ensemble.test.gbm(x = NULL, p = NULL, a = NULL, an = 1000, excludep = FALSE, ext = NULL, 
    k = 4, 
    TrainData = NULL,
    VIF = FALSE, COR = FALSE,
    SINK = FALSE, PLOTS = FALSE, 
    species.name = "Species001",
    Yweights = "BIOMOD", 
    layer.drops = NULL, factors = NULL, 
    GBMSTEP.gbm.x = 2:(ncol(TrainData.orig)), 
    complexity = c(3:6), learning = c(0.005, 0.002, 0.001), 
    GBMSTEP.bag.fraction = 0.5, GBMSTEP.step.size = 100)

ensemble.test.nnet(x = NULL, p = NULL, a = NULL, an = 1000, excludep = FALSE, 
    ext = NULL, k = 4, 
    TrainData = NULL,
    VIF = FALSE, COR = FALSE,
    SINK = FALSE, PLOTS = FALSE, 
    species.name = "Species001",
    Yweights = "BIOMOD", 
    layer.drops = NULL, factors = NULL, 
    formulae.defaults = TRUE, maxit = 100, 
    NNET.formula = NULL, 
    sizes = c(2, 4, 6, 8), decays = c(0.1, 0.05, 0.01, 0.001) )

ensemble.drop1(x = NULL, p = NULL, a = NULL, an = 1000, excludep = FALSE, ext = NULL, 
    k = 0, pt = NULL, at = NULL, CIRCLES.at = FALSE, CIRCLES.d = 100000, 
    TrainData = NULL, TestData = NULL,
    VIF = FALSE, COR = FALSE,
    SINK = FALSE, 
    species.name = "Species001",
    difference = FALSE, 
    ENSEMBLE.best = 0, ENSEMBLE.min = 0.7, ENSEMBLE.exponent = 1,
    input.weights = NULL, 
    MAXENT = 1, GBM = 1, GBMSTEP = 1, RF = 1, GLM = 1, GLMSTEP = 1, 
    GAM = 1, GAMSTEP = 1, MGCV = 1, MGCVFIX = 0, EARTH = 1, 
    RPART = 1, NNET = 1, FDA = 1, SVM = 1, SVME = 1, 
    BIOCLIM = 1, DOMAIN = 1, MAHAL = 1, 
    PROBIT = FALSE,
    Yweights = "BIOMOD", 
    layer.drops = NULL, factors = NULL, dummy.vars = NULL, 
    maxit = 100,
    MAXENT.a = NULL, MAXENT.an = 10000, MAXENT.BackData = NULL, 
    MAXENT.path = paste(getwd(), "/models/maxent_", species.name,  sep=""), 
    GBM.n.trees = 2001, 
    GBMSTEP.tree.complexity = 5, GBMSTEP.learning.rate = 0.005, 
    GBMSTEP.bag.fraction = 0.5, GBMSTEP.step.size = 100, 
    RF.ntree = 751, 
    GLM.family = binomial(link = "logit"), 
    GLMSTEP.steps = 1000, GLMSTEP.scope = NULL, GLMSTEP.k = 2, 
    GAM.family = binomial(link = "logit"), 
    GAMSTEP.steps = 1000, GAMSTEP.scope = NULL, GAMSTEP.pos = 1, 
    MGCV.select = FALSE, 
    EARTH.glm = list(family = binomial(link = "logit"), maxit = maxit), 
    RPART.xval = 50, 
    NNET.size = 8, NNET.decay = 0.01, 
    MAHAL.shape = 1)

ensemble.weights(weights = c(0.9, 0.8, 0.7, 0.5), 
    best = 0, min.weight = 0, 
    exponent = 1, digits = 4)

ensemble.strategy(TrainData = NULL, TestData = NULL, 
    verbose = FALSE,
    ENSEMBLE.best = c(4:10), ENSEMBLE.min = c(0.7),
    ENSEMBLE.exponent = c(1, 2, 4, 6, 8) )

ensemble.formulae(x, factors = NULL, dummy.vars = NULL)

ensemble.threshold(eval, threshold.method="spec_sens", threshold.sensitivity=0.9, 
    threshold.PresenceAbsence=FALSE, Pres, Abs) 

}

\arguments{
  \item{x}{ RasterStack object (\code{\link[raster]{stack}}) containing all layers that correspond to explanatory variables  }
  \item{p}{ presence points used for calibrating the suitability models, typically available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{a}{ background points used for calibrating the suitability models (except for \code{\link[dismo]{maxent}}), typically available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{an}{ number of background points for calibration to be selected with \code{\link[dismo]{randomPoints}} in case argument \code{a} is missing}
  \item{excludep}{ parameter that indicates (if \code{TRUE}) that presence points will be excluded from the background points; see also \code{\link[dismo]{randomPoints}} }
  \item{ext}{ an Extent object to limit the selection of background points to a sub-region of \code{x}, typically provided as c(lonmin, lonmax, latmin, latmax); see also \code{\link[dismo]{randomPoints}} and \code{\link[raster]{extent}} }
  \item{k}{ If larger than 1, the number of groups to split between calibration (k-1) and evaluation (1) data sets (for example, \code{k = 4} results in 3/4 of presence and background points to be used for calibrating the models, and 1/4 of presence and background points to be used for evaluating the models). For \code{ensemble.test.splits}, \code{ensemble.test.gbm} and \code{ensemble.test.nnet}, this procedure is repeated \code{k} times (k-fold cross-validation). See also \code{\link[dismo]{kfold}}.  }
  \item{pt}{ presence points used for evaluating the suitability models, available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{at}{ background points used for evaluating the suitability models, available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{CIRCLES.at}{ If \code{TRUE}, then new background points that will be used for evaluationg the suitability models will be selected (\code{\link[dismo]{randomPoints}}) in circular neighbourhoods (created with \code{\link[dismo]{circles}}) around presence locations (p and pt).  }
  \item{CIRCLES.d}{ Radius in m of circular neighbourhoods (created with \code{\link[dismo]{circles}}) around presence locations (p and pt).  }
  \item{TrainData}{ dataframe with first column 'pb' describing presence (1) and absence (0) and other columns containing explanatory variables; see also \code{\link[dismo]{prepareData}}. In case that this dataframe is provided, then locations p and a are not used. For the maximum entropy model (\code{\link[dismo]{maxent}}), a different dataframe is used for calibration; see parameter \code{MAXENT.TrainData}. }
  \item{TestData}{ dataframe with first column 'pb' describing presence (1) and absence (0) and other columns containing explanatory variables; see also \code{\link[dismo]{prepareData}}. In case that this dataframe is provided, then locations pt and at are not used. For \code{ensemble.strategy}, this dataframe should be a dataframe that contains predictions for various models - such dataframe can be provided by the \code{ensemble.test} or \code{\link{ensemble.raster}} functions.  }
  \item{VIF}{ Estimate the variance inflation factors based on a linear model calibrated on the training data  (if \code{TRUE}). Only background locations will be used and the response variable 'pb' will be replaced by a random variable. See also \code{\link[car]{vif}}.  }
  \item{COR}{ Provide information on the correlation between the numeric explanatory variables  (if \code{TRUE}). See also \code{\link[stats]{cor}}.  }
  \item{SINK}{ Append the results to a text file in subfolder 'outputs' (if \code{TRUE}). The name of file is based on argument \code{species.name}. In case the file already exists, then results are appended. See also \code{\link[base]{sink}}.  }
  \item{PLOTS}{ Plot the results of evaluation for the various suitability models (if \code{TRUE}). See also \code{\link[dismo]{evaluate}}.  }
  \item{threshold.method}{ Method to calculate the threshold between predicted absence and presence; possibilities include \code{spec_sens} (highest sum of the true positive rate and the true negative rate), \code{kappa} (highest kappa value), \code{no_omission} (highest threshold that corresponds to no omission), \code{prevalence} (modeled prevalence is closest to observed prevalence) and \code{equal_sens_spec} (equal true positive rate and true negative rate). See \code{\link[dismo]{threshold}}. Options specific to the BiodiversityR implementation are: \code{threshold.mean} (resulting in calculating the mean value of \code{spec_sens}, \code{equal_sens_spec} and \code{prevalence}) and \code{threshold.min} (resulting in calculating the minimum value of \code{spec_sens}, \code{equal_sens_spec} and \code{prevalence}).   }
  \item{threshold.sensitivity}{ Sensitivity value for \code{threshold.method = 'sensitivity'}. See \code{\link[dismo]{threshold}}.  }
  \item{threshold.PresenceAbsence}{ If \code{TRUE} calculate thresholds with the PresenceAbsence package. See \code{\link[PresenceAbsence]{optimal.thresholds}}.  }
  \item{evaluations.keep}{ Keep the results of evaluations (if \code{TRUE}). See also \code{\link[dismo]{evaluate}}. }
  \item{models.list}{ list with 'old' model objects such as \code{MAXENT} or \code{RF}. }
  \item{models.keep}{ store the details for each suitability modelling algorithm (if \code{TRUE}). (This may be particularly useful when projecting to different possible future climates.)  }
  \item{models.save}{ Save the list with model details to a file (if \code{TRUE}). The filename will be \code{species.name} with extension \code{.models}; this file will be saved in subfolder of \code{models}. When loading this file, model results will be available as \code{ensemble.models}. }
  \item{species.name}{ Name by which the model details will be saved to a file; see also argument \code{models.save} }
  \item{data.keep}{ Keep the data for each k-fold cross-validation run (if \code{TRUE}). }
  \item{AUC.weights}{ If \code{TRUE}, then AUC values are used as weights for the ensemble model, scaled to sum to 1 for the ensemble model through \code{ensemble.weights}. If \code{FALSE}, then input weights are scaled to sum to 1 for the ensemble model, but AUC values are not considered. }
  \item{ENSEMBLE.tune}{ Determine weights for the ensemble model based on AUC values through internal cross-validation procedures (if \code{TRUE}). See details. }
  \item{ENSEMBLE.best}{ The number of individual suitability models to be used in the consensus suitability map (based on a weighted average). In case this parameter is smaller than 1 or larger than the number of positive input weights of individual models, then all individual suitability models with positive input weights are included in the consensus suitability map. In case a vector is provided, \code{ensemble.strategy} is called internally to determine weights for the ensemble model. }
  \item{ENSEMBLE.min}{ The minimum input weight (typically corresponding to AUC values) for a model to be included in the ensemble. In case a vector is provided, function \code{ensemble.strategy} is called internally to determine weights for the ensemble model.  }
  \item{ENSEMBLE.exponent}{ Exponent applied to AUC values to convert AUC values into weights (for example, an exponent of 2 converts input weights of 0.7, 0.8 and 0.9 into 0.7^2=0.49, 0.8^2=0.64 and 0.9^2=0.81). See details. }
  \item{ENSEMBLE.weight.min}{ The minimum output weight for models included in the ensemble, applying to weights that sum to one. Note that \code{ENSEMBLE.min} typically refers to input AUC values.   }
  \item{input.weights}{ array with numeric values for the different modelling algorithms; if \code{NULL} then values provided by parameters such as \code{MAXENT} and \code{GBM} will be used. As an alternative, the output from \code{ensemble.test.splits} can be used. }
  \item{MAXENT}{ number: if larger than 0, then a maximum entropy model (\code{\link[dismo]{maxent}}) will be fitted among ensemble }
  \item{GBM}{ number: if larger than 0, then a boosted regression trees model (\code{\link[gbm]{gbm}}) will be fitted among ensemble }
  \item{GBMSTEP}{ number: if larger than 0, then a stepwise boosted regression trees model (\code{\link[dismo]{gbm.step}}) will be fitted among ensemble }
  \item{RF}{ number: if larger than 0, then a random forest model (\code{\link[randomForest]{randomForest}}) will be fitted among ensemble }
  \item{GLM}{ number: if larger than 0, then a generalized linear model (\code{\link[stats]{glm}}) will be fitted among ensemble }
  \item{GLMSTEP}{ number: if larger than 0, then a stepwise generalized linear model (\code{\link[MASS]{stepAIC}}) will be fitted among ensemble }
  \item{GAM}{ number: if larger than 0, then a generalized additive model (\code{\link[gam]{gam}}) will be fitted among ensemble }
  \item{GAMSTEP}{ number: if larger than 0, then a stepwise generalized additive model (\code{\link[gam]{step.gam}}) will be fitted among ensemble }
  \item{MGCV}{ number: if larger than 0, then a generalized additive model (\code{\link[mgcv]{gam}}) will be fitted among ensemble }
  \item{MGCVFIX}{ number: if larger than 0, then a generalized additive model with fixed d.f. regression splines (\code{\link[mgcv]{gam}}) will be fitted among ensemble }
  \item{EARTH}{ number: if larger than 0, then a multivariate adaptive regression spline model (\code{\link[earth]{earth}}) will be fitted among ensemble }
  \item{RPART}{ number: if larger than 0, then a recursive partioning and regression tree model (\code{\link[rpart]{rpart}}) will be fitted among ensemble }
  \item{NNET}{ number: if larger than 0, then an artificial neural network model (\code{\link[nnet]{nnet}}) will be fitted among ensemble }
  \item{FDA}{ number: if larger than 0, then a flexible discriminant analysis model (\code{\link[mda]{fda}}) will be fitted among ensemble }
  \item{SVM}{ number: if larger than 0, then a support vector machine model (\code{\link[kernlab]{ksvm}}) will be fitted among ensemble }
  \item{SVME}{ number: if larger than 0, then a support vector machine model (\code{\link[e1071]{svm}}) will be fitted among ensemble }
  \item{BIOCLIM}{ number: if larger than 0, then the BIOCLIM algorithm (\code{\link[dismo]{bioclim}}) will be fitted among ensemble }
  \item{DOMAIN}{ number: if larger than 0, then the DOMAIN algorithm (\code{\link[dismo]{domain}}) will be fitted among ensemble }
  \item{MAHAL}{ number: if larger than 0, then the Mahalanobis algorithm (\code{\link[dismo]{mahal}}) will be fitted among ensemble }
  \item{PROBIT}{ If \code{TRUE}, then subsequently to the fitting of the individual algorithm (e.g. maximum entropy or GAM) a generalized linear model (\code{\link[stats]{glm}}) with probit link \code{family=binomial(link="probit")} will be fitted to transform the predictions, using the previous predictions as explanatory variable. This transformation results in all model predictions to be probability estimates.  }
  \item{GEODIST}{ number: if larger than 0, then the geoDist algorithm (\code{\link[dismo]{geoDist}}) will be fitted among ensemble (note that this algorithm does not use environmental layers, and is based only on the distance from presence points used to calibrate this algorithm) }
  \item{Yweights}{ chooses how cases of presence and background (absence) are weighted; \code{"BIOMOD"} results in equal weighting of all presence and all background cases, \code{"equal"} results in equal weighting of all cases. The user can supply a vector of weights similar to the number of cases in the calibration data set. }
  \item{layer.drops}{ vector that indicates which layers should be removed from RasterStack \code{x}. This argument is especially useful for the \code{ensemble.drop1} function. See also \code{\link[raster]{addLayer}}.  }
  \item{factors}{ vector that indicates which variables are factors; see also \code{\link[dismo]{prepareData}}  }
  \item{dummy.vars}{ vector that indicates which variables are dummy variables (influences formulae suggestions) }
  \item{formulae.defaults}{ Suggest formulae for most of the models (if \code{TRUE}). See also \code{\link{ensemble.formulae}}.   }
  \item{maxit}{ Maximum number of iterations for some of the models. See also \code{\link[stats]{glm.control}}, \code{\link[gam]{gam.control}}, \code{\link[mgcv]{gam.control}} and \code{\link[nnet]{nnet}}.  }
  \item{MAXENT.a}{ background points used for calibrating the maximum entropy model (\code{\link[dismo]{maxent}}), typically available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}. Ignored if \code{MAXENT.BackData} is provided.  }
  \item{MAXENT.an}{ number of background points for calibration to be selected with \code{\link[dismo]{randomPoints}} in case argument \code{MAXENT.a} is missing}
  \item{MAXENT.BackData}{ dataframe containing explanatory variables for the background locations. This information will be used for calibrating the maximum entropy model (\code{\link[dismo]{maxent}}). When used with the \code{ensemble.test.splits} function, the same background locations will be used for each of the cross-validation runs; this is based on the assumption that a large number (~10000) of background locations are used.  }
  \item{MAXENT.path}{ path to the directory where output files of the maximum entropy model are stored; see also \code{\link[dismo]{maxent}}   }
  \item{GBM.formula}{ formula for the boosted regression trees algorithm; see also \code{\link[gbm]{gbm}}  }
  \item{GBM.n.trees}{ total number of trees to fit for the boosted regression trees model; see also \code{\link[gbm]{gbm}}  }
  \item{GBMSTEP.gbm.x}{ indices of column numbers with explanatory variables for stepwise boosted regression trees; see also \code{\link[dismo]{gbm.step}}  }
  \item{GBMSTEP.tree.complexity}{ complexity of individual trees for stepwise boosted regression trees; see also \code{\link[dismo]{gbm.step}}  }
  \item{GBMSTEP.learning.rate}{ weight applied to individual trees for stepwise boosted regression trees; see also \code{\link[dismo]{gbm.step}} }
  \item{GBMSTEP.bag.fraction}{ proportion of observations used in selecting variables for stepwise boosted regression trees; see also \code{\link[dismo]{gbm.step}}  }
  \item{GBMSTEP.step.size}{ number of trees to add at each cycle for stepwise boosted regression trees (should be small enough to result in a smaller holdout deviance than the initial number of trees [50]); see also \code{\link[dismo]{gbm.step}}  }
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
  \item{GAMSTEP.pos}{ parameter expected to be set to 1 to allow for fitting of the stepwise generalized additive model  }
  \item{MGCV.formula}{ formula for the generalized additive model; see also \code{\link[mgcv]{gam}}  }
  \item{MGCV.select}{ if \code{TRUE}, then the smoothing parameter estimation that is part of fitting can completely remove terms from the model; see also \code{\link[mgcv]{gam}}  }
  \item{MGCVFIX.formula}{ formula for the generalized additive model with fixed d.f. regression splines; see also \code{\link[mgcv]{gam}} (the default formulae sets "s(..., fx = TRUE, ...)"; see also \code{\link[mgcv]{s}}) }
  \item{EARTH.formula}{ formula for the multivariate adaptive regression spline model; see also \code{\link[earth]{earth}}  }
  \item{EARTH.glm}{ list of arguments to pass on to \code{\link[stats]{glm}}; see also \code{\link[earth]{earth}}  }
  \item{RPART.formula}{ formula for the recursive partioning and regression tree model; see also \code{\link[rpart]{rpart}}  }
  \item{RPART.xval}{ number of cross-validations for the recursive partioning and regression tree model; see also \code{\link[rpart]{rpart.control}}  }
  \item{NNET.formula}{ formula for the artificial neural network model; see also \code{\link[nnet]{nnet}}  }
  \item{NNET.size}{ number of units in the hidden layer for the artificial neural network model; see also \code{\link[nnet]{nnet}}  }
  \item{NNET.decay}{ parameter of weight decay for the artificial neural network model; see also \code{\link[nnet]{nnet}}  }
  \item{FDA.formula}{ formula for the flexible discriminant analysis model; see also \code{\link[mda]{fda}}  }
  \item{SVM.formula}{ formula for the support vector machine model; see also \code{\link[kernlab]{ksvm}}  }
  \item{SVME.formula}{ formula for the support vector machine model; see also \code{\link[e1071]{svm}}  }
  \item{MAHAL.shape}{ parameter that influences the transformation of output values of \code{\link[dismo]{mahal}}. See details section.  }
  \item{RASTER.format}{ Format of the raster files that will be generated for the GEODIST model. See \code{\link[raster]{writeFormats}} and \code{\link[raster]{writeRaster}}.  }

  \item{complexity}{ vector with values of complexity of individual trees (\code{tree.complexity}) for boosted regression trees; see also \code{\link[dismo]{gbm.step}}  }
  \item{learning}{ vector with values of weights applied to individual trees (\code{learning.rate}) for boosted regression trees; see also \code{\link[dismo]{gbm.step}}  }

  \item{sizes}{ vector with values of number of units in the hidden layer for the artificial neural network model; see also \code{\link[nnet]{nnet}}  }
  \item{decays}{ vector with values of weight decay for the artificial neural network model; see also \code{\link[nnet]{nnet}}  }

  \item{difference}{ if \code{TRUE}, then AUC values of the models with all variables are subtracted from the models where one explanatory variable was excluded. After subtraction, positive values indicate that the model without the explanatory variable has a higher AUC than the model with all variables.   }

  \item{weights}{ input weights for the \code{ensemble.weights} function }
  \item{best}{ The number of final weights. In case this parameter is smaller than 1 or larger than the number of positive input weights of individual models, then this parameter is ignored. }
  \item{min.weight}{ The minimum input weight to be included in the output.  }
  \item{exponent}{ Exponent applied to AUC values to convert AUC values into weights (for example, an exponent of 2 converts input weights of 0.7, 0.8 and 0.9 into 0.7^2=0.49, 0.8^2=0.64 and 0.9^2=0.81). See details. }
  \item{digits}{ Number of number of decimal places in the output weights; see also \code{\link[base]{round}}.  }

  \item{verbose}{ If \code{TRUE}, then provide intermediate results for \code{ensemble.strategy}) }

  \item{eval}{ ModelEvaluation object obtained by \code{\link[dismo]{evaluate}} }
  \item{Pres}{ Suitabilities (probabilities) at presence locations }
  \item{Abs}{ Suitabilities (probabilities) at background locations }

}

\details{
The basic function \code{ensemble.test} first calibrates individual suitability models based on presence locations \code{p} and background locations \code{a}, then evaluates these suitability models based on presence locations \code{pt} and background locations \code{at}. While calibrating and testing individual models, results obtained via the \code{\link[dismo]{evaluate}} function are shown in the GUI, and possibly plotted (\code{PLOTS}) or saved (\code{evaluations.keep}).

As an alternative to providing presence locations \code{p}, models can be calibrated with data provided in \code{TrainData}. In case that both \code{p} and \code{TrainData} are provided, then models will be calibrated with \code{TrainData}.

Calibration of the maximum entropy (MAXENT) algorithm is not based on background locations \code{a}, but based on background locations \code{MAXENT.a} instead. However, to compare evaluations with evaluations of other algorithms, during evaluations of the MAXENT algorithm, presence locations \code{p} and background locations \code{a} are used (and not background locations \code{MAXENT.a}). 

As the Mahalanobis function (\code{\link[dismo]{mahal}}) does not always provide values within a range of 0 - 1, the output values are rescaled by first subtracting the value of \code{1 - MAHAL.shape} from each prediction, followed by calculating the absolute value, followed by calculating the reciprocal value and finally multiplying this reciprocal value with \code{MAHAL.shape}. As this rescaling method does not estimate probabilities, inclusion in the calculation of a (weighted) average of ensemble probabilities may be problematic (the same applies to other distance-based methods). 

With parameter \code{ENSEMBLE.best}, the subset of best models (evaluated by the individual AUC values) can be selected and only those models will be used for calculating the ensemble model (in other words, weights for models not included in the ensemble will be set to zero). It is possible to further increase the contribution to the ensemble model for models with higher AUC values through parameter \code{ENSEMBLE.exponent}. With \code{ENSEMBLE.exponent = 2}, AUC values of 0.7, 0.8 and 0.9 are converted into weights of 0.7^2=0.49, 0.8^2=0.64 and 0.9^2=0.81). With \code{ENSEMBLE.exponent = 4}, AUC values of 0.7, 0.8 and 0.9 are converted into weights of 0.7^4=0.2401, 0.8^4=0.4096 and 0.9^2=0.6561). 

\code{ENSEMBLE.tune} will result in an internal procedure whereby the best selection of parameter values for \code{ENSEMBLE.min}, \code{ENSEMBLE.best} or \code{ENSEMBLE.exponent} can be identified. Through a factorial procedure, the ensemble model with best AUC for a specific combination of parameter values is identified. The procedure also provides the weights that correspond to the best ensemble.  

Function \code{ensemble.test.splits} splits the presence and background locations in a user-defined (\code{k}) number of subsets (i.e. k-fold cross-validation), then sequentially calibrates individual suitability models with (\code{k-1}) combined subsets and evaluates those with the remaining one subset, whereby each subset is used once for evaluation in the user-defined number (\code{k}) of runs. For example, \code{k = 4} results in splitting the locations in 4 subsets, then using one of these subsets in turn for evaluations (see also \code{\link[dismo]{kfold}}). Note that for the maximum entropy (MAXENT) algorithm, the same background data will be used in each cross-validation run (this is based on the assumption that a large number (~10000) of background locations are used).

Among the output from function \code{ensemble.test.splits} are suggested weights for an ensemble model (\code{output.weights} and \code{output.weights.AUC}), and information on the respective AUC values of the ensemble model with the suggested weights for each of the (\code{k}) subsets. Suggested weights \code{output.weights} are calculated as the average of the weights of the different algorithms (submodels) of the \code{k} ensembles. Suggested weights \code{output.weights.AUC} are calculated as the average of the AUC of the different algorithms of the for the \code{k} runs. 

Function \code{ensemble.test.gbm} allows to test various combinations of parameters \code{tree.complexity} and \code{learning.rate} for the \code{\link[dismo]{gbm.step}} model.

Function \code{ensemble.test.nnet} allows to test various combinations of parameters \code{size} and \code{decay} for the \code{\link[nnet]{nnet}} model.

Function \code{ensemble.drop1} allows to test the effects of leaving out each of the explanatory variables, and comparing these results with the "full" model. Note that option of \code{difference = TRUE} may result in positive values, indicating that the model without the explanatory variable having larger AUC than the "full" model. A procedure is included to estimate the deviance of a model based on the fitted values, using -2 * (sum(x*log(x)) + sum((1-x)*log(1-x))) where x is a vector of the fitted values for a respective model. (It was checked that this procedure results in similar deviance estimates for the null and 'full' models for glm, but note that it is not certain whether deviance can be calculated in a similar way for other submodels.)

Function \code{ensemble.formulae} provides suggestions for formulae that can be used for \code{ensemble.test} and \code{\link{ensemble.raster}}. This function is always used internally by the \code{ensemble.drop1} function.

The \code{ensemble.weights} function is used internally by the \code{ensemble.test} and \code{\link{ensemble.raster}} functions, using the input weights for the different suitability modelling algorithms. Ties between input weights result in the same output weights.

The \code{ensemble.strategy} function is used internally by the \code{ensemble.test} function, using the train and test data sets with predictions of the different suitability modelling algorithms and different combinations of parameters \code{ENSEMBLE.best}, \code{ENSEMBLE.min} and \code{ENSEMBLE.exponent}. The final ensemble model is based on the parameters that generate the best AUC.

The \code{ensemble.threshold} function is used internally by the \code{ensemble.test}, \code{ensemble.mean} and \code{ensemble.plot} functions. \code{threshold.mean} and \code{threshold.min} result in calculating the mean or minimum value of threshold methods that resulted in better results in a study by Liu et al. (Ecography 28: 385-393. 2005) with threshold availabe in \code{\link[dismo]{threshold}} (\code{prevalence}, \code{spec_sens} and \code{equal_spec_sens}) or \code{\link[PresenceAbsence]{optimal.thresholds}} (\code{ObsPrev}, \code{MeanProb}, \code{MaxSens+Spec}, \code{Sens=Spec} and \code{MinROCdist}).

}

\value{
Function \code{ensemble.test} (potentially) returns a list with results from evaluations (via \code{\link[dismo]{evaluate}}) of calibration and test runs of individual suitability models.

Function \code{ensemble.test.splits} returns a matrix with, for each individual suitability model, the AUC of each run and the average AUC over the runs. Models are sorted by the average AUC. The average AUC for each model can be used as input weights for the \code{\link{ensemble.raster}} function.

Functions \code{ensemble.test.gbm} and \code{ensemble.test.nnet} return a matrix with, for each combination of model parameters, the AUC of each run and the average AUC. Models are sorted by the average AUC.
}

\references{
Buisson L, Thuiller W, Casajus N, Lek S and Grenouillet G. 2010. Uncertainty in ensemble forecasting of species distribution. Global Change Biology 16: 1145-1157

Liu C, Berry PM, Dawson TP and Pearson RC. 2005. Selecting thresholds of occurrence in the prediction of species distributions. Ecography 28: 385-393
}



\author{ Roeland Kindt (World Agroforestry Centre) }

\seealso{
\code{\link{ensemble.raster}}
}

\examples{
\dontrun{
# based on examples in the dismo package

# get predictor variables
library(dismo)
predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
    pattern='grd', full.names=TRUE)
predictors <- stack(predictor.files)
# subset based on Variance Inflation Factors
predictors <- subset(predictors, subset=c("bio5", "bio6", 
    "bio16", "bio17", "biome"))
predictors
predictors@title <- "base"

# presence points
presence_file <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
pres <- read.table(presence_file, header=TRUE, sep=',')[,-1]

# the kfold function randomly assigns data to groups; 
# groups are used as calibration (1/5) and training (4/5) data
groupp <- kfold(pres, 5)
pres_train <- pres[groupp !=  1, ]
pres_test <- pres[groupp ==  1, ]

# choose background points
ext <- extent(-90, -32, -33, 23)
background <- randomPoints(predictors, n=1000, ext=ext, extf=1.00)
colnames(background)=c('lon', 'lat')
groupa <- kfold(background, 5)
backg_train <- background[groupa != 1, ]
backg_test <- background[groupa == 1, ]

# formulae for random forest and generalized linear model
# compare with: ensemble.formulae(predictors, factors=c("biome"))

rfformula <- as.formula(pb ~ bio5+bio6+bio16+bio17)

glmformula <- as.formula(pb ~ bio5 + I(bio5^2) + I(bio5^3) + 
    bio6 + I(bio6^2) + I(bio6^3) + bio16 + I(bio16^2) + I(bio16^3) + 
    bio17 + I(bio17^2) + I(bio17^3) )

# fit four ensemble models (RF, GLM, BIOCLIM, DOMAIN)
ensemble.nofactors <- ensemble.test(x=predictors, p=pres_train, a=backg_train, 
    pt=pres_test, at=backg_test,
    species.name="Bradypus",
    MAXENT=0, GBM=0, GBMSTEP=0, RF=1, GLM=1, GLMSTEP=0, GAM=0, 
    GAMSTEP=0, MGCV=0, MGCVFIX=0,EARTH=0, RPART=0, NNET=0, FDA=0, 
    SVM=0, SVME=0, BIOCLIM=1, DOMAIN=1, MAHAL=0,
    Yweights="BIOMOD", factors="biome",
    PLOTS=FALSE, evaluations.keep=TRUE,
    RF.formula=rfformula,
    GLM.formula=glmformula)

# fit four ensemble models (RF, GLM, BIOCLIM, DOMAIN) using default formulae
# variable 'biome' is not included as explanatory variable
# results are provided in a file in the 'outputs' subfolder of the working
# directory
ensemble.nofactors <- ensemble.test(x=predictors,
    p=pres_train, a=backg_train, 
    pt=pres_test, at=backg_test,
    layer.drops="biome",
    species.name="Bradypus",
    SINK=TRUE,
    MAXENT=0, GBM=0, GBMSTEP=0, RF=1, GLM=1, GLMSTEP=0, GAM=0, 
    GAMSTEP=0, MGCV=0, MGCVFIX=0, EARTH=0, RPART=0, NNET=0, FDA=0, 
    SVM=0, SVME=0, BIOCLIM=1, DOMAIN=1, MAHAL=0,
    Yweights="BIOMOD", factors="biome",
    PLOTS=FALSE, evaluations.keep=TRUE,
    formulae.defaults=TRUE)    

# after fitting the individual algorithms (submodels),
# transform predictions with a probit link.
ensemble.nofactors <- ensemble.test(x=predictors,
    p=pres_train, a=backg_train, 
    pt=pres_test, at=backg_test,
    layer.drops="biome",
    species.name="Bradypus",
    SINK=TRUE,
    ENSEMBLE.min=0.6,
    MAXENT=0, GBM=0, GBMSTEP=0, RF=1, GLM=1, GLMSTEP=0, GAM=0, 
    GAMSTEP=0, MGCV=0, MGCVFIX=0, EARTH=0, RPART=0, NNET=0, FDA=0, 
    SVM=0, SVME=0, BIOCLIM=1, DOMAIN=1, MAHAL=0,
    PROBIT=TRUE,
    Yweights="BIOMOD", factors="biome",
    PLOTS=FALSE, evaluations.keep=TRUE,
    formulae.defaults=TRUE)    


# instead of providing presence and background locations, provide data.frames
# because 'biome' is a factor, RasterStack and extent need to be provided
# to check for levels in the Training and Testing data set
TrainData1 <- prepareData(x=predictors, p=pres_train, b=backg_train, 
    factors=c("biome"), xy=FALSE)
TestData1 <- prepareData(x=predictors, p=pres_test, b=backg_test, 
    factors=c("biome"), xy=FALSE)
ensemble.factors1 <- ensemble.test(x=predictors, ext=ext, 
    TrainData=TrainData1, TestData=TestData1,
    p=pres_train, a=backg_train, 
    pt=pres_test, at=backg_test,
    species.name="Bradypus",
    SINK=TRUE,
    MAXENT=1, GBM=1, GBMSTEP=0, RF=1, GLM=1, GLMSTEP=1, GAM=1, 
    GAMSTEP=1, MGCV=1, MGCVFIX=1, EARTH=1, RPART=1, NNET=1, FDA=1, 
    SVM=1, SVME=1, BIOCLIM=1, DOMAIN=1, MAHAL=0,
    Yweights="BIOMOD", factors="biome",
    PLOTS=FALSE, evaluations.keep=TRUE)

# compare different methods of calculating ensembles
ensemble.factors2 <- ensemble.test(x=predictors, ext=ext,
    TrainData=TrainData1, TestData=TestData1,
    p=pres_train, a=backg_train, 
    pt=pres_test, at=backg_test,
    species.name="Bradypus",
    SINK=TRUE,
    MAXENT=1, GBM=1, GBMSTEP=0, RF=1, GLM=1, GLMSTEP=1, GAM=1, 
    GAMSTEP=1, MGCV=1, MGCVFIX=1, EARTH=1, RPART=1, NNET=1, FDA=1, 
    SVM=1, SVME=1, BIOCLIM=1, DOMAIN=1, MAHAL=0,
    ENSEMBLE.best=c(4:10), ENSEMBLE.exponent=c(1, 2, 4, 6, 8),
    Yweights="BIOMOD", factors="biome",
    PLOTS=FALSE, evaluations.keep=TRUE)

# test performance of different suitability models
# data are split in 4 subsets, each used once for evaluation
ensemble.nofactors2 <- ensemble.test.splits(x=predictors, ext=ext,
    p=pres, a=background, k=4, 
    layer.drops=c("biome"),
    species.name="Bradypus",
    SINK=TRUE,
    MAXENT=1, GBM=1, GBMSTEP=0, RF=1, GLM=1, GLMSTEP=1, GAM=1, 
    GAMSTEP=1, MGCV=1, MGCVFIX=1, EARTH=1, RPART=1, NNET=1, FDA=1, 
    SVM=1, SVME=1, BIOCLIM=1, DOMAIN=1, MAHAL=0,
    ENSEMBLE.best=0, ENSEMBLE.exponent=c(1, 2, 4, 6, 8),
    ENSEMBLE.min=0.7,
    Yweights="BIOMOD", factors="biome",
    PLOTS=FALSE, formulae.defaults=TRUE,
    GBMSTEP.learning.rate=0.002)
ensemble.nofactors2

# test the result of leaving out one of the variables from the model
# note that positive differences indicate that the model without the variable 
# has higher AUC than the full model
ensemble.variables <- ensemble.drop1(x=predictors, ext=ext,
    p=pres, a=background, k=5,
    layer.drops=c("bio6", "bio1", "bio12"),
    species.name="Bradypus",
    SINK=TRUE,
    difference=TRUE,
    VIF=TRUE, 
    MAXENT=0, GBM=1, GBMSTEP=0, RF=1, GLM=1, GLMSTEP=1, GAM=1, 
    GAMSTEP=1, MGCV=1, MGCVFIX=1, EARTH=1, RPART=1, NNET=1, FDA=1, 
    SVM=1, SVME=1, BIOCLIM=0, DOMAIN=0, MAHAL=0,
    ENSEMBLE.best=0, ENSEMBLE.exponent=c(1, 2, 4, 6, 8),
    ENSEMBLE.min=0.7,
    Yweights="BIOMOD", factors="biome",
    GBMSTEP.learning.rate=0.002)
ensemble.variables

}
}
