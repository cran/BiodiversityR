\name{ensemble.raster}
\alias{ensemble.raster}
\alias{ensemble.batch}
\alias{ensemble.weights}

\title{
Suitability mapping based on ensembles of modelling algorithms: consensus mapping
}

\description{
The basic function \code{ensemble.raster} creates two consensus raster layers, one based on a (weighted) average of different suitability modelling algorithms, and a second one documenting the number of modelling algorithms that predict presence of the focal organisms. Modelling algorithms include maximum entropy (MAXENT), boosted regression trees, random forests, generalized linear models (including stepwise selection of explanatory variables), generalized additive models (including stepwise selection of explanatory variables), multivariate adaptive regression splines, regression trees, artificial neural networks, flexible discriminant analysis, support vector machines, the BIOCLIM algorithm, the DOMAIN algorithm and the Mahalonobis algorithm. These sets of functions were developed in parallel with the \code{biomod2} package, especially for inclusion of the maximum entropy algorithm, but also to allow for a more direct integration with the BiodiversityR package, more direct handling of model formulae and greater focus on mapping. Researchers and students of species distribution are strongly encouraged to familiarize themselves with all the options of the \code{biomod2} and \code{dismo} packages. 
}

\usage{
ensemble.raster(x, p, a = NULL, an = 1000, ext = NULL, k = 0, pt = NULL, at = NULL, xn = x, 
    layer.drops = NULL,
    models.keep = FALSE, 
    RASTER.species.name = "Species001", RASTER.stack.name = xn@title, 
    RASTER.format = "raster", RASTER.datatype = "INT2S", RASTER.NAflag=-32767, 
    RASTER.models.overwrite = TRUE, 
    threshold.method = "spec_sens",
    ENSEMBLE.decay = 1, ENSEMBLE.multiply = TRUE, ENSEMBLE.best = 0, 
    input.weights = NULL,
    MAXENT = 1, GBM = 1, GBMSTEP = 1, RF = 1, GLM = 1, GLMSTEP = 1, GAM = 1, GAMSTEP = 1, MGCV = 1, MGCVFIX = 0,
    EARTH = 1, RPART = 1, NNET = 1, FDA = 1, SVM = 1, SVME = 1, BIOCLIM = 1, DOMAIN = 1, MAHAL = 1, 
    Yweights = "BIOMOD", factors = NULL, dummy.vars = NULL,
    evaluation.strip = TRUE, 
    formulae.defaults = TRUE, maxit = 100,
    MAXENT.OLD = NULL, 
    GBM.formula = NULL, GBM.n.trees = 2000, GBM.OLD = NULL,
    GBMSTEP.gbm.x = 2:(1 + nlayers(x)), GBMSTEP.tree.complexity = 5, GBMSTEP.learning.rate = 0.005, 
    GBMSTEP.bag.fraction = 0.5, GBMSTEP.step.size = 100, GBMSTEP.OLD = NULL, 
    RF.formula = NULL, RF.ntree = 750, RF.mtry = log(nlayers(x)), RF.OLD = NULL, 
    GLM.formula = NULL, GLM.family = binomial(link = "logit"), GLM.OLD = NULL, 
    GLMSTEP.steps = 1000, STEP.formula = NULL, GLMSTEP.scope = NULL, 
    GLMSTEP.k = 2, GLMSTEP.OLD = NULL, GAM.formula = NULL, GAM.family = binomial(link = "logit"), GAM.OLD = NULL, 
    GAMSTEP.steps = 1000, GAMSTEP.scope = NULL, GAMSTEP.OLD = NULL, 
    MGCV.formula = NULL, MGCV.select = FALSE, MGCV.OLD = NULL, 
    MGCVFIX.formula = NULL, MGCVFIX.OLD = NULL, 
    EARTH.formula = NULL, EARTH.glm = list(family = binomial(link = "logit"), maxit = maxit), EARTH.OLD = NULL, 
    RPART.formula = NULL, RPART.xval = 50, RPART.OLD = NULL, 
    NNET.formula = NULL, NNET.size = 8, NNET.decay = 0.01, NNET.OLD = NULL, 
    FDA.formula = NULL, FDA.OLD = NULL, 
    SVM.formula = NULL, SVM.OLD= NULL, 
    SVME.formula = NULL, SVME.OLD= NULL, 
    BIOCLIM.OLD = NULL, DOMAIN.OLD = NULL, MAHAL.OLD = NULL)

ensemble.batch(x, species.presence, species.absence = NULL, an = 1000, ext = NULL, k = 0, k.splits = 5, xn = c(x), 
    layer.drops = NULL,
    digits = 2,
    RASTER.format = "raster", RASTER.datatype = "INT2S", RASTER.NAflag = -32767, 
    threshold.method = "spec_sens", 
    ENSEMBLE.decay = 1, ENSEMBLE.multiply = TRUE, ENSEMBLE.best = 0, 
    input.weights = NULL, 
    MAXENT = 1, GBM = 1, GBMSTEP = 1, RF = 1, GLM = 1, GLMSTEP = 1, GAM = 1, GAMSTEP = 1, MGCV = 1, MGCVFIX = 0,
    EARTH = 1, RPART = 1, NNET = 1, FDA = 1, SVM = 1, SVME = 1, BIOCLIM = 1, DOMAIN = 1, MAHAL = 1, 
    GEODIST = 0, 
    Yweights = "BIOMOD", factors = NULL, dummy.vars = NULL,
    formulae.defaults = TRUE, maxit = 100,
    GBM.formula = NULL, GBM.n.trees = 2000, 
    GBMSTEP.gbm.x = 2:(1 + nlayers(x)), GBMSTEP.tree.complexity = 5, GBMSTEP.learning.rate = 0.005, 
    GBMSTEP.bag.fraction = 0.5, GBMSTEP.step.size = 100,
    RF.formula = NULL, RF.ntree = 750, RF.mtry = log(nlayers(x)), 
    GLM.formula = NULL, GLM.family = binomial(link = "logit"), 
    GLMSTEP.steps = 1000, STEP.formula = NULL, GLMSTEP.scope = NULL, GLMSTEP.k = 2, 
    GAM.formula = NULL, GAM.family = binomial(link = "logit"), 
    GAMSTEP.steps = 1000, GAMSTEP.scope = NULL, 
    MGCV.formula = NULL, MGCV.select = FALSE, 
    MGCVFIX.formula = NULL, 
    EARTH.formula = NULL, EARTH.glm = list(family = binomial(link = "logit"), maxit = maxit), 
    RPART.formula = NULL, RPART.xval = 50, 
    NNET.formula = NULL, NNET.size = 8, NNET.decay = 0.01, FDA.formula = NULL, 
    SVM.formula = NULL, SVME.formula = NULL)

ensemble.weights(weights = c(0.9, 0.8, 0.7, 0.5), decay = 1.5, multiply = TRUE, scale = TRUE, best = 0)
}

\arguments{
  \item{x}{ RasterStack object (\code{\link[raster]{stack}}) containing all layers that correspond to explanatory variables. This RasterStack object is used for calibrating the suitability models.  }
  \item{p}{ presence points used for calibrating the suitability models, typically available in 2-column (x, y) or (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{species.presence}{ presence points used for calibrating the suitability models, available in 3-column (species, x, y) or (species, lon, lat) dataframe  }
  \item{a}{ background points used for calibrating the suitability models, typically available in 2-column (x, y) or (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{species.absence}{ background points used for calibrating the suitability models, either available in a 3-column (species, x, y) or (species, lon, lat), or available in a 2-column (x, y) or (lon, lat) dataframe. In case of a 2-column dataframe, the same background locations will be used for all species.  }
  \item{an}{ number of background points for calibration to be selected with \code{\link[dismo]{randomPoints}} in case argument \code{a} or \code{species.absence} is missing  }
  \item{ext}{ an Extent object to limit the prediction to a sub-region of \code{xn} and the selection of background points to a sub-region of \code{x}, typically provided as c(lonmin, lonmax, latmin, latmax); see also \code{\link[raster]{predict}}, \code{\link[dismo]{randomPoints}} and \code{\link[raster]{extent}} }
  \item{k}{ If larger than 1, the mumber of groups to split between calibration (k-1) and evaluation (1) data sets (for example, \code{k=5} results in 4/5 of presence and background points to be used for calibrating the models, and 1/5 of presence and background points to be used for evaluating the models). See also \code{\link[dismo]{kfold}}.  }
  \item{k.splits}{ If larger than 1, the number of splits for the \code{\link{ensemble.test.splits}} step in batch processing. See also \code{\link[dismo]{kfold}}.  }
  \item{pt}{ presence points used for evaluating the suitability models, typically available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}}  }
  \item{at}{ background points used for calibrating the suitability models, typicall available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{xn}{ RasterStack object (\code{\link[raster]{stack}}) containing all layers that correspond to explanatory variables. This RasterStack object is used for creating the suitability maps corresponding to individual suitability models and consensus models. For batch processing (\code{ensemble.batch}), several RasterStack objects can be provided in a format as c(stack1, stack2, stack3); these will be used sequentially. See also \code{\link[raster]{predict}}.  }
  \item{layer.drops}{ vector that indicates which layers should be removed from RasterStack \code{x}. See also \code{\link[raster]{addLayer}}.  }
  \item{models.keep}{ store the details for each suitability modelling algorithm (if \code{TRUE}). (This may be particularly useful when projecting to different possible future climates.)  }
  \item{RASTER.species.name}{ First part of the names of the raster files that will be generated, potentially an identifier for the modelled species (or organism). }
  \item{RASTER.stack.name}{ Last part of the names of the raster files that will be generated, potentially an identifier for the RasterStack used for predictions (i.e. parameter \code{xn}). }
  \item{RASTER.format}{ Format of the raster files that will be generated. See \code{\link[raster]{writeFormats}} and \code{\link[raster]{writeRaster}}.  }
  \item{RASTER.datatype}{ Format of the raster files that will be generated. See \code{\link[raster]{dataType}} and \code{\link[raster]{writeRaster}}.  }
  \item{RASTER.NAflag}{ Value that is used to store missing data. See \code{\link[raster]{writeRaster}}.  }
  \item{RASTER.models.overwrite}{ Overwrite the raster files that correspond to each suitability modelling algorithm (if \code{TRUE}). (Overwriting actually implies that raster files are created or overwritten that start with "working_").  }
  \item{threshold.method}{ Method to calculate the threshold between predicted absence and presence; one of \code{spec_sens} (highest sum of sum of the true positive rate and the true negative rate), \code{kappa} (highest kappa value) or \code{no_omission} (highest threshold that corresponds to no omission). See See \code{\link[dismo]{threshold}}.   }
  \item{ENSEMBLE.decay}{ Ratio by which the input weights (assumed to correspond to evaluation indices of suitability models) are multiplied after they were sorted from worst to best model based on their input weights. For example, with 4 models and a decay of 2, the input weight of the best model (i.e. the model with the highest input weight) is multiplied by 2*2*2, the input weight of the second best model is multiplied by 2*2, the input weight of the third best model is multiplied by 2, and the worst model is not multiplied.  }
  \item{ENSEMBLE.multiply}{ If \code{TRUE}, then the input weights (assumed to correspond to evaluation indices of suitability models) are multiplied by the \code{ENSEMBLE.decay} to obtain the weighting used by the consensus algorithm (based on a weighted average). If \code{FALSE}, then the weighting used by the consensus algorithm directly corresponds to the \code{ENSEMBLE.decay}.  }
  \item{ENSEMBLE.best}{ The number of individual suitability models to be used in the consensus suitability map (based on a weighted average). In case this parameter is smaller than 1 or larger than the number of positive input weights of individual models, then all individual suitability models with positive input weights are included in the consensus suitability map.  }
  \item{input.weights}{ array with numeric values for the different modelling algorithms; if \code{NULL} then values provided by parameters such as \code{MAXENT} and \code{GBM} will be used. As an alternative, the output from \code{ensemble.test.splits} can be used. }
  \item{MAXENT}{ Input weight for a maximum entropy model (\code{\link[dismo]{maxent}}). (Only weights > 0 will be used.)  }
  \item{GBM}{ Input weight for a boosted regression trees model (\code{\link[gbm]{gbm}}). (Only weights > 0 will be used.)  }
  \item{GBMSTEP}{ Input weight for a stepwise boosted regression trees model (\code{\link[dismo]{gbm.step}}). (Only weights > 0 will be used.)  }
  \item{RF}{ Input weight for a random forest model (\code{\link[randomForest]{randomForest}}). (Only weights > 0 will be used.)  }
  \item{GLM}{ Input weight for a generalized linear model (\code{\link[stats]{glm}}). (Only weights > 0 will be used.)  }
  \item{GLMSTEP}{ Input weight for a stepwise generalized linear model (\code{\link[MASS]{stepAIC}}). (Only weights > 0 will be used.)  }
  \item{GAM}{ Input weight for a generalized additive model (\code{\link[gam]{gam}}). (Only weights > 0 will be used.)  }
  \item{GAMSTEP}{ Input weight for a stepwise generalized additive model (\code{\link[gam]{step.gam}}). (Only weights > 0 will be used.)  }
  \item{MGCV}{ Input weight for a generalized additive model (\code{\link[mgcv]{gam}}). (Only weights > 0 will be used.)  }
  \item{MGCVFIX}{ number: if larger than 0, then a generalized additive model with fixed d.f. regression splines (\code{\link[mgcv]{gam}}) will be fitted among ensemble }
  \item{EARTH}{ Input weight for a multivariate adaptive regression spline model (\code{\link[earth]{earth}}). (Only weights > 0 will be used.)  }
  \item{RPART}{ Input weight for a recursive partioning and regression tree model (\code{\link[rpart]{rpart}}). (Only weights > 0 will be used.)  }
  \item{NNET}{ Input weight for an artificial neural network model (\code{\link[nnet]{nnet}}). (Only weights > 0 will be used.)  }
  \item{FDA}{ Input weight for a flexible discriminant analysis model (\code{\link[mda]{fda}}). (Only weights > 0 will be used.)  }
  \item{SVM}{ Input weight for a support vector machine model (\code{\link[kernlab]{ksvm}}). (Only weights > 0 will be used.)  }
  \item{SVME}{ Input weight for a support vector machine model (\code{\link[e1071]{svm}}). (Only weights > 0 will be used.) }
  \item{BIOCLIM}{ Input weight for the BIOCLIM algorithm (\code{\link[dismo]{bioclim}}). (Only weights > 0 will be used.)  }
  \item{DOMAIN}{ Input weight for the DOMAIN algorithm (\code{\link[dismo]{domain}}). (Only weights > 0 will be used.)  }
  \item{MAHAL}{ Input weight for the Mahalonobis algorithm (\code{\link[dismo]{mahal}}). (Only weights > 0 will be used.)  }
  \item{GEODIST}{ number: if larger than 0, then the geoDist algorithm (\code{\link[dismo]{geoDist}}) will be fitted among ensemble, but only for calculation of AUC values (note that this algorithm does not use environmental layers, and is based only on the distance from presence points used to calibrate this algorithm) }
  \item{Yweights}{ chooses how cases of presence and background (absence) are weighted; \code{"BIOMOD"} results in equal weighting of all presence and all background cases, \code{"equal"} results in equal weighting of all cases. The user can supply a vector of weights similar to the number of cases in the calibration data set. }
  \item{factors}{ vector that indicates which variables are factors; see also \code{\link[dismo]{prepareData}}  }
  \item{dummy.vars}{ vector that indicates which variables are dummy variables (influences formulae suggestions) }
  \item{evaluation.strip}{ Obtain information that can be plotted by the \code{\link{evaluation.strip.plot}} function (if \code{TRUE})  }
  \item{formulae.defaults}{ Suggest formulae for most of the models (if \code{TRUE}). See also \code{\link{ensemble.formulae}}.   }
  \item{maxit}{ Maximum number of iterations for some of the models. See also \code{\link[stats]{glm.control}}, \code{\link[gam]{gam.control}}, \code{\link[mgcv]{gam.control}} and \code{\link[nnet]{nnet}}.  }
  \item{MAXENT.OLD}{ if not \code{NULL}, then do not calibrate a new maximum entropy model, but use this model instead  }
  \item{GBM.formula}{ formula for the boosted regression trees algorithm; see also \code{\link[gbm]{gbm}}  }
  \item{GBM.n.trees}{ total number of trees to fit for the boosted regression trees model; see also \code{\link[gbm]{gbm}}  }
  \item{GBM.OLD}{ if not \code{NULL}, then do not calibrate a new boosted regression tree model, but use this model instead }
  \item{GBMSTEP.gbm.x}{ indices of column numbers with explanatory variables for stepwise boosted regression trees; see also \code{\link[dismo]{gbm.step}}  }
  \item{GBMSTEP.tree.complexity}{ complexity of individual trees for stepwise boosted regression trees; see also \code{\link[dismo]{gbm.step}}  }
  \item{GBMSTEP.learning.rate}{ weight applied to individual trees for stepwise boosted regression trees; see also \code{\link[dismo]{gbm.step}} }
  \item{GBMSTEP.bag.fraction}{ proportion of observations used in selecting variables for stepwise boosted regression trees; see also \code{\link[dismo]{gbm.step}}  }
  \item{GBMSTEP.step.size}{ number of trees to add at each cycle for stepwise boosted regression trees (should be small enough to result in a smaller holdout deviance than the initial number of trees [50]); see also \code{\link[dismo]{gbm.step}}  }
  \item{GBMSTEP.OLD}{ if not \code{NULL}, then do not calibrate a new stepwise boosted regression tree model, but use this model instead }
  \item{RF.formula}{ formula for the random forest algorithm; see also \code{\link[randomForest]{randomForest}}  }
  \item{RF.ntree}{ number of trees to grow for random forest algorithm; see also \code{\link[randomForest]{randomForest}}  }
  \item{RF.mtry}{ number of variables randomly sampled as candidates at each split for random forest algorithm; see also \code{\link[randomForest]{randomForest}}  }
  \item{RF.OLD}{ if not \code{NULL}, then do not calibrate a new random forest model, but use this model instead }
  \item{GLM.formula}{ formula for the generalized linear model; see also \code{\link[stats]{glm}}  }
  \item{GLM.family}{ description of the error distribution and link function for the generalized linear model; see also \code{\link[stats]{glm}}  }
  \item{GLM.OLD}{ if not \code{NULL}, then do not calibrate a new generalized linear model, but use this model instead }
  \item{GLMSTEP.steps}{ maximum number of steps to be considered for stepwise generalized linear model; see also \code{\link[MASS]{stepAIC}}  }
  \item{STEP.formula}{ formula for the "starting model" to be considered for stepwise generalized linear model; see also \code{\link[MASS]{stepAIC}}  }
  \item{GLMSTEP.scope}{ range of models examined in the stepwise search; see also \code{\link[MASS]{stepAIC}}  }
  \item{GLMSTEP.k}{ multiple of the number of degrees of freedom used for the penalty (only k = 2 gives the genuine AIC); see also \code{\link[MASS]{stepAIC}}  }
  \item{GLMSTEP.OLD}{ if not \code{NULL}, then do not calibrate a new stepwise generalized linear model, but use this model instead }
  \item{GAM.formula}{ formula for the generalized additive model; see also \code{\link[gam]{gam}}  }
  \item{GAM.family}{ description of the error distribution and link function for the generalized additive model; see also \code{\link[gam]{gam}}  }
  \item{GAM.OLD}{ if not \code{NULL}, then do not calibrate a new generalized additive model, but use this model instead }
  \item{GAMSTEP.steps}{ maximum number of steps to be considered in the stepwise generalized additive model; see also \code{\link[gam]{step.gam}}  }
  \item{GAMSTEP.scope}{ range of models examined in the step-wise search n the stepwise generalized additive model; see also \code{\link[gam]{step.gam}}  }
  \item{GAMSTEP.OLD}{ if not \code{NULL}, then do not calibrate a new stepwise generalized additive model, but use this model instead }
  \item{MGCV.formula}{ formula for the generalized additive model; see also \code{\link[mgcv]{gam}}  }
  \item{MGCV.select}{ if \code{TRUE}, then the smoothing parameter estimation that is part of fitting can completely remove terms from the model; see also \code{\link[mgcv]{gam}}  }
  \item{MGCV.OLD}{ if not \code{NULL}, then do not calibrate a new generalized additive model, but use this model instead }
  \item{MGCVFIX.formula}{ formula for the generalized additive model with fixed d.f. regression splines; see also \code{\link[mgcv]{gam}} (the default formulae sets "s(..., fx=TRUE, ...)"; see also \code{\link[mgcv]{s}}) }
  \item{MGCVFIX.OLD}{ if not \code{NULL}, then do not calibrate a new generalized additive model with fixed d.f. regression splines, but use this model instead }
  \item{EARTH.formula}{ formula for the multivariate adaptive regression spline model; see also \code{\link[earth]{earth}}  }
  \item{EARTH.glm}{ list of arguments to pass on to \code{\link[stats]{glm}}; see also \code{\link[earth]{earth}}  }
  \item{EARTH.OLD}{ if not \code{NULL}, then do not calibrate a new multivariate adaptive regression spline model, but use this model instead }
  \item{RPART.formula}{ formula for the recursive partioning and regression tree model; see also \code{\link[rpart]{rpart}}  }
  \item{RPART.xval}{ number of cross-validations for the recursive partioning and regression tree model; see also \code{\link[rpart]{rpart.control}}  }
  \item{RPART.OLD}{ if not \code{NULL}, then do not calibrate a new recursive partioning and regression tree model, but use this model instead }
  \item{NNET.formula}{ formula for the artificial neural network model; see also \code{\link[nnet]{nnet}}  }
  \item{NNET.size}{ number of units in the hidden layer for the artificial neural network model; see also \code{\link[nnet]{nnet}}  }
  \item{NNET.decay}{ parameter of weight decay for the artificial neural network model; see also \code{\link[nnet]{nnet}}  }
  \item{NNET.OLD}{ if not \code{NULL}, then do not calibrate a new artificial neural network model, but use this model instead }
  \item{FDA.formula}{ formula for the flexible discriminant analysis model; see also \code{\link[mda]{fda}}  }
  \item{FDA.OLD}{ if not \code{NULL}, then do not calibrate a new flexible discriminant model, but use this model instead }
  \item{SVM.formula}{ formula for the support vector machine model; see also \code{\link[kernlab]{ksvm}}  }
  \item{SVM.OLD}{ if not \code{NULL}, then do not calibrate a new support vector machine model, but use this model instead }
  \item{SVME.formula}{ formula for the support vector machine model; see also \code{\link[e1071]{svm}}  }
  \item{SVME.OLD}{ if not \code{NULL}, then do not calibrate a new support vector machine model, but use this model instead }
  \item{BIOCLIM.OLD}{ if not \code{NULL}, then do not calibrate a new BIOCLIM model, but use this model instead }
  \item{DOMAIN.OLD}{ if not \code{NULL}, then do not calibrate a new DOMAIN model, but use this model instead }
  \item{MAHAL.OLD}{ if not \code{NULL}, then do not calibrate a new Mahalonobis model, but use this model instead }
  \item{digits}{ number of decimal place after multiplying AUC results by 100 (percentages); see also \code{\link[base]{round}}  }
  \item{weights}{ input weights for the \code{ensemble.weights} function }
  \item{decay}{ Ratio by which the input weights are multiplied after they were sorted. For example, with 4 input weights and a decay of 2, the largest input weight is multiplied by 2*2*2, the second largest input weight is multiplied by 2*2, the third largest input weight is multiplied by 2, and the smallest input weight is not multiplied.  }
  \item{multiply}{ If \code{TRUE}, then the input weights are multiplied by the \code{decay} to obtain the output weighting. If \code{FALSE}, then the output weighting directly corresponds to the \code{ENSEMBLE.decay}.  }
  \item{best}{ The number of final weights. In case this parameter is smaller than 1 or larger than the number of positive input weights of individual models, then this parameter is ignored.  }
  \item{scale}{ scale the output weights so that the sum of output weights equals 1 (if \code{TRUE}, default for \code{ensemble.raster}) }
}

\details{
The basic function \code{ensemble.raster} fits individual suitability models for all models with positive input weights. In subfolder "models" of the working directory, suitability maps for the individual suitability modelling algorithms are stored. In subfolder "ensembles", a consensus suitability map based on a weighted average of individual suitability models is stored (the name of this map contains "_ENSEMBLE_"). In subfolder "ensembles", a presence-absence (1-0) map is also provided (based on the consensus suitability map and the defined method for calculating the threshold); part of the name for this map contains "_PRESENCE_". Also in subfolder "ensembles", a consensus suitability map (with part of the name containing "_COUNT") based on the number of individual suitability models that predict presence of the focal organism is stored.

Note that values in suitability maps are integer values that were calculated by multiplying probabilities by 1000 (see also \code{\link[base]{trunc}}).

As the Mahalanobis function does not always provide values within a range of 0 - 1, the output values are rescaled by subtracting the value of 2 from each prediction, calculating the absolute value and finally the reciprocal value. As this rescaling method does not estimate probabilies, inclusion in the calculation of a (weighted) average of ensemble probabilities may be problematic (the same applies to other distance-based methods). 

The \code{ensemble.batch} function allows for batch processing, using presence (possibly absence) locations for several species. The function also allows predictions for more than one rasterStack object. The function first calls the \code{ensemble.test.splits} function to determine the average AUC (parameter \code{k.splits} determines the number of splits); thereafter the average AUC is used as input weight in the \code{ensemble.raster} function. Note that more information is provided when exploring all the options of \code{ensemble.raster} by a species-by-species and raster-by-raster process.   

When reporting results and as part of file names, names for \code{rasterStack} used for calibrations or predictions are obtained from the "title" slot. In case that these slots do not contain a name, then names such as "base" and "stack2" are given for these slots to differentiate between the \code{rasterStack}.

The \code{ensemble.weights} function is used internally by the \code{ensemble.raster} function, using the input weights for the different suitability modelling algorithms. Ties between input weights result in the same output weights.
}

\value{
The basic function \code{ensemble.raster} mainly results in the creation of raster layers that correspond to fitted probabilities of presence of individual suitability models (in folder "models") and consensus models (in folder "ensembles"), and the number of suitability models that predict presence (in folder "ensembles"). Prediction of presence is based on a threshold usually defined by maximizing the sum of the true presence and true absence rates (see \code{threshold.method} and also \code{\link[dismo]{ModelEvaluation}}).

If desired by the user, the \code{ensemble.raster} function also saves details of fitted suitability models or data that can be plotted with the \code{\link{evaluation.strip.plot}} function.

}

\references{
Buisson L, Thuiller W, Casajus N, Lek S and Grenouillet G. 2010. Uncertainty in ensemble forecasting of species distribution. Global Change Biology 16: 1145-1157
}



\author{ Roeland Kindt (World Agroforestry Centre) and Eike Luedeling (World Agroforestry Centre) }


\seealso{
\code{\link{evaluation.strip.plot}}, \code{\link{ensemble.test}}, \code{\link{ensemble.test.splits}}
}

\examples{
\dontrun{
# based on examples in the dismo package
# simplified example of ensemble modelling with 4 modeling algorithms

# get predictor variables
library(dismo)
predictors <- stack(list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), pattern='grd', full.names=TRUE ))
predictors


# presence points
presence <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
pres <- read.table(presence, header=TRUE, sep=',')[,-1]

# choose background points
ext = extent(-90, -32, -33, 23)
background <- randomPoints(predictors, n=1000, ext=ext, extf = 1.00)


# formulae for random forest and generalized linear model
rfformula <- as.formula(pb ~ bio1+bio5+bio6+bio7+bio8+bio12+bio16+bio17)

glmformula <- as.formula(pb~bio1 + I(bio1^2) + I(bio1^3) + bio12 + I(bio12^2) + I(bio12^3) + 
    bio16 + I(bio16^2) + I(bio16^3) + bio17 + I(bio17^2) + I(bio17^3) + 
    bio5 + I(bio5^2) + I(bio5^3) + bio6 + I(bio6^2) + I(bio6^3) + 
    bio7 + I(bio7^2) + I(bio7^3) + bio8 + I(bio8^2) + I(bio8^3) )

# if desired, change working directory where subfolders of "models" and "ensembles" will be created

# fit four ensemble models
# weights are based on previous run of ensemble.test.splits function
# all presence points are utilized
ensemble.nofactors <- ensemble.raster(x=predictors, p=pres, a=background, k=0, ext=ext,
    layer.drops=c("biome"),
    MAXENT=0, GBM=0, GBMSTEP=0, RF=0.8538, GLM=0.7948, GLMSTEP=0, GAM=0, GAMSTEP=0,
    MGCV=0, EARTH=0, RPART=0, NNET=0, FDA=0, SVM=0, SVME=0, BIOCLIM=0.7031, DOMAIN=0.7265, MAHAL=0,
    Yweights="BIOMOD", factors=NULL,
    models.keep=TRUE, evaluation.strip=TRUE,
    RF.formula=rfformula,
    GLM.formula=glmformula)

evaluation.strip.plot(ensemble.nofactors$evaluation.strip, model="ENSEMBLE", type="o", col="red")

# split data in 5 subsets, use 1/5 for evaluation and 4/5 for calibration
# GAMSTEP option is not available at the moment (problems with assignments)

ensemble.nofactors <- ensemble.raster(x=predictors, p=pres, a=background, k=5, ext=ext,
    layer.drops=c("biome"),
    MAXENT=1, GBM=1, GBMSTEP=1, RF=1, GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=0, MGCV=1, 
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, BIOCLIM=1, DOMAIN=1, MAHAL=1,
    Yweights="BIOMOD",
    models.keep=TRUE, evaluation.strip=TRUE,
    formulae.defaults=TRUE)

evaluation.strip.plot(ensemble.nofactors$evaluation.strip, model="ENSEMBLE", type="o", col="red")
}

}