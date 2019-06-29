\name{ensemble.batch}
\alias{ensemble.batch}
\alias{ensemble.mean}
\alias{ensemble.plot}

\title{
Suitability mapping based on ensembles of modelling algorithms: batch processing
}

\description{
The main function allows for batch processing of different species and different environmental RasterStacks. The function makes internal calls to \code{\link{ensemble.calibrate.weights}}, \code{\link{ensemble.calibrate.models}} and \code{\link{ensemble.raster}}. 
}

\usage{
ensemble.batch(x = NULL, xn = c(x), 
    species.presence = NULL, species.absence = NULL, 
    presence.min = 20, thin.km = 0.1,
    an = 1000, excludep = FALSE, target.groups = FALSE,
    get.block = FALSE, block.default = runif(1) > 0.5, get.subblocks = FALSE,
    SSB.reduce = FALSE, CIRCLES.d = 250000, 
    k.splits = 4, k.test = 0, 
    n.ensembles = 1,
    VIF.max = 10, VIF.keep = NULL, 
    SINK = FALSE, CATCH.OFF = FALSE, 
    RASTER.format = "raster", RASTER.datatype = "INT2S", RASTER.NAflag = -32767, 
    KML.out = FALSE, KML.maxpixels = 100000, KML.blur = 10,
    models.save = FALSE,
    threshold.method = "spec_sens", threshold.sensitivity = 0.9, 
    threshold.PresenceAbsence = FALSE,
    ENSEMBLE.best = 0, ENSEMBLE.min = 0.7, ENSEMBLE.exponent = 1,
    ENSEMBLE.weight.min = 0.05,
    input.weights = NULL, 
    MAXENT = 1, MAXNET = 1, MAXLIKE = 1, GBM = 1, GBMSTEP = 0, RF = 1, CF = 1,
    GLM = 1, GLMSTEP = 1, GAM = 1, GAMSTEP = 1, MGCV = 1, MGCVFIX = 0, 
    EARTH = 1, RPART = 1, NNET = 1, FDA = 1, SVM = 1 , SVME = 1, GLMNET = 1,
    BIOCLIM.O = 0, BIOCLIM = 1, DOMAIN = 1, MAHAL = 1, MAHAL01 = 1,
    PROBIT = FALSE,
    Yweights = "BIOMOD", 
    layer.drops = NULL, factors = NULL, dummy.vars = NULL,
    formulae.defaults = TRUE, maxit = 100,
    MAXENT.a = NULL, MAXENT.an = 10000, 
    MAXENT.path = paste(getwd(), "/models/maxent", sep=""),
    MAXNET.classes = "default", MAXNET.clamp = FALSE, MAXNET.type = "cloglog",
    MAXLIKE.formula = NULL, MAXLIKE.method = "BFGS",
    GBM.formula = NULL, GBM.n.trees = 2001, 
    GBMSTEP.tree.complexity = 5, GBMSTEP.learning.rate = 0.005, 
    GBMSTEP.bag.fraction = 0.5, GBMSTEP.step.size = 100,
    RF.formula = NULL, RF.ntree = 751, RF.mtry = floor(sqrt(raster::nlayers(x))),
    CF.formula = NULL, CF.ntree = 751, CF.mtry = floor(sqrt(raster::nlayers(x))), 
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
    SVM.formula = NULL, SVME.formula = NULL,
    GLMNET.nlambda = 100, GLMNET.class = FALSE,
    BIOCLIM.O.fraction = 0.9,
    MAHAL.shape = 1)

ensemble.mean(RASTER.species.name = "Species001", RASTER.stack.name = "base",
    positive.filters = c("grd", "_ENSEMBLE_"), negative.filters = c("xml"), 
    RASTER.format = "raster", RASTER.datatype = "INT2S", RASTER.NAflag = -32767,
    KML.out = FALSE, KML.maxpixels = 100000, KML.blur = 10,
    abs.breaks = 6, pres.breaks = 6, sd.breaks = 9,
    p = NULL, a = NULL,
    pt = NULL, at = NULL,
    threshold = -1,
    threshold.method = "spec_sens", threshold.sensitivity = 0.9, 
    threshold.PresenceAbsence = FALSE)

ensemble.plot(RASTER.species.name = "Species001", RASTER.stack.name = "base",
    plot.method=c("suitability", "presence", "count", 
        "consensussuitability", "consensuspresence", "consensuscount", "consensussd"),  
    dev.new.width = 7, dev.new.height = 7,
    main = paste(RASTER.species.name, " ", plot.method, 
        " for ", RASTER.stack.name, sep=""),
    positive.filters = c("grd"), negative.filters = c("xml"), 
    p=NULL, a=NULL,
    threshold = -1,
    threshold.method = "spec_sens", threshold.sensitivity = 0.9, 
    threshold.PresenceAbsence = FALSE,
    abs.breaks = 6, abs.col = NULL,
    pres.breaks = 6, pres.col = NULL,
    sd.breaks = 9, sd.col = NULL,
    absencePresence.col = NULL,
    count.col = NULL,
    maptools.boundaries = TRUE, maptools.col = "dimgrey", ...)
}

\arguments{
  \item{x}{ RasterStack object (\code{\link[raster]{stack}}) containing all layers to calibrate an ensemble.  }
  \item{xn}{ RasterStack object (\code{\link[raster]{stack}}) containing all layers that correspond to explanatory variables of an ensemble calibrated earlier with \code{x}. Several RasterStack objects can be provided in a format as c(stack1, stack2, stack3); these will be used sequentially. See also \code{\link[raster]{predict}}.   }
  \item{species.presence}{ presence points used for calibrating the suitability models, available in 3-column (species, x, y) or (species, lon, lat) dataframe  }
  \item{species.absence}{ background points used for calibrating the suitability models, either available in a 3-column (species, x, y) or (species, lon, lat), or available in a 2-column (x, y) or (lon, lat) dataframe. In case of a 2-column dataframe, the same background locations will be used for all species.  }
  \item{presence.min}{ minimum number of presence locations for the organism (if smaller, no models are fitted).  }
  \item{thin.km}{ Threshold for minimum distance (km) between presence point locations for focal species for model calibrations in each run. A new data set is randomly selected via \code{\link{ensemble.spatialThin}} in each of ensemble run.  }
  \item{an}{ number of background points for calibration to be selected with \code{\link[dismo]{randomPoints}} in case argument \code{a} or \code{species.absence} is missing  }
  \item{excludep}{ parameter that indicates (if \code{TRUE}) that presence points will be excluded from the background points; see also \code{\link[dismo]{randomPoints}} }
  \item{target.groups}{ Parameter that indicates (if \code{TRUE}) that the provided background points (argument \code{a}) represent presence points from a target group sensu Phillips et al. 2009 (these are species that are all collected or observed using the same methods or equipment). Setting the parameter to \code{TRUE} results in selecting the centres of cells of the target groups as background points, while avoiding to select the same cells twice. Via argument \code{excludep}, it is possible to filter out cells with presence observations (argument \code{p}).  }
  \item{get.block}{ if \code{TRUE}, instead of creating k-fold cross-validation subsets randomly (\code{\link[dismo]{kfold}}), create 4 subsets of presence and background locations with \code{\link[ENMeval]{get.block}}.  }
  \item{block.default}{ if \code{FALSE}, instead of making the first division of presence point locations along the y-coordinates (latitude) as in \code{\link[ENMeval]{get.block}}, make the first division along the x-coordinates (longitude).  }
  \item{get.subblocks}{ if \code{TRUE}, then 4 subsets of presence and background locations are generated in a checkerboard configuration by applying \code{\link[ENMeval]{get.block}} to each of the 4 blocks generated by \code{\link[ENMeval]{get.block}} in a first step. }
  \item{SSB.reduce}{ If \code{TRUE}, then new background points that will be used for evaluationg the suitability models will be selected (\code{\link[dismo]{randomPoints}}) in circular neighbourhoods (created with \code{\link[dismo]{circles}}) around presence locations (p and pt). The abbreviation of SSB refers to spatial sorting bias; see also \code{\link[dismo]{ssb}}.  }
  \item{CIRCLES.d}{ Radius in m of circular neighbourhoods (created with \code{\link[dismo]{circles}}) around presence locations (p and pt).  }
  \item{k}{ If larger than 1, the mumber of groups to split between calibration (k-1) and evaluation (1) data sets (for example, \code{k=5} results in 4/5 of presence and background points to be used for calibrating the models, and 1/5 of presence and background points to be used for evaluating the models). See also \code{\link[dismo]{kfold}}.  }
  \item{k.splits}{ If larger than 1, the number of splits for the \code{\link{ensemble.calibrate.weights}} step in batch processing. See also \code{\link[dismo]{kfold}}.  }
  \item{k.test}{ If larger than 1, the mumber of groups to split between calibration (k-1) and evaluation (1) data sets when calibrating the final models (for example, \code{k=5} results in 4/5 of presence and background points to be used for calibrating the models, and 1/5 of presence and background points to be used for evaluating the models). See also \code{\link[dismo]{kfold}}.  }
  \item{n.ensembles}{ If larger than 1, the number of different ensembles generated per species in batch processing. }
  \item{VIF.max}{ Maximum Variance Inflation Factor of variables; see \code{\link[BiodiversityR]{ensemble.VIF}}.  }
  \item{VIF.keep}{ character vector with names of the variables to be kept; see \code{\link[BiodiversityR]{ensemble.VIF}}.  }
  \item{SINK}{ Append the results to a text file in subfolder 'outputs' (if \code{TRUE}). The name of file is based on species names. In case a file already exists, then results are appended. See also \code{\link[base]{sink}}.  }
  \item{CATCH.OFF}{ Disable calls to function \code{\link[base]{tryCatch}}.  }
  \item{RASTER.format}{ Format of the raster files that will be generated. See \code{\link[raster]{writeFormats}} and \code{\link[raster]{writeRaster}}.  }
  \item{RASTER.datatype}{ Format of the raster files that will be generated. See \code{\link[raster]{dataType}} and \code{\link[raster]{writeRaster}}.  }
  \item{RASTER.NAflag}{ Value that is used to store missing data. See \code{\link[raster]{writeRaster}}.  }
  \item{KML.out}{ if \code{FALSE}, then no kml layers (layers that can be shown in Google Earth) are produced. If \code{TRUE}, then kml files will be saved in a subfolder 'kml'.   }
  \item{KML.maxpixels}{ Maximum number of pixels for the PNG image that will be displayed in Google Earth. See also \code{\link[raster]{KML}}.  }
  \item{KML.blur}{ Integer that results in increasing the size of the PNG image by \code{KML.blur^2}, which may help avoid blurring of isolated pixels. See also \code{\link[raster]{KML}}.  }
  \item{models.save}{ Save the list with model details to a file (if \code{TRUE}). The filename will be \code{species.name} with extension \code{.models}; this file will be saved in subfolder of \code{models}. When loading this file, model results will be available as \code{ensemble.models}. }
  \item{threshold.method}{ Method to calculate the threshold between predicted absence and presence; possibilities include \code{spec_sens} (highest sum of the true positive rate and the true negative rate), \code{kappa} (highest kappa value), \code{no_omission} (highest threshold that corresponds to no omission), \code{prevalence} (modeled prevalence is closest to observed prevalence) and \code{equal_sens_spec} (equal true positive rate and true negative rate). See \code{\link[dismo]{threshold}}. Options specific to the BiodiversityR implementation are: \code{threshold.mean} (resulting in calculating the mean value of \code{spec_sens}, \code{equal_sens_spec} and \code{prevalence}) and \code{threshold.min} (resulting in calculating the minimum value of \code{spec_sens}, \code{equal_sens_spec} and \code{prevalence}).   }
  \item{threshold.sensitivity}{ Sensitivity value for \code{threshold.method = 'sensitivity'}. See \code{\link[dismo]{threshold}}.   }
  \item{threshold.PresenceAbsence}{ If \code{TRUE} calculate thresholds with the PresenceAbsence package. See \code{\link[PresenceAbsence]{optimal.thresholds}}.  }
  \item{ENSEMBLE.best}{ The number of individual suitability models to be used in the consensus suitability map (based on a weighted average). In case this parameter is smaller than 1 or larger than the number of positive input weights of individual models, then all individual suitability models with positive input weights are included in the consensus suitability map. In case a vector is provided, \code{ensemble.strategy} is called internally to determine weights for the ensemble model. }
  \item{ENSEMBLE.min}{ The minimum input weight (typically corresponding to AUC values) for a model to be included in the ensemble. In case a vector is provided, function \code{ensemble.strategy} is called internally to determine weights for the ensemble model.  }
  \item{ENSEMBLE.exponent}{ Exponent applied to AUC values to convert AUC values into weights (for example, an exponent of 2 converts input weights of 0.7, 0.8 and 0.9 into 0.7^2=0.49, 0.8^2=0.64 and 0.9^2=0.81). See details. }
  \item{ENSEMBLE.weight.min}{ The minimum output weight for models included in the ensemble, applying to weights that sum to one. Note that \code{ENSEMBLE.min} typically refers to input AUC values.   }
  \item{input.weights}{ array with numeric values for the different modelling algorithms; if \code{NULL} then values provided by parameters such as \code{MAXENT} and \code{GBM} will be used. As an alternative, the output from \code{ensemble.calibrate.weights} can be used. }
  \item{MAXENT}{ Input weight for a maximum entropy model (\code{\link[dismo]{maxent}}). (Only weights > 0 will be used.)  }
  \item{MAXNET}{ number: if larger than 0, then a maximum entropy model (\code{\link[maxnet]{maxnet}}) will be fitted among ensemble }
  \item{MAXLIKE}{ Input weight for a maxlike model (\code{\link[maxlike]{maxlike}}). (Only weights > 0 will be used.)  }
  \item{GBM}{ Input weight for a boosted regression trees model (\code{\link[gbm]{gbm}}). (Only weights > 0 will be used.)  }
  \item{GBMSTEP}{ Input weight for a stepwise boosted regression trees model (\code{\link[dismo]{gbm.step}}). (Only weights > 0 will be used.)  }
  \item{RF}{ Input weight for a random forest model (\code{\link[randomForest]{randomForest}}). (Only weights > 0 will be used.)  }
  \item{CF}{ number: if larger than 0, then a random forest model (\code{\link[party]{cforest}}) will be fitted among ensemble }
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
  \item{GLMNET}{ Input weight for a GLM with lasso or elasticnet regularization (\code{\link[glmnet]{glmnet}}). (Only weights > 0 will be used.) }
  \item{BIOCLIM.O}{ Input weight for the original BIOCLIM algorithm (\code{\link[BiodiversityR]{ensemble.bioclim}}). (Only weights > 0 will be used.) }
  \item{BIOCLIM}{ Input weight for the BIOCLIM algorithm (\code{\link[dismo]{bioclim}}). (Only weights > 0 will be used.)  }
  \item{DOMAIN}{ Input weight for the DOMAIN algorithm (\code{\link[dismo]{domain}}). (Only weights > 0 will be used.)  }
  \item{MAHAL}{ Input weight for the Mahalonobis algorithm (\code{\link[dismo]{mahal}}). (Only weights > 0 will be used.)  }
  \item{MAHAL01}{ Input weight for the Mahalanobis algorithm (\code{\link[dismo]{mahal}}), using a transformation method afterwards whereby output is within the range between 0 and 1. (Only weights > 0 will be used.) }
  \item{PROBIT}{ If \code{TRUE}, then subsequently to the fitting of the individual algorithm (e.g. maximum entropy or GAM) a generalized linear model (\code{\link[stats]{glm}}) with probit link \code{family=binomial(link="probit")} will be fitted to transform the predictions, using the previous predictions as explanatory variable. This transformation results in all model predictions to be probability estimates.  }
  \item{Yweights}{ chooses how cases of presence and background (absence) are weighted; \code{"BIOMOD"} results in equal weighting of all presence and all background cases, \code{"equal"} results in equal weighting of all cases. The user can supply a vector of weights similar to the number of cases in the calibration data set. }
  \item{layer.drops}{ vector that indicates which layers should be removed from RasterStack \code{x}. See also \code{\link[raster]{addLayer}}.  }
  \item{factors}{ vector that indicates which variables are factors; see also \code{\link[dismo]{prepareData}}  }
  \item{dummy.vars}{ vector that indicates which variables are dummy variables (influences formulae suggestions) }
  \item{formulae.defaults}{ Suggest formulae for most of the models (if \code{TRUE}). See also \code{\link{ensemble.formulae}}.   }
  \item{maxit}{ Maximum number of iterations for some of the models. See also \code{\link[stats]{glm.control}}, \code{\link[gam]{gam.control}}, \code{\link[mgcv]{gam.control}} and \code{\link[nnet]{nnet}}.  }
  \item{MAXENT.a}{ background points used for calibrating the maximum entropy model (\code{\link[dismo]{maxent}}), typically available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}.  }
  \item{MAXENT.an}{ number of background points for calibration to be selected with \code{\link[dismo]{randomPoints}} in case argument \code{MAXENT.a} is missing. When used with the \code{ensemble.batch} function, the same background locations will be used for each of the species runs; this implies that for each species, presence locations are not excluded from the background data for this function.   }
  \item{MAXENT.path}{ path to the directory where output files of the maximum entropy model are stored; see also \code{\link[dismo]{maxent}}   }
  \item{MAXNET.classes}{ continuous feature classes, either "default" or any subset of "lqpht" (linear, quadratic, product, hinge, threshold). Note that the "default" option chooses feature classes based on the number of presence locations as "l" (< 10 locations), "lq" (10 - 14 locations), "lqh" (15 - 79 locations) or "lqph" (> 79 locations). See also \code{\link[maxnet]{maxnet}}.   }
  \item{MAXNET.clamp}{ restrict predictors and features to the range seen during model training; see also \code{\link[maxnet]{predict.maxnet}}   }
  \item{MAXNET.type}{ type of response required; see also \code{\link[maxnet]{predict.maxnet}}   }
  \item{MAXLIKE.formula}{ formula for the maxlike algorithm; see also \code{\link[maxlike]{maxlike}}  }
  \item{MAXLIKE.method}{ method for the maxlike algorithm; see also \code{\link[stats]{optim}}  }
  \item{GBM.formula}{ formula for the boosted regression trees algorithm; see also \code{\link[gbm]{gbm}}  }
  \item{GBM.n.trees}{ total number of trees to fit for the boosted regression trees model; see also \code{\link[gbm]{gbm}}  }
  \item{GBMSTEP.tree.complexity}{ complexity of individual trees for stepwise boosted regression trees; see also \code{\link[dismo]{gbm.step}}  }
  \item{GBMSTEP.learning.rate}{ weight applied to individual trees for stepwise boosted regression trees; see also \code{\link[dismo]{gbm.step}} }
  \item{GBMSTEP.bag.fraction}{ proportion of observations used in selecting variables for stepwise boosted regression trees; see also \code{\link[dismo]{gbm.step}}  }
  \item{GBMSTEP.step.size}{ number of trees to add at each cycle for stepwise boosted regression trees (should be small enough to result in a smaller holdout deviance than the initial number of trees [50]); see also \code{\link[dismo]{gbm.step}}  }
  \item{RF.formula}{ formula for the random forest algorithm; see also \code{\link[randomForest]{randomForest}}  }
  \item{RF.ntree}{ number of trees to grow for random forest algorithm; see also \code{\link[randomForest]{randomForest}}  }
  \item{RF.mtry}{ number of variables randomly sampled as candidates at each split for random forest algorithm; see also \code{\link[randomForest]{randomForest}}  }
  \item{CF.formula}{ formula for random forest algorithm; see also \code{\link[party]{cforest}}  }
  \item{CF.ntree}{ number of trees to grow in a forest; see also \code{\link[party]{cforest_control}}  }
  \item{CF.mtry}{ number of input variables randomly sampled as candidates at each node for random forest like algorithms; see also \code{\link[party]{cforest_control}}   }
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
  \item{MGCVFIX.formula}{ formula for the generalized additive model with fixed d.f. regression splines; see also \code{\link[mgcv]{gam}} (the default formulae sets "s(..., fx=TRUE, ...)"; see also \code{\link[mgcv]{s}}) }
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
  \item{GLMNET.nlambda}{ The number of \code{lambda} values; see also \code{\link[glmnet]{glmnet}}  }
  \item{GLMNET.class}{ Use the predicted class to calculate the mean predictions of GLMNET; see also \code{\link[glmnet]{predict.glmnet}}  }
  \item{BIOCLIM.O.fraction}{ Fraction of range representing the optimal limits, default value of 0.9 as in the original BIOCLIM software (\code{\link[BiodiversityR]{ensemble.bioclim}}).  }
  \item{MAHAL.shape}{ parameter that influences the transformation of output values of \code{\link[dismo]{mahal}}.  }

  \item{RASTER.species.name}{ First part of the names of the raster files, expected to identify the modelled species (or organism). }
  \item{RASTER.stack.name}{ Last part of the names of the raster files, expected to identify the predictor stack used. }
  \item{positive.filters}{ vector that indicates parts of filenames for files that will be included in the calculation of the mean probability values  }
  \item{negative.filters}{ vector that indicates parts of filenames for files that will not be included in the calculation of the mean probability values  }
  \item{abs.breaks}{ Number of breaks in the colouring scheme for absence (only applies to \code{suitability} mapping).   }
  \item{pres.breaks}{ Number of breaks in the colouring scheme for presence (only applies to \code{suitability} mapping).  }
  \item{sd.breaks}{ Number of breaks in the colouring scheme for standard deviation (only applies to \code{sd} mapping).  }
  \item{p}{ presence points used for calibrating the suitability models, typically available in 2-column (x, y) or (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{a}{ background points used for calibrating the suitability models, typically available in 2-column (x, y) or (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{pt}{ presence points used for evaluating the suitability models, typically available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}}  }
  \item{at}{ background points used for calibrating the suitability models, typicall available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{threshold}{ Threshold value that will be used to distinguish between presence and absence. If < 0, then a threshold value will be calculated from the provided presence \code{p} and absence \code{a} locations.  }

  \item{plot.method}{ Choice of maps to be plotted: \code{suitability} plots suitability maps, \code{presence} plots presence-absence maps, \code{count} plots count maps (count of number of algorithms or number of ensembles predicting presence) and \code{sd} plots standard deviation maps.  }
  \item{dev.new.width}{ Width for new graphics device (\code{\link[grDevices]{dev.new}}). If < 0, then no new graphics device is opened.  }
  \item{dev.new.height}{ Heigth for new graphics device (\code{\link[grDevices]{dev.new}}). If < 0, then no new graphics device is opened.  }
  \item{main}{ main title for the plots.  }
  \item{abs.col}{ specify colours for absence (see examples on how not to plot areas where the species is predicted absent)   }
  \item{pres.col}{ specify colours for presence  }
  \item{sd.col}{ specify colours for standard deviation   }
  \item{absencePresence.col}{ specify colours for absence - presence maps (see examples on how not to plot areas where the species is predicted absent)   }
  \item{count.col}{ specify colours for number of algorithms or ensembles (see examples on how not to plot areas where the species is predicted absent)   }
  \item{maptools.boundaries}{ If \code{TRUE}, then plot approximate country boundaries \code{\link[maptools]{wrld_simpl}}   }
  \item{maptools.col}{ Colour for approximate country boundaries plotted via \code{\link[maptools]{wrld_simpl}}   }
  \item{...}{ Other items passed to function \code{\link[raster]{plot}}. }

}

\details{
This function allows for batch processing of different species and different environmental RasterStacks. The function makes internal calls to \code{\link{ensemble.spatialThin}}, \code{\link{ensemble.VIF}}, \code{\link{ensemble.calibrate.weights}}, \code{\link{ensemble.calibrate.models}} and \code{\link{ensemble.raster}}. 

Different ensemble runs allow for different random selection of k-fold subsets, background locations or spatial thinning of presence locations.

\code{\link{ensemble.calibrate.weights}} results in a cross-validation procedure whereby the data set is split in calibration and testing subsets and the best weights for the ensemble model are determined (including the possibility for weights = 0).

\code{\link{ensemble.calibrate.models}} is the step whereby models are calibrated using all the available presence data.

\code{\link{ensemble.raster}} is the final step whereby raster layers are produced for the ensemble model.

Function \code{ensemble.mean} results in raster layers that are based on the summary of several ensemble layers: the new ensemble has probability values that are the mean of the probabilities of the different raster layers, the presence-absence threshold is derived for this new ensemble layer, whereas the count reflects the number of ensemble layers where presence was predicted. Note the assumption that input probabilities are scaled between 0 and 1000 (as the output from \code{\link{ensemble.raster}}), whereas thresholds are based on actual probabilities (scaled between 0 and 1). After the mean probability has been calculated, the niche overlap (\code{\link[dismo]{nicheOverlap}}) with the different input layers is calculated. 

Function \code{ensemble.plot} plots suitability, presence-absence or count maps. In the case of suitability maps, the presence-absence threshold needs to be provide as suitabilities smaller than the threshold will be coloured red to orange, whereas suitabilities larger than the threshold will be coloured light blue to dark blue.

}

\value{
The function finally results in ensemble raster layers for each species, including the fitted values for the ensemble model, the estimated presence-absence and the count of the number of submodels prediction presence and absence. 

}

\references{
Kindt R. 2018. Ensemble species distribution modelling with transformed suitability values. Environmental Modelling & Software 100: 136-145. \url{https://doi.org/10.1016/j.envsoft.2017.11.009}

Buisson L, Thuiller W, Casajus N, Lek S and Grenouillet G. 2010. Uncertainty in ensemble forecasting of species distribution. Global Change Biology 16: 1145-1157

Phillips SJ, Dudik M, Elith J et al. 2009. Sample selection bias and presence-only distribution models: implications for background and pseudo-absence data. Ecological Applications 19: 181-197.
}



\author{ Roeland Kindt (World Agroforestry Centre), Eike Luedeling (World Agroforestry Centre) and Evert Thomas (Bioversity International) }


\seealso{
\code{\link{ensemble.calibrate.weights}}, \code{\link{ensemble.calibrate.models}}, \code{\link{ensemble.raster}}
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
pres <- read.table(presence_file, header=TRUE, sep=',')
pres[,1] <- rep("Bradypus", nrow(pres))

# choose background points
background <- randomPoints(predictors, n=1000, extf = 1.00)

# north and south for new predictions (as if new climates)
ext2 <- extent(-90, -32, 0, 23)
predictors2 <- crop(predictors, y=ext2)
predictors2 <- stack(predictors2)
predictors2@title <- "north"

ext3 <- extent(-90, -32, -33, 0)
predictors3 <- crop(predictors, y=ext3)
predictors3 <- stack(predictors3)
predictors3@title <- "south"

# fit 3 ensembles with batch processing, choosing the best ensemble model based on the 
# average weights of 4-fold split of calibration and testing data
# final models use all available presence data and average weights determined by the 
# ensemble.calibrate.weights function (called internally)
# batch processing can handle several species by using 3-column species.presence and 
# species.absence data sets
# note that these calculations can take a while

ensemble.nofactors <- ensemble.batch(x=predictors, 
    xn=c(predictors, predictors2, predictors3),
    species.presence=pres, 
    species.absence=background, 
    k.splits=4, k.test=0, 
    n.ensembles=3, 
    SINK=TRUE, 
    layer.drops=c("biome"),
    ENSEMBLE.best=0, ENSEMBLE.exponent=c(1, 2, 3), 
    ENSEMBLE.min=0.7,
    MAXENT=0, MAXNET=1, MAXLIKE=1, GBM=1, GBMSTEP=0, RF=1, CF=1,
    GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=1, MGCV=1, MGCVFIX=1, 
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, GLMNET=1,
    BIOCLIM.O=1, BIOCLIM=1, DOMAIN=1, MAHAL=0, MAHAL01=1,
    PROBIT=TRUE,
    Yweights="BIOMOD",
    formulae.defaults=TRUE)

# summaries for the 3 ensembles for the species
# summaries are based on files in folders ensemble/suitability, 
# ensemble/presence and ensemble/count
# ensemble.mean is used internally in ensemble.batch

ensemble.mean(RASTER.species.name="Bradypus", RASTER.stack.name="base",
    p=pres, a=background)

# plot mean suitability without specifying colours
plot1 <- ensemble.plot(RASTER.species.name="Bradypus", RASTER.stack.name="base",
    plot.method="consensussuitability",
    p=pres, a=background, abs.breaks=4, pres.breaks=9)
plot1

# only colour the areas where species is predicted to be present
# option is invoked by having no absence breaks
# same colourscheme as \url{http://www.worldagroforestry.org/atlas-central-america}
LAatlascols <- grDevices::colorRampPalette(c("#FFFF80", "#38E009","#1A93AB", "#0C1078"))
plot2 <- ensemble.plot(RASTER.species.name="Bradypus", RASTER.stack.name="base",
    plot.method="consensussuitability",
    p=pres, a=background, abs.breaks=0, pres.breaks=9, pres.col=LAatlascols(8))
plot2

# only colour the areas where species is predicted to be present
# option is invoked by only setting one colour for absence-presence
plot3 <- ensemble.plot(RASTER.species.name="Bradypus", RASTER.stack.name="base",
    plot.method="consensuspresence",
    absencePresence.col=c("#90EE90"))

# only colour presence area by specifying colours > 0
plot4 <- ensemble.plot(RASTER.species.name="Bradypus", RASTER.stack.name="base",
    plot.method="consensuscount",
    count.col=LAatlascols(3))



}

}