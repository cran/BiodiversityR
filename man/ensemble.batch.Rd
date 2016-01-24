\name{ensemble.batch}
\alias{ensemble.batch}
\alias{ensemble.mean}

\title{
Suitability mapping based on ensembles of modelling algorithms: batch processing
}

\description{
The main function allows for batch processing of different species and different environmental RasterStacks. The function makes internal calls to \code{\link{ensemble.test.splits}}, \code{\link{ensemble.test}} and \code{\link{ensemble.raster}}. 
}

\usage{
ensemble.batch(x = NULL, xn = c(x), ext = NULL, 
    species.presence = NULL, species.absence = NULL, 
    presence.min = 20,
    an = 1000, excludep = FALSE, CIRCLES.at = FALSE, CIRCLES.d = 100000, 
    k.splits = 4, k.test = 0, 
    n.ensembles = 1, 
    SINK = FALSE,
    RASTER.format = "raster", RASTER.datatype = "INT2S", RASTER.NAflag = -32767, 
    KML.out = FALSE, KML.maxpixels = 100000, KML.blur = 10,
    models.save = FALSE,
    threshold.method = "spec_sens", threshold.sensitivity = 0.9, 
    threshold.PresenceAbsence = FALSE,
    ENSEMBLE.best = 0, ENSEMBLE.min = 0.7, ENSEMBLE.exponent = 1,
    input.weights = NULL, 
    MAXENT = 1, GBM = 1, GBMSTEP = 1, RF = 1, GLM = 1, GLMSTEP = 1, GAM = 1, 
    GAMSTEP = 1, MGCV = 1, MGCVFIX = 0, EARTH = 1, RPART = 1, NNET = 1, 
    FDA = 1, SVM = 1, SVME = 1, BIOCLIM = 1, DOMAIN = 1, MAHAL = 1, 
    PROBIT = FALSE, AUC.weights = TRUE,
    Yweights = "BIOMOD", 
    layer.drops = NULL, factors = NULL, dummy.vars = NULL,
    formulae.defaults = TRUE, maxit = 100,
    MAXENT.a = NULL, MAXENT.an = 10000, MAXENT.BackData = NULL, 
    MAXENT.path = paste(getwd(), "/models/maxent", sep=""),
    GBM.formula = NULL, GBM.n.trees = 2001, 
    GBMSTEP.gbm.x = 2:(1 + raster::nlayers(x)), 
    GBMSTEP.tree.complexity = 5, GBMSTEP.learning.rate = 0.005, 
    GBMSTEP.bag.fraction = 0.5, GBMSTEP.step.size = 100,
    RF.formula = NULL, RF.ntree = 751, RF.mtry = floor(sqrt(raster::nlayers(x))), 
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
    MAHAL.shape = 1)

ensemble.mean(RASTER.species.name = "Species001", RASTER.stack.name = "base",
    positive.filters = c("grd", "_ENSEMBLE_"), negative.filters = c("xml"), 
    RASTER.format = "raster", RASTER.datatype = "INT2S", RASTER.NAflag = -32767,
    KML.out = FALSE, KML.maxpixels = 100000, KML.blur = 10,
    p = NULL, a = NULL,
    pt = NULL, at = NULL,
    threshold = -1,
    threshold.method = "spec_sens", threshold.sensitivity = 0.9, 
    threshold.PresenceAbsence = FALSE)

}

\arguments{
  \item{x}{ RasterStack object (\code{\link[raster]{stack}}) containing all layers to calibrate an ensemble.  }
  \item{xn}{ RasterStack object (\code{\link[raster]{stack}}) containing all layers that correspond to explanatory variables of an ensemble calibrated earlier with \code{x}. Several RasterStack objects can be provided in a format as c(stack1, stack2, stack3); these will be used sequentially. See also \code{\link[raster]{predict}}.   }
  \item{ext}{ an Extent object to limit the prediction to a sub-region of \code{xn} and the selection of background points to a sub-region of \code{x}, typically provided as c(lonmin, lonmax, latmin, latmax); see also \code{\link[raster]{predict}}, \code{\link[dismo]{randomPoints}} and \code{\link[raster]{extent}} }
  \item{species.presence}{ presence points used for calibrating the suitability models, available in 3-column (species, x, y) or (species, lon, lat) dataframe  }
  \item{species.absence}{ background points used for calibrating the suitability models, either available in a 3-column (species, x, y) or (species, lon, lat), or available in a 2-column (x, y) or (lon, lat) dataframe. In case of a 2-column dataframe, the same background locations will be used for all species.  }
  \item{presence.min}{ minimum number of presence locations for the organism (if smaller, no models are fitted).  }
  \item{an}{ number of background points for calibration to be selected with \code{\link[dismo]{randomPoints}} in case argument \code{a} or \code{species.absence} is missing  }
  \item{excludep}{ parameter that indicates (if \code{TRUE}) that presence points will be excluded from the background points; see also \code{\link[dismo]{randomPoints}} }
  \item{CIRCLES.at}{ If \code{TRUE}, then new background points that will be used for evaluationg the suitability models will be selected (\code{\link[dismo]{randomPoints}}) in circular neighbourhoods (created with \code{\link[dismo]{circles}}) around presence locations (p and pt).  }
  \item{CIRCLES.d}{ Radius in m of circular neighbourhoods (created with \code{\link[dismo]{circles}}) around presence locations (p and pt).  }
  \item{k}{ If larger than 1, the mumber of groups to split between calibration (k-1) and evaluation (1) data sets (for example, \code{k=5} results in 4/5 of presence and background points to be used for calibrating the models, and 1/5 of presence and background points to be used for evaluating the models). See also \code{\link[dismo]{kfold}}.  }
  \item{k.splits}{ If larger than 1, the number of splits for the \code{\link{ensemble.test.splits}} step in batch processing. See also \code{\link[dismo]{kfold}}.  }
  \item{k.test}{ If larger than 1, the mumber of groups to split between calibration (k-1) and evaluation (1) data sets when calibrating the final models (for example, \code{k=5} results in 4/5 of presence and background points to be used for calibrating the models, and 1/5 of presence and background points to be used for evaluating the models). See also \code{\link[dismo]{kfold}}.  }
  \item{n.ensembles}{ If larger than 1, the number of different ensembles generated per species in batch processing. }
  \item{SINK}{ Append the results to a text file in subfolder 'outputs' (if \code{TRUE}). The name of file is based on species names. In case a file already exists, then results are appended. See also \code{\link[base]{sink}}.  }
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
  \item{ENSEMBLE.best}{ The number of individual suitability models to be used in the consensus suitability map (based on a weighted average). In case this parameter is smaller than 1 or larger than the number of positive input weights of individual models, then all individual suitability models with positive input weights are included in the consensus suitability map.  }
  \item{ENSEMBLE.min}{ The minimum input weight (assumed to correspond to evaluation indices of suitability models) for a model to be included in the ensemble.  }
  \item{ENSEMBLE.exponent}{ Exponent applied to AUC values to convert AUC values into weights (for example, an exponent of 2 converts input weights of 0.7, 0.8 and 0.9 into 0.7^2=0.49, 0.8^2=0.64 and 0.9^2=0.81). See \code{\link{ensemble.test}} for details.  }
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
  \item{PROBIT}{ If \code{TRUE}, then subsequently to the fitting of the individual algorithm (e.g. maximum entropy or GAM) a generalized linear model (\code{\link[stats]{glm}}) with probit link \code{family=binomial(link="probit")} will be fitted to transform the predictions, using the previous predictions as explanatory variable. This transformation results in all model predictions to be probability estimates.  }
  \item{AUC.weights}{ If \code{TRUE}, then use the average of the AUC for the different submodels in the different crossvalidation runs as weights for the 'full' ensemble model. See \code{\link{ensemble.test.splits}} for details.  }
  \item{Yweights}{ chooses how cases of presence and background (absence) are weighted; \code{"BIOMOD"} results in equal weighting of all presence and all background cases, \code{"equal"} results in equal weighting of all cases. The user can supply a vector of weights similar to the number of cases in the calibration data set. }
  \item{layer.drops}{ vector that indicates which layers should be removed from RasterStack \code{x}. See also \code{\link[raster]{addLayer}}.  }
  \item{factors}{ vector that indicates which variables are factors; see also \code{\link[dismo]{prepareData}}  }
  \item{dummy.vars}{ vector that indicates which variables are dummy variables (influences formulae suggestions) }
  \item{formulae.defaults}{ Suggest formulae for most of the models (if \code{TRUE}). See also \code{\link{ensemble.formulae}}.   }
  \item{maxit}{ Maximum number of iterations for some of the models. See also \code{\link[stats]{glm.control}}, \code{\link[gam]{gam.control}}, \code{\link[mgcv]{gam.control}} and \code{\link[nnet]{nnet}}.  }
  \item{MAXENT.a}{ background points used for calibrating the maximum entropy model (\code{\link[dismo]{maxent}}), typically available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}. Ignored if \code{MAXENT.BackData} is provided.  }
  \item{MAXENT.an}{ number of background points for calibration to be selected with \code{\link[dismo]{randomPoints}} in case argument \code{MAXENT.a} is missing. When used with the \code{ensemble.batch} function, the same background locations will be used for each of the species runs; this implies that for each species, presence locations are not excluded from the background data for this function.   }
  \item{MAXENT.BackData}{ dataframe containing explanatory variables for the background locations. This information will be used for calibrating the maximum entropy model (\code{\link[dismo]{maxent}}). When used with the \code{ensemble.batch} function, the same background locations will be used for each of the cross-validation runs; this is based on the assumption that a large number (~10000) of background locations are used.  }
  \item{MAXENT.path}{ path to the directory where output files of the maximum entropy model are stored; see also \code{\link[dismo]{maxent}}   }
  \item{GBM.formula}{ formula for the boosted regression trees algorithm; see also \code{\link[gbm]{gbm}}  }
  \item{GBM.n.trees}{ total number of trees to fit for the boosted regression trees model; see also \code{\link[gbm]{gbm}}  }
  \item{GBMSTEP.gbm.x}{ indices of column numbers with explanatory variables for stepwise boosted regression trees; see also \code{\link[dismo]{gbm.step}}  }
  \item{GBMSTEP.tree.complexity}{ complexity of individual trees for stepwise boosted regression trees; see also \code{\link[dismo]{gbm.step}}  }
  \item{GBMSTEP.learning.rate}{ weight applied to individual trees for stepwise boosted regression trees; see also \code{\link[dismo]{gbm.step}} }
  \item{GBMSTEP.bag.fraction}{ proportion of observations used in selecting variables for stepwise boosted regression trees; see also \code{\link[dismo]{gbm.step}}  }
  \item{GBMSTEP.step.size}{ number of trees to add at each cycle for stepwise boosted regression trees (should be small enough to result in a smaller holdout deviance than the initial number of trees [50]); see also \code{\link[dismo]{gbm.step}}  }
  \item{RF.formula}{ formula for the random forest algorithm; see also \code{\link[randomForest]{randomForest}}  }
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
  \item{MAHAL.shape}{ parameter that influences the transformation of output values of \code{\link[dismo]{mahal}}. See details section.  }

  \item{RASTER.species.name}{ First part of the names of the raster files, expected to identify the modelled species (or organism). }
  \item{RASTER.stack.name}{ Last part of the names of the raster files, expected to identify the predictor stack used. }
  \item{positive.filters}{ vector that indicates parts of filenames for files that will be included in the calculation of the mean probability values  }
  \item{negative.filters}{ vector that indicates parts of filenames for files that will not be included in the calculation of the mean probability values  }
  \item{p}{ presence points used for calibrating the suitability models, typically available in 2-column (x, y) or (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{a}{ background points used for calibrating the suitability models, typically available in 2-column (x, y) or (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{pt}{ presence points used for evaluating the suitability models, typically available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}}  }
  \item{at}{ background points used for calibrating the suitability models, typicall available in 2-column (lon, lat) dataframe; see also \code{\link[dismo]{prepareData}} and \code{\link[raster]{extract}}  }
  \item{threshold}{ Threshold value that will be used to distinguish between presence and absence. If < 0, then a threshold value will be calculated from the provided presence \code{p} and absence \code{a} locations.  }

}

\details{
This function allows for batch processing of different species and different environmental RasterStacks. The function makes internal calls to \code{\link{ensemble.test.splits}}, \code{\link{ensemble.test}} and \code{\link{ensemble.raster}}. 

\code{\link{ensemble.test.splits}} results in a cross-validation procedure whereby the data set is split in calibration and testing subsets and the best weights for the ensemble model are determined (including the possibility for weights = 0).

\code{\link{ensemble.test}} is the step whereby models are calibrated using all the available presence data.

\code{\link{ensemble.raster}} is the final step whereby raster layers are produced for the ensemble model.

Function \code{ensemble.mean} results in raster layers that are based on the summary of several ensemble layers: the new ensemble has probability values that are the mean of the probabilities of the different raster layers, the presence-absence threshold is derived for this new ensemble layer, whereas the count reflects the number of ensemble layers where presence was predicted. Note the assumption that input probabilities are scaled between 0 and 1000 (as the output from \code{\link{ensemble.raster}}), whereas thresholds are based on actual probabilities (scaled between 0 and 1).

}

\value{
The function finally results in ensemble raster layers for each species, including the fitted values for the ensemble model, the estimated presence-absence and the count of the number of submodels prediction presence and absence. 

}

\references{
Buisson L, Thuiller W, Casajus N, Lek S and Grenouillet G. 2010. Uncertainty in ensemble forecasting of species distribution. Global Change Biology 16: 1145-1157
}



\author{ Roeland Kindt (World Agroforestry Centre), Eike Luedeling (World Agroforestry Centre) and Evert Thomas (Bioversity International) }


\seealso{
\code{\link{ensemble.test.splits}}, \code{\link{ensemble.test}}, \code{\link{ensemble.raster}}
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
ext <- extent(-90, -32, -33, 23)
background <- randomPoints(predictors, n=1000, ext=ext, extf = 1.00)

# fit 3 ensembles with batch processing, choosing the best ensemble model based on the 
# average AUC of 4-fold split of calibration and testing data
# final models use all available presence data and average weights determined by the 
# ensemble.test.splits function (called internally)
# batch processing can handle several species by using 3-column species.presence and 
# species.absence data sets
# note that these calculations can take a while

ensemble.nofactors <- ensemble.batch(x=predictors, ext=ext,
    species.presence=pres, 
    species.absence=background, 
    k.splits=4, k.test=0, 
    n.ensembles=3, 
    SINK=TRUE,
    layer.drops=c("biome"),
    ENSEMBLE.best=0, ENSEMBLE.exponent=c(1, 2, 4, 6, 8), 
    ENSEMBLE.min=0.7,
    MAXENT=1, GBM=1, GBMSTEP=0, RF=1, GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=0, MGCV=1, 
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, BIOCLIM=1, DOMAIN=1, MAHAL=0,
    Yweights="BIOMOD",
    formulae.defaults=TRUE)

# summaries for the 3 ensembles for the species
# summaries are based on files in folders ensemble, ensemble/presence and 
# ensemble/count
pres1 <- pres[,-1]
ensemble.mean(RASTER.species.name="Bradypus", 
    p=pres1, a=background, 
    KML.out=T)

}

}