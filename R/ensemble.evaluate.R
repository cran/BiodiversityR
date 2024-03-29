`ensemble.SEDI` <- function(
    TPR, FPR, small=1e-9) 
{
# Equation 2 in Ferro, C.A. and D.B. Stephenson (2011) Extremal Dependence Indices: Improved Verification Measures for Deterministic Forecasts of Rare 
#     Binary Events. Wea. Forecasting, 26, 699-713, https://doi.org/10.1175/WAF-D-10-05030.1
# Equation 7 in Wunderlich RF, Lin Y-P, Anthony J, Petway JR (2019) Two alternative evaluation metrics to replace the true skill statistic in the assessment 
#     of species distribution models. Nature Conservation 35: 97-116. https://doi.org/10.3897/natureconservation.35.33918
# Zeroes are substituted by small number (1e-9) as discussed by Wunderlich et al. and as implemented in https://github.com/RWunderlich/SEDI/blob/master/R/sedi.R

    TPR[TPR==0] <- small
    FPR[FPR==0] <- small
    TNR <- 1-FPR
    FNR <- 1-TPR
    TNR[TNR==0] <- small
    FNR[FNR==0] <- small
    s <- (log(FPR) - log(TPR) - log(TNR) + log(FNR)) / (log(FPR) + log(TPR) + log(TNR) + log(FNR))
    return(s)
}


# Inspired by: 
# https://github.com/kerickson22/SDMs_for_rare_species_modeling/blob/main/code/00b_Constants.R

#computeTjurR2 = function(Y, predY) {
#
#    R2 = mean(predY[which(Y == 1)]) - mean(predY[which(Y == 0)])
#    return(R2)
#}

`ensemble.Tjur` <- function(
    eval
)
{
    result <- mean(eval@presence) - mean(eval@absence)
    return(result)
}


`ensemble.evaluate` <- function(
    eval, fixed.threshold=NULL, eval.train=NULL) 
{
    if (inherits(eval, "ModelEvaluation") == FALSE) {stop(paste("Please provide a ModelEvaluation object", "\n", sep = ""))}
    if (is.null(fixed.threshold) == T) {
        if (is.null(eval.train) == T) {
            stop(paste("Please provide a ModelEvaluation object for the training data", "\n", sep = ""))
        }else{
            fixed.threshold <- eval.train@t[which.max(eval.train@TPR + eval.train@TNR)]
            cat(paste("Calculated fixed threshold of ", fixed.threshold, " corresponding to highest sum of sensitivity and specificity", "\n", sep = ""))
        }
    }
    result <- as.numeric(rep(NA, 9))
    names(result) <- c("AUC", "TSS", "SEDI", "TSS.fixed", "SEDI.fixed", "FNR.fixed", "MCR.fixed", "AUCdiff", "Tjur")
    result["AUC"] <- eval@auc
    tss <- eval@TPR - eval@FPR
    result["TSS"] <- max(tss)
    sedi <- ensemble.SEDI(eval@TPR, eval@FPR)
    result["SEDI"] <- max(sedi)
    result["TSS.fixed"] <- tss[which(eval@t >= fixed.threshold)][1]
    result["SEDI.fixed"] <- sedi[which(eval@t >= fixed.threshold)][1]
    result["FNR.fixed"] <- eval@FNR[which(eval@t >= fixed.threshold)][1]
    result["MCR.fixed"] <- eval@MCR[which(eval@t >= fixed.threshold)][1]
    if (is.null(eval.train) == T) {
        cat(paste("Please provide a ModelEvaluation object for calculating AUCdiff", "\n", sep = ""))
    }else{
       result["AUCdiff"] <- eval.train@auc - eval@auc
    }
    result["Tjur"] <- ensemble.Tjur(eval)
    return(result)
}

