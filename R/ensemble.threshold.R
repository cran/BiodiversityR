`ensemble.threshold` <- function(
    eval, threshold.method="spec_sens", threshold.sensitivity=0.9, 
    threshold.PresenceAbsence=FALSE, Pres, Abs) 
{
    if (threshold.method == "threshold.min") {threshold.method <- "threshold2013.min"}
    if (threshold.method == "threshold.mean") {threshold.method <- "threshold2013.mean"}
    if (threshold.method %in% c("threshold2013.min", "threshold2013.mean", "threshold2005.min", "threshold2005.mean")) {threshold.PresenceAbsence <- TRUE}

    if (threshold.PresenceAbsence == T) {
        Pres2 <- cbind(rep(1, length(Pres)), Pres)
        Abs2 <- cbind(rep(0, length(Abs)), Abs)
        data1 <- rbind(Pres2, Abs2)
        data2 <- cbind(seq(1:nrow(data1)), data1)
        auc.value <- PresenceAbsence::auc(data2, st.dev=F)
        cat(paste("\n", "AUC from PresenceAbsence package (also used to calculate threshold): ", auc.value, "\n", sep = ""))
        if (threshold.method=="kappa") {threshold.method <- "MaxKappa"}
        if (threshold.method=="spec_sens") {threshold.method <- "MaxSens+Spec"}
        if (threshold.method=="prevalence") {threshold.method <- "ObsPrev"}
        if (threshold.method=="equal_sens_spec") {threshold.method <- "Sens=Spec"}
        if (threshold.method=="sensitivity") {threshold.method <- "ReqSens"}
        req.sens <- threshold.sensitivity
        if (threshold.method=="no_omission") {
            threshold.method <- "ReqSens"
            req.sens <- 1.0
        }
        result <- PresenceAbsence::optimal.thresholds(data2, threshold=seq(from=0, to=1, by=0.0001), req.sens=req.sens)
        result2 <- as.numeric(result[, 2])
        names(result2) <- result[, 1]

# threshold2005.min and threshold2005.mean build on results from Liu et al. 2005. Ecography 28: 385-393
# this study of selecting best thresholds recommended following approaches
# 4. prevalence (ObsPrev)
# 5. average probability (MeanProb)
# 7. sensitivity-specificity maximization (MaxSens+Spec)
# 8. sensitivity-specificity equality (Sens=Spec)
# 9. ROC-based (MinROCdist)

        if (threshold.method %in% c("threshold2005.min", "threshold2005.mean")) {
            t1 <- result2[["MaxSens+Spec"]]
            cat(paste("\n", "Threshold (method: MaxSens+Spec): ", sep = ""))
            print(t1)
            t2 <- result2[["Sens=Spec"]]
            cat(paste("Threshold (method: Sens=Spec): ", sep = ""))
            print(t2)         
            t3 <- result2[["ObsPrev"]]
            cat(paste("Threshold (method: ObsPrev): ", sep = ""))
            print(t3)
            t4 <- result2[["MeanProb"]]
            cat(paste("Threshold (method: MeanProb): ", sep = ""))
            print(t4)
            t5 <- result2[["MinROCdist"]]
            cat(paste("Threshold (method: MinROCdist): ", sep = ""))
            print(t5)
            thresholds <- as.numeric(c(t1, t2, t3, t4, t5))
            thresholds <- thresholds[thresholds > 0]
            if (threshold.method == "threshold2005.min") {return(min(thresholds))}
            if (threshold.method == "threshold2005.mean") {return(mean(thresholds))}
        }

# threshold2013.min and threshold2013.mean build on results from Liu et al. 2013. Journal of Biogeography 40: 778-789
# this study of selecting best thresholds recommended following approaches
# 7. maximizing the sum of sensitivity and specificity, max SSS (MaxSens+Spec)
# 1. training data prevalence, trainPrev (ObsPrev)
# 2. mean predicted value for a set of random points over the whole study area, meanPred (MeanProb)
# when species prevalence is high, however, max SSS is the superior method (page 786)

        if (threshold.method %in% c("threshold2013.min", "threshold2013.mean")) {
            t1 <- result2[["MaxSens+Spec"]]
            cat(paste("\n", "Threshold (method: MaxSens+Spec): ", sep = ""))
            print(as.numeric(t1))
            t3 <- result2[["ObsPrev"]]
            cat(paste("Threshold (method: ObsPrev): ", sep = ""))
            print(as.numeric(t3))
            t4 <- result2[["MeanProb"]]
            cat(paste("Threshold (method: MeanProb): ", sep = ""))
            print(as.numeric(t4))
            thresholds <- as.numeric(c(t1, t3, t4))
            thresholds <- thresholds[thresholds > 0]
            if (threshold.method == "threshold2013.min") {return(min(thresholds))}
            if (threshold.method == "threshold2013.mean") {return(mean(thresholds))}
        }
        return(as.numeric(result2[[threshold.method]]))
    }

    if (threshold.PresenceAbsence == F) {
        result <- dismo::threshold(eval, sensitivity=threshold.sensitivity)        
        return(result[[threshold.method]])
    }
}

