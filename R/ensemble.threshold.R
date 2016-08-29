`ensemble.threshold` <- function(
    eval, threshold.method="spec_sens", threshold.sensitivity=0.9, 
    threshold.PresenceAbsence=FALSE, Pres, Abs) 
{
    if (threshold.PresenceAbsence == T){        
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
        result <- PresenceAbsence::optimal.thresholds(data2, threshold=seq(from=0, to=1, by=0.005), req.sens=req.sens)
        result2 <- as.numeric(result[, 2])
        names(result2) <- result[, 1]

# threshold.min and threshold.mean build on results from Liu et al. 2005. Ecography 28: 385-393
# this study of selecting best thresholds recommended following approaches
# 4. prevalence
# 5. average probability
# 7. sensitivity-specificity maximization 
# 8. sensitivity-specificity equality
# 9. ROC-based 
# approaches 4 (prevalence), 7 (spec_sens) and 8 (equal_spec_sens) are available from dismo::threshold
# all approaches 4 (ObsPrev), 5 (MeanProb), 7 (MaxSens+Spec), 8 (Sens=Spec) and 9 (MinROCdist) are available from PresenceAbsence::optimal.thresholds

        if (threshold.method == "threshold.min") {
            t1 <- result2[["MaxSens+Spec"]]
            t2 <- result2[["Sens=Spec"]]            
            t3 <- result2[["ObsPrev"]]
            t4 <- result2[["MeanProb"]]
            t5 <- result2[["MinROCdist"]]
            thresholds <- as.numeric(c(t1, t2, t3, t4, t5))
            thresholds <- thresholds[thresholds > 0]
            return(min(thresholds))
        }
        if (threshold.method == "threshold.mean") {
            t1 <- result2[["MaxSens+Spec"]]
            t2 <- result2[["Sens=Spec"]]            
            t3 <- result2[["ObsPrev"]]
            t4 <- result2[["MeanProb"]]
            t5 <- result2[["MinROCdist"]]
            thresholds <- as.numeric(c(t1, t2, t3, t4, t5))
            thresholds <- thresholds[thresholds > 0]
            return(mean(thresholds))
            }
        return(as.numeric(result2[[threshold.method]]))
    }else{
        result <- dismo::threshold(eval, sensitivity=threshold.sensitivity)        
        if (threshold.method == "threshold.min") {
            t1 <- result[["spec_sens"]]
            t2 <- result[["equal_sens_spec"]]            
            t3 <- result[["prevalence"]]
            thresholds <- as.numeric(c(t1, t2, t3))
            thresholds <- thresholds[thresholds > 0]
            return(min(thresholds))
        }
        if (threshold.method == "threshold.mean") {
            t1 <- result[["spec_sens"]]
            t2 <- result[["equal_sens_spec"]]            
            t3 <- result[["prevalence"]]
            thresholds <- as.numeric(c(t1, t2, t3))
            thresholds <- thresholds[thresholds > 0]
            return(mean(thresholds))
        }
        return(result[[threshold.method]])
    }
}

