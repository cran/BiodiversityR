`ordicoeno` <-
function(x,ordiplot,axis=1,...) {
    if (!require(mgcv)) {stop("Requires package mgcv")}
    ordiscore <- scores(ordiplot,display="sites")[,axis]
    original <- cbind(x,ordiscore)
    sorted <- original
    seq <- order(ordiscore)
    sorted[1:nrow(original),] <- original[seq,]
    edfs <- array(NA,dim=c(ncol(x)))
    names(edfs) <- colnames(x)
    palette(rainbow(ncol(x)))
    gamresult <- gam(sorted[,1]~s(ordiscore),data=sorted)
    edfs[1] <- summary(gamresult)$edf
    newdata <- data.frame(seq(min(sorted$ordiscore), max(sorted$ordiscore), length = 1000))
    colnames(newdata) <- "ordiscore"
    gamresult2 <- predict(gamresult,newdata)
    plot(newdata$ordiscore,gamresult2,type="l",ylim=c(0,max(x)),
        col=1,pch=1,xlab="site score on ordination axis",ylab="species values",...)    
    for (i in 2:ncol(x)) {
        gamresult <- gam(sorted[,i]~s(ordiscore),data=sorted)
        gamresult2 <- predict(gamresult,newdata)
        edfs[i] <- summary(gamresult)$edf
        points(newdata$ordiscore,gamresult2,type="l",pch=19,col=i,...)
    }
    palette("default")
    cat("edfs from GAM models for each species...\n")
    return(edfs)
}

