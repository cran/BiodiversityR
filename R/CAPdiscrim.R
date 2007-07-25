`CAPdiscrim` <-
function(formula,data,dist="bray",axes=4,m=0,permutations=0) {
    if (!require(MASS)) {stop("Requires package MASS")}
    CAPresult=function(points,y,group,axes=4,m=1,eig) {
        lda1 <- lda(y[,group]~points[,1:m],CV=T)
        lda2 <- lda(y[,group]~points[,1:m])
        matches <- (lda1$class == y[,group])
        correct <- sum(matches) / length(matches) * 100
        lda3 <- predict(lda2,y[,group])
        rownames(lda3$x) <- rownames(points)
        tot <- sum(eig)
        varm <- sum(eig[1:m])/tot*100
        result <- list(PCoA=points[,1:axes],m=m,tot=tot,varm=varm,group=y[,group],CV=lda1$class,percent=correct,
            x=lda3$x,F=(lda2$svd)^2)
        return(result)
    }
    x <- eval(as.name((all.vars(formula)[1])))
    group <- all.vars(formula)[2]
    y <- data
    if (inherits(x, "dist")) {
        distmatrix <- x
        x <- data.frame(as.matrix(x))
    }else{
        distmatrix <- vegdist(x, method = dist)
    }
    pcoa <- cmdscale(distmatrix, k=nrow(x)-1, eig=T, add=F)    
    points <- pcoa$points
    rownames(points) <- rownames(x)
    eig <- pcoa$eig
    if (m==0) {
        tot <- sum(eig)
        correct <- -1
        for (i in 1:(nrow(x)-1)) {
            if (sum(eig[1:i]) < tot) {
                result1 <- CAPresult(points=points,y=y,group=group,axes=axes,m=i,eig=eig)
                if (result1$percent > correct) {
                    correct <- result1$percent
                    result <- result1
                } 
            }            
        }
    }else{
        result <- CAPresult(points=points,y=y,group=group,axes=axes,m=m,eig=eig)
    }
    if (permutations>0) {
        permutations <- permutations-1
        permresult <- numeric(permutations)
        y1 <- y
        n <- nrow(y1)
        for (j in 1: permutations){
            y1[1:n,] <- y[sample(n),]
            if (m==0) {
                correct <- -1
                for (i in 1:(nrow(x)-1)) {
                    if (sum(eig[1:i]) < tot) {
                        result1 <- CAPresult(points=points,y=y1,group=group,axes=axes,m=i,eig=eig)
                        if (result1$percent > correct) {
                            correct <- result1$percent
                            resultp <- result1
                        } 
                    }            
                }
            }else{
                resultp <- CAPresult(points=points,y=y1,group=group,axes=axes,m=m,eig=eig)
            }
            permresult[j] <- resultp$percent
        }
        signi <- sum(permresult > result$percent)
        signi <- (1+signi)/(1+permutations)
        cat("Percentage of correct classifications was", result$percent, "\n") 
        cat("Significance of this percentage was", signi, "\n\n")
        permresult <- summary(permresult)
    }else{
        signi <- NA
        permresult <- NULL
    }
    m <- result$m
    if (m>1) {
        result1 <- summary(manova(points[,1:m]~y[,group]))
    }else{
        result1 <- summary(lm(points[,1]~y[,group]))
    }
    result2 <- list(PCoA=result$PCoA,m=m,tot=result$tot,varm=result$varm,group=result$group,CV=result$CV,
        percent=result$percent,x=result$x,F=result$F,manova=result1,signi=signi,permutations=permresult)        
    return(result2)
}


