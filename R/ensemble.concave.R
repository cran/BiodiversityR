`ensemble.concave.hull` <- function(
        baseline.data,
        change.data,
        complete.cases=TRUE,
        VIF=TRUE, VIF.max=20, VIF.silent=TRUE,
        method=c("rda", "pca", "prcomp"),
        ax1=1, ax2=2,
        concavity=2.5,
        buffer.dist=NA,
        ggplot=TRUE
    )
    {
        
        if (complete.cases==TRUE) {
            baseline.data <- baseline.data[complete.cases(baseline.data), ]
            change.data <- change.data[complete.cases(change.data), ]   
        }
        
        comm1 <- rbind(baseline.data, change.data)
        
        if (VIF==TRUE) {
            VIF.vars <- BiodiversityR::ensemble.VIF.dataframe(baseline.data, 
                                                              VIF.max=VIF.max, 
                                                              car=FALSE,
                                                              silent=VIF.silent)$vars.included
            cat(paste("VIF selection:", "\n"))
            print(VIF.vars)
            comm1 <- comm1[, VIF.vars]
        }
        
        env1 <- data.frame(climate=c(rep("baseline", nrow(baseline.data)),
                                     rep("change", nrow(change.data))))
        env1$climate <- factor(env1$climate)
        
        method <- method[1]
        
        if (method == "rda") { 
            rda1 <- vegan::rda(comm1 ~ climate, data=env1, scale=TRUE)
            envfit.result <- envfit(rda1, env=comm1, choices=c(ax1, ax2))
            envfit.df <- data.frame(variable = names(envfit.result$vectors$r),
                                    R2 = envfit.result$vectors$r, 
                                    pvals =envfit.result$vectors$pvals)
            envfit.df <- envfit.df[envfit.df$pvals <= 0.05, ]
            envfit.df <- envfit.df[envfit.df$R2 >= 0.50, ]
            envfit.vars <- envfit.df$variable
            cat(paste("Variables selected via envfit:", "\n"))
            print(envfit.vars)
            comm2 <- comm1[, which(names(comm1) %in% envfit.vars)]
            rda1 <- vegan::rda(comm2, scale=TRUE)
            rda1 <- BiodiversityR::caprescale(rda1)
        }
        if (method == "pca") { 
            rda1 <- vegan::rda(comm1, scale=TRUE)
            rda1 <- BiodiversityR::caprescale(rda1)
        }
        if (method == "prcomp") { 
            rda1 <- stats::prcomp(comm1, center=TRUE, scale.=TRUE)
        }
        
        rda1.s <- scores(rda1, choices=c(ax1, ax2),
                         scaling="sites", display="sites")
        rda1.b <- as.matrix(rda1.s[env1$climate == "baseline", ])
        rda1.c <- as.matrix(rda1.s[env1$climate == "change", ])
        
        baseline.hull <- concaveman::concaveman(rda1.b, concavity=concavity)
        baseline.hull <- sf::st_make_valid(sf::st_polygon(list(baseline.hull)))
        baseline.area <- sf::st_area(baseline.hull)
        
        change.hull <- concaveman::concaveman(rda1.c, concavity=concavity)
        change.hull <- sf::st_make_valid(sf::st_polygon(list(change.hull)))
        change.area <- sf::st_area(change.hull)
        
        novel.hull <- sf::st_difference(change.hull, baseline.hull)
        novel.area <- sf::st_area(novel.hull)
        
        rda1.base <- data.frame(rda1.b)
        names(rda1.base) <- c("X", "Y")  
        
        rda1.df <- data.frame(rda1.c)
        names(rda1.df) <- c("X", "Y")
        rda1.sf <- sf::st_as_sf(rda1.df, coords=c("X", "Y")) 
        
        intersect.hull <- sf::st_intersection(change.hull, baseline.hull)
        intersect.area <- sf::st_area(intersect.hull)
        
        if (is.na(buffer.dist) == TRUE) {
            range.x <- (max(rda1.df$X)-min(rda1.df$X))/10000
            range.y <- (max(rda1.df$Y)-min(rda1.df$Y))/10000
            buffer.dist <- max(c(range.x, range.y))
        }
        
        intersect.buffer <- sf::st_buffer(intersect.hull, 
                                      dist=buffer.dist)
        
        rda1.novel <- sf::st_intersects(rda1.sf, intersect.buffer, sparse=FALSE)
        #  rda1.novel <- st_contains_properly(novel.hull, rda1.sf, sparse=FALSE)
        
        rda1.df$novel <- as.numeric(rda1.novel) == 0
        
        out1 <- list(rda.object=rda1,
                     method=method,
                     baseline.hull=baseline.hull,
                     baseline.area=baseline.area,
                     change.hull=change.hull,
                     change.area=change.area,
                     overlap.hull=intersect.hull,
                     overlap.area=intersect.area,
                     novel.hull=novel.hull,
                     novel.area=novel.area,
                     buffer.dist=buffer.dist,
                     change.points=rda1.df,
                     baseline.points=rda1.base)
        
        if (ggplot==TRUE) {
            
            ggplot.out <- ggplot() +
                ggplot2::geom_sf(data=baseline.hull,
                        colour="green", fill=ggplot2::alpha("green", 0.4), size=0.7) +
                ggplot2::geom_sf(data=change.hull,
                        colour="blue", fill=ggplot2::alpha("blue", 0.4), size=0.7) +
                ggplot2::geom_sf(data=novel.hull,
                        colour="orange", fill=NA, size=0.7)
            
            if (nrow(rda1.df[rda1.df$novel == FALSE, ]) > 0) {
                points.in <- sf::st_as_sf(rda1.df[rda1.df$novel == FALSE, ],
                                      coords=c("X", "Y"))
                ggplot.out <- ggplot.out +
                    ggplot2::geom_sf(data=points.in,
                            colour="black", size=1.0)
            }else{
                cat(paste("Note: No future points inside the baseline hull", "\n"))
            }
            
            if (nrow(rda1.df[rda1.df$novel == TRUE, ]) > 0) {
                points.out <- sf::st_as_sf(rda1.df[rda1.df$novel == TRUE, ],
                                       coords=c("X", "Y"))
                ggplot.out <- ggplot.out +
                    ggplot2::geom_sf(data=points.out,
                            colour="red", size=1.0)
            }else{
                cat(paste("Note: No future points outside the baseline hull", "\n"))
            }    
            out1 <- as.list(c(out1, list(ggplot.out=ggplot.out)))
        }
        
        return(out1)
    }

`ensemble.concave.venn` <- function(
    x,
    candidate.data,
    concavity=x$concavity,
    buffer.dist=x$buffer.dist,
    ggplot=TRUE,
    show.candidate.points=TRUE
)
{
    
    if (x$method == "rda" || x$method == "pca"){
        rda2.c2 <- predict(x$rda.object, newdata=candidate.data, 
                           scaling="sites", type="wa", model="CA")    
        rda2.c <- rda2.c2[, c(1:2)]
    }
    if (x$method == "prcomp"){
        rda2.c1 <- predict(x$rda.object, newdata=candidate.data, 
                           center=TRUE, .scale=TRUE)    
        rda2.c <- scores(rda2.c1[, c(1:2), drop=FALSE])
    }
    
    rda2.c <- rda2.c[complete.cases(rda2.c), ]
    rda2.c <- as.matrix(rda2.c)
    
    cand.hull <- concaveman::concaveman(rda2.c, concavity=concavity)
    cand.hull <- sf::st_make_valid(sf::st_polygon(list(cand.hull)))
    cand.area <- sf::st_area(cand.hull)
    
    cand.buffer <- sf::st_buffer(cand.hull, 
                             dist=buffer.dist)
    
    rda1.df <- x$change.points
    rda1.sf <- sf::st_as_sf(rda1.df, coords=c("X", "Y")) 
    
    rda1.candidate.in <- sf::st_intersects(rda1.sf, cand.buffer, sparse=FALSE)
    #  rda1.novel <- st_contains_properly(novel.hull, rda1.sf, sparse=FALSE)
    
    rda1.df$candidate.in <- as.numeric(rda1.candidate.in) == 1 
    
    rda2.df <- data.frame(rda2.c)
    names(rda2.df) <- c("X", "Y")
    
    out1 <- list(change.points=rda1.df,
                 candidate.points=rda2.df, 
                 candidate.hull=cand.hull,
                 candidate.area=cand.area)
    
    if (ggplot==TRUE) {  
        baseline.hull <- x$baseline.hull
        change.hull <- x$change.hull
        ggplot.out <- ggplot() +
            ggplot2::geom_sf(data=baseline.hull,
                    colour="green", fill=ggplot2::alpha("green", 0.4), size=0.7) +
            ggplot2::geom_sf(data=change.hull,
                    colour="blue", fill=ggplot2::alpha("blue", 0.4), size=0.7) +
            
            ggplot2::geom_sf(data=cand.hull,
                    colour="orange", fill=NA, size=0.7)
        
        rda1.fut <- rda1.df[rda1.df$novel == TRUE, ] 
        
        if (nrow(rda1.fut[rda1.fut$candidate.in == TRUE, ]) > 0) {
            points.in <- sf::st_as_sf(rda1.fut[rda1.fut$candidate.in == TRUE, ],
                                  coords=c("X", "Y"))
            ggplot.out <- ggplot.out +
                ggplot2::geom_sf(data=points.in,
                        colour="black", size=1.0)
        }else{
            cat(paste("Note: No future points inside the candidate hull", "\n"))
        }
        
        if (nrow(rda1.fut[rda1.fut$candidate.in == FALSE, ]) > 0) {
            points.out <- sf::st_as_sf(rda1.fut[rda1.fut$candidate.in == FALSE, ],
                                   coords=c("X", "Y"))
            ggplot.out <- ggplot.out +
                ggplot2::geom_sf(data=points.out,
                        colour="red", size=1.0)
        }else{
            cat(paste("Note: No future points outside the candidate hull", "\n"))
        }      
        
        if (show.candidate.points == TRUE) {
            rda2.sf <- sf::st_as_sf(rda2.df, coords=c("X", "Y"))
            ggplot.out <- ggplot.out + 
                ggplot2::geom_sf(data=rda2.sf,
                        colour="orange", size=1.0) 
        }   
        
        out1 <- as.list(c(out1, list(ggplot.out=ggplot.out)))
    }  
    
    return(out1)  
    
}

`ensemble.concave.union` <- function(
    x,
    candidate.venns,
    buffer.dist=x$buffer.dist,
    ggplot=TRUE,
    show.candidate.points=TRUE
)
{
    cand.hulls <- vector("list", length(candidate.venns))
    for (i in 1:length(candidate.venns)) {
        cand.hulls[[i]] <- candidate.venns[[i]]$candidate.hull
    }
    cand.hulls <- sf::st_make_valid(sf::st_sfc(cand.hulls))
    cand.hull <- sf::st_make_valid(sf::st_union(cand.hulls, by_feature=FALSE))
    cand.area <- sf::st_area(cand.hull)
    
    cand.buffer <- sf::st_buffer(cand.hull, 
                             dist=buffer.dist)
    
    rda1.df <- x$change.points
    rda1.sf <- sf::st_as_sf(rda1.df, coords=c("X", "Y")) 
    
    rda1.candidate.in <- sf::st_intersects(rda1.sf, cand.buffer, sparse=FALSE)
    
    rda1.df$candidate.in <- as.numeric(rda1.candidate.in) == 1 
    
    out1 <- list(change.points=rda1.df,
                 candidate.hull=cand.hull,
                 candidate.area=cand.area)
    
    if (ggplot==TRUE) {  
        baseline.hull <- x$baseline.hull
        change.hull <- x$change.hull
        ggplot.out <- ggplot() +
            ggplot2::geom_sf(data=baseline.hull,
                    colour="green", fill=ggplot2::alpha("green", 0.4)) +
            ggplot2::geom_sf(data=change.hull,
                    colour="blue", fill=ggplot2::alpha("blue", 0.4)) +
            
            ggplot2::geom_sf(data=cand.hull,
                    colour="orange", fill=NA, size=0.7)
        
        rda1.fut <- rda1.df[rda1.df$novel == TRUE, ] 
        
        if (nrow(rda1.fut[rda1.fut$candidate.in == TRUE, ]) > 0) {
            points.in <- sf::st_as_sf(rda1.fut[rda1.fut$candidate.in == TRUE, ],
                                  coords=c("X", "Y"))
            ggplot.out <- ggplot.out +
                ggplot2::geom_sf(data=points.in,
                        colour="black", size=1.0)
        }else{
            cat(paste("Note: No future points inside the candidate hull", "\n"))
        }
        
        if (nrow(rda1.fut[rda1.fut$candidate.in == FALSE, ]) > 0) {
            points.out <- sf::st_as_sf(rda1.fut[rda1.fut$candidate.in == FALSE, ],
                                   coords=c("X", "Y"))
            ggplot.out <- ggplot.out +
                ggplot2::geom_sf(data=points.out,
                        colour="red", size=1.0)
        }else{
            cat(paste("Note: No future points outside the candidate hull", "\n"))
        }      
        
        if (show.candidate.points == TRUE) {
            rda2.df <- candidate.venns[[1]]$candidate.points
            for (i in 2:length(candidate.venns)) {
                rda2.df <- rbind(rda2.df, candidate.venns[[i]]$candidate.points)
            }
            rda2.sf <- sf::st_as_sf(rda2.df, coords=c("X", "Y"))
            ggplot.out <- ggplot.out + 
                ggplot2::geom_sf(data=rda2.sf,
                        colour="orange", size=1.0) 
        }
        
        out1 <- as.list(c(out1, list(ggplot.out=ggplot.out)))
    }  
    
    return(out1)  
    
}

ensemble.outliers <- function(x,
    ID.var=NULL, bioc.vars=NULL,
    fence.k=2.5, n_min=5) 
{
    y <- x
    y$count <- rep(0, nrow(y))
  
    if (is.null(ID.var)) {
        y$ID <- row.names(y)
        ID.var <- "ID"
    }else{
        x <- x[, which(names(x) != ID.var), drop=FALSE]
    }

    if (is.null(bioc.vars) == FALSE) {x <- x[, bioc.vars, drop=FALSE]}

    for (i in 1:ncol(x)) {    
        iqr <- IQR(x[, i], na.rm=TRUE)
        lb <- summary(x[, i], na.rm=TRUE)[2] - iqr * fence.k
        ub <- summary(x[, i], na.rm=TRUE)[5] + iqr * fence.k
# error in the function when only one variable  
#  outl <- as.numeric(x[, i] <= lb | x[, i] >= ub)
        outl <- as.numeric(x[, i] < lb | x[, i] > ub)   
        y$count <- y$count + outl
    }

    y$outlier <- y$count >= n_min
  
  return(y[, c(ID.var, "count", "outlier")])
}
