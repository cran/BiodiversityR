`stackcommunitydataset` <-
    function(comm, remove.zeroes=FALSE, order.sites=FALSE, order.species=FALSE)
{
    x <- data.frame(comm)
    site.names <- rownames(x)
    site.names <- as.character(levels(factor(site.names, levels=unique(site.names))))
    if (order.sites == T) {site.names <- as.character(levels(factor(site.names)))}
    n.sites <- length(site.names)
    species.names <- names(x)
    species.names <- as.character(levels(factor(species.names, levels=unique(species.names))))
    if (order.species == T) {    species.names <- as.character(levels(factor(species.names)))}
    n.species <- length(species.names)
    result <- data.frame(array(0, dim=c(n.sites*n.species, 3)))
    names(result) <- c("sites", "species", "abundance")
    result[, "sites"] <- rep(site.names, each=n.species)
    result[, "species"] <- rep(species.names, times=n.sites)
    for (r in 1:nrow(x)) {
        for (c in 1:ncol(x)) {
            i <- result[, "sites"] == site.names[r] & result[, "species"] == species.names[c]
            if (is.na(x[r, c]) == F) {result[i, "abundance"] <- result[i, "abundance"] + x[r, c]}
         }
    }
    if (remove.zeroes == T) {result <- result[which(result[, "abundance"] > 0), ]}
    return(result)          
}

