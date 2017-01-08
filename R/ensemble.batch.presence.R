`ensemble.batch.presence` <- function(
    p=NULL, species.name="Species001"
)
{
    p <- as.matrix(p)
    p <- cbind(rep(species.name, nrow(p)), p)
    p <- data.frame(p)
    names(p) <- c("species", "x", "y")
    return(p)
}

