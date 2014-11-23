`bioenv.numeric` <-
    function(...)
{
    UseMethod("bioenv.numeric")
}


`bioenv.numeric.default` <-
function(comm, env, method = "spearman", index = "bray", as.numeric = c(),
       upto = ncol(env), trace = FALSE, partial = NULL, ...) {
    env2 <- env
    if (is.null(as.numeric) == F) {
        for (i in 1:length(as.numeric)) {
            if(any(names(env) == as.numeric[i])) { 
                env2[, as.numeric[i]] <- as.numeric(env[, as.numeric[i]])
            }
        }
    }
    vars <- names(env2)
    for (i in 1:length(vars)) {
        focal.var <- which(names(env2)==vars[i])
        if (is.numeric(env2[, focal.var]) == F) {env2 <- env2[, -focal.var]}
    }
    env <- env2
    result <- bioenv(comm=comm, env=env, method=method, index=index,
        upto=upto, trace=trace, partial=partial, ...)
    return(result)
}

`summary.bioenv.numeric` <-
    function(object, ...)
{
    x <- object$models
    nam <- object$names
    size <- seq(1:length(x))
    cor <- unlist(lapply(x, function(tmp) tmp$est))
    pars <- unlist(lapply(x, function(tmp) paste(nam[tmp$best], collapse=" ")))
    out <- list(size = size, correlation = cor, variables = pars)
    class(out) <- "summary.bioenv"
    out
}
