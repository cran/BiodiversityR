`bioenv.numeric` <-
function(comm, env, method = "spearman", index = "bray", as.numeric = NULL,
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
    if (length(names(env)) < 2) {warning("Not enough numeric variables in environmental data set")}
    result <- bioenv(comm=comm, env=env, method=method, index=index,
        upto=upto, trace=trace, partial=partial, ...)
    return(result)
}

