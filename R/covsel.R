covsel <- function(X, Y, nvar = NULL, Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL) {
  
    X <- .mat(X)
    varnames <- colnames(X)
    zdim <- dim(X)
    n <- zdim[1]
    p <- zdim[2]
    Y <- .mat(Y, "y")
    q <- dim(Y)[2]
  
    if(is.null(nvar)) 
        nvar <- p
  
    if(is.null(weights))
        weights <- rep(1 / n, n)
    weights <- .mweights(weights)
  
    xmeans <- .colmeans(X, weights = weights)
    ymeans <- .colmeans(Y, weights = weights)
    
    if(Xscaling == "none") {xscales <- rep(1, p)}
    if(Yscaling == "none") {yscales <- rep(1, q)}
    if(Xscaling == "sd") {xscales <- sqrt(.colvars(X, weights = weights))}
    if(Yscaling == "sd") {yscales <- sqrt(.colvars(Y, weights = weights))}
    if(Xscaling == "pareto") {xscales <- sqrt(sqrt(.colvars(X, weights = weights)))}
    if(Yscaling == "pareto") {yscales <- sqrt(sqrt(.colvars(Y, weights = weights)))}
    
    X <- scale(X, center = xmeans, scale = xscales)
    Y <- scale(Y, center = ymeans, scale = yscales)
  
    xsstot <- sum(weights * X * X)
    ysstot <- sum(weights * Y * Y)
    yss <- xss <- selvar <- vector(length = nvar)
    for(i in seq_len(nvar)) {
        z <- rowSums(crossprod(weights * X, Y)^2)
        selvar[i] <- which(z == max(z))
        u <- X[, selvar[i], drop = FALSE]
        P <- tcrossprod(u) %*% diag(weights) / sum(weights * u * u)
        # Same as
        #P <- u %*% solve(t(u) %*% D %*% u) %*% t(u) %*% D
        #P <- u %*% t(u) %*% D / sum(d * u * u)
        #P <- crossprod(u %*% t(u), D) / sum(d * u * u)
        ## End
        ## The deflated X and Y are centered matrix (with metric D)
        X <- X - P %*% X 
        Y <- Y - P %*% Y
        xss[i] <- sum(weights * X * X)
        yss[i] <- sum(weights * Y * Y)
        }
    cumpvarx <- 1 - xss / xsstot
    cumpvary <- 1 - yss / ysstot
    sel <- data.frame(sel = selvar, 
                      cumpvarx = cumpvarx, cumpvary = cumpvary,
                      selnames = sapply(selvar, function(j) varnames[j]))
  
    structure(
      list(sel = sel, weights = weights),# xmeans = xmeans, ymeans = ymeans, xscales = xscales, yscales = yscales)
      class = c("Covsel"))
    }


