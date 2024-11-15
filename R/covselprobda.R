.covselprobda <- function(X, y, nvar = NULL, fun, prior = c("unif", "prop"), Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL) {
  
  X <- .mat(X)
  varnames <- colnames(X)
  zdim <- dim(X)
  n <- zdim[1]
  p <- zdim[2]
  if(is.factor(y))
    y <- as.character(y)
  zY <- dummy(y)
  q <- dim(zY$Y)[2]
  
  if(is.null(nvar)) 
    nvar <- p
  
  if(is.null(weights))
    weights <- rep(1 / n, n)
  weights <- .mweights(weights)
  
  covselout <- covsel(X = X, Y = zY$Y, nvar = nvar, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights)
  
  fmlist <- lapply(1:nvar, function(i) fun(X[,covselout$sel$sel[1:i], drop = FALSE], y, prior = prior))
  names(fmlist) <- paste0("nvar=",1:nvar)
  
  structure(
    list(sel = covselout, fm = fmlist, lev = zY$lev, ni = zY$ni, weights = weights),
    class = c("Covselprobda"))
}

covsellda <- function(X, y, nvar = NULL, prior = c("unif", "prop"), Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL) {
  .covselprobda(X, y, nvar = NULL, fun = rchemo::lda, prior = c("unif", "prop"), Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL)
}
  
covselqda <- function(X, y, nvar = NULL, prior = c("unif", "prop"), Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL) {
  .covselprobda(X, y, nvar = NULL, fun = rchemo::qda, prior = c("unif", "prop"), Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL)
}

predict.Covselprobda <- function(object, X, ..., nvar = NULL) {
  X <- .mat(X)
  
  if(is.null(nvar)){
    nvar <- length(object$sel$sel$sel)
  }
  le_nvar <- length(nvar)
  
  predlist <- lapply(1:le_nvar, function(i) predict(object$fm[[nvar[i]]], X[,which(colnames(X) %in% object$sel$sel$selnames[1:nvar[i]]),drop=FALSE]))
  res <- list()
  res$posterior <- lapply(1:length(predlist), function(i) predlist[[i]]$posterior)
  res$pred <- lapply(1:length(predlist), function(i) predlist[[i]]$pred)
  names(res$posterior) <- names(res$pred) <- paste("nvar=", nvar, sep = "")
  
  return(res)
}

  