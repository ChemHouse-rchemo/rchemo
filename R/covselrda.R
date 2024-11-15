covselrda <- function(X, y, nvar = NULL, Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL) {
  
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
  
  lmrlist <- lapply(1:nvar, function(i) lmr(X[,covselout$sel$sel[1:i], drop = FALSE], zY$Y, weights = weights))
  names(lmrlist) <- paste0("nvar=",1:nvar)
  
  structure(
    list(sel = covselout, fm = lmrlist, lev = zY$lev, ni = zY$ni, weights = weights),
    class = c("Covselrda"))
}


predict.Covselrda <- function(object, X, ..., nvar = NULL) {
  X <- .mat(X)
  
  if(is.null(nvar)){
    nvar <- length(object$sel$sel$sel)
  }
  le_nvar <- length(nvar)
  
  res <- list()
  res$posterior <- lapply(1:le_nvar, function(i) predict(object$fm[[nvar[i]]], X[,which(colnames(X) %in% object$sel$sel$selnames[1:nvar[i]]),drop=FALSE])$pred)
  res$pred <- lapply(1:length(res$posterior), function(i) matrix(.replace_bylev(apply(res$posterior[[i]], FUN = .findmax, MARGIN = 1), object$lev), ncol = 1))
  names(res$posterior) <- names(res$pred) <- paste("nvar=", nvar, sep = "")
  
  return(res)
}

