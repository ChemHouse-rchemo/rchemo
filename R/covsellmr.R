covsellmr <- function(X, Y, nvar = NULL, Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL) {
  
  X <- .mat(X)
  varnames <- colnames(X)
  zdim <- dim(X)
  n <- zdim[1]
  p <- zdim[2]
  Y <- .mat(Y, "x")
  q <- dim(Y)[2]
  
  if(is.null(nvar)) 
    nvar <- p
  
  if(is.null(weights))
    weights <- rep(1 / n, n)
  weights <- .mweights(weights)
  
  covselout <- covsel(X = X, Y = Y, nvar = nvar, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights)
  
  lmrlist <- lapply(1:nvar, function(i) lmr(X[,covselout$sel$sel[1:i], drop = FALSE], Y, weights = weights))
  names(lmrlist) <- paste0("nvar=",1:nvar)
  
  structure(
    list(sel = covselout, fm = lmrlist, weights = weights),
    class = c("Covsellmr"))
}
  

predict.Covsellmr <- function(object, X, ..., nvar = NULL) {
  X <- .mat(X)

  if(is.null(nvar)){
    nvar <- length(object$sel$sel$sel)
  }
  
  le_nvar <- length(nvar)

  pred <- lapply(1:le_nvar, function(i) predict(object$fm[[nvar[i]]], X[,which(colnames(X) %in% object$sel$sel$selnames[1:nvar[i]]),drop=FALSE])$pred)
  names(pred) <- paste("nvar=", nvar, sep = "")

  return(pred)
}
