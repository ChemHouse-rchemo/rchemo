consensuspca <- function(Xlist, blockscaling = TRUE, weights = NULL, nlv, Xscaling = c("none", "pareto", "sd")[1], algo = c("svd","eigen","eigenk","nipals","nipalsna","sph")[1], gs = TRUE,
                         tol = .Machine$double.eps^0.5, maxit = 200) {
  Xlist <- lapply(1:length(Xlist), function(i) .mat(Xlist[[i]]))

  if(is.null(weights))
    weights <- rep(1, nrow(Xlist[[1]]))
  weights <- .mweights(weights)
  
  xmeanslist <- lapply(1:length(Xlist), function(i) .colmeans(Xlist[[i]], weights = weights))

  if((length(Xscaling)=1) & (length(Xlist)>1)){Xscaling = rep(Xscaling, length(Xlist))}
  
  xscaleslist <- list()
  for(i in 1:length(Xlist)){
    if(Xscaling[i] == "none"){xscaleslist[[i]] <- rep(1, ncol(Xlist[[i]]))}
    if(Xscaling[i] == "pareto"){xscaleslist[[i]] <- sqrt(sqrt(.colvars(Xlist[[i]], weights = weights)))}
    if(Xscaling[i] == "sd"){xscaleslist[[i]] <- sqrt(.colvars(Xlist[[i]], weights = weights))}
  }
  Xlist <- lapply(1:length(Xlist), function(i) scale(Xlist[[i]], center = xmeanslist[[i]], scale = xscaleslist[[i]]))
  
  if(blockscaling==TRUE){
    Xblockscaled <- blockscal(Xtrain = Xlist, weights = weights)
    Xlist <- Xblockscaled$Xtrain
    Xnorms <- Xblockscaled$disp
  }else{
    Xnorms <- NA
  }

  Xconc <- do.call("cbind",Xlist)
  # zdim <- dim(Xconc)
  # n <- zdim[1]
  # p <- zdim[2]

  if(algo=="svd"){respca <- pcasvd(Xconc, weights = NULL, nlv)} 
  if(algo=="eigen"){respca <- pcaeigen(Xconc, weights = NULL, nlv)} 
  if(algo=="eigenk"){respca <- pcaeigenk(Xconc, weights = NULL, nlv)} 
  if(algo=="nipals"){respca <- pcanipals(Xconc, weights = NULL, nlv, gs = gs, tol = tol)} 
  if(algo=="nipalsna"){respca <- pcanipalsna(Xconc, nlv, gs = gs, tol = tol)} 
  if(algo=="sph"){respca <- pcasph(Xconc, weights = NULL, nlv)} 
  
  Tk <- lapply(1:length(Xlist), function(i) Xlist[[i]]%*%t(Xlist[[i]])%*%respca$T)

  structure(
    list(T = respca$T, Tk = Tk, P = respca$P, sv = respca$sv, eig = respca$eig,
         xmeans = xmeanslist, xscales = xscaleslist, weights = weights, blockscaling = blockscaling, Xnorms = Xnorms, niter = NULL, conv = NULL),
    class = c("Consensuspca"))
}


summary.Consensuspca <- function(object, X, ...) {
  
  Xinit <- X
  p <- dim(object$P)[1]
  nlv <- dim(object$T)[2]
  TT <- object$weights * object$T * object$T
  tt <- colSums(TT)
  
  X <- lapply(1:length(X), function(i) scale(.mat(X[[i]]), center = object$xmeans[[i]], scale = object$xscales[[i]]))
  if(object$blockscaling==TRUE){X <- blockscal(Xtrain = X, weights = object$weights)$Xtrain}
  Xconc <- do.call("cbind",X)
  
  sstot <- sum(object$weights * Xconc * Xconc, na.rm = TRUE)
  pvar <- tt / sstot
  cumpvar <- cumsum(pvar)
  explvar <- data.frame(pc = seq(nlv), var = tt, pvar = pvar, cumpvar = cumpvar)
  row.names(explvar) <- seq(nlv)
  
  Tknorm <- object$Tk
  for(j in 1:length(X)){
    for(i in 1:nlv){Tknorm[[j]][,i] <- Tknorm[[j]][,i,drop=FALSE]/norm(Tknorm[[j]][,i,drop=FALSE], type = "e")}
  }
  Tnorm <- object$T
  for(i in 1:nlv){Tnorm[,i] <- object$T[,i,drop=FALSE]/norm(object$T[,i,drop=FALSE], type = "e")}
  contr.block <- matrix(NA, ncol = nlv, nrow = length(X), dimnames = list(paste0("X",1:length(X)),paste0("pc",1:nlv)))
  for(i in 1:nlv){
    contr.block[,i] <- sapply(1:length(X), function(j) nrow(X[[1]])*cov(Tknorm[[j]][,i],Tnorm[,i]))
    contr.block[,i] <- contr.block[,i]/sum(contr.block[,i])*100
  }
  
  contr.ind <- data.frame(.scale(TT, center = rep(0, nlv), scale = tt))
  
  cor.circle <- contr.var <- coord.var <- NULL
  # xvars <- .colvars(X, weights = object$weights) equivalent to object$xscaleslist
  # zX <- .scale(X, center = rep(0, p), scale = sqrt(xvars)) equivalent to Xconc
  zT <- .scale(object$T, center = rep(0, nlv), scale = sqrt(tt))
  cor.circle <- data.frame(t(object$weights * Xconc) %*% zT)#data.frame(t(object$weights * zX) %*% zT)
  coord.var <- data.frame(crossprod(do.call("cbind", Xinit), object$weights * zT))#data.frame(crossprod(X, object$weights * zT))
  z <- coord.var^2
  contr.var <- data.frame(.scale(z, rep(0, nlv), colSums(z)))
  row.names(cor.circle) <- row.names(contr.var) <- row.names(coord.var) <- row.names(object$P)
  
  list(explvar = explvar, contr.block = contr.block, contr.ind = contr.ind, 
       contr.var = contr.var, coord.var = coord.var, cor.circle = cor.circle)    
}

transform.Consensuspca <- function(object, X, ..., nlv = NULL) {
  a <- dim(object$T)[2]
  if(is.null(nlv))
    nlv <- a
  else 
    nlv <- min(nlv, a)
  
  X <- lapply(1:length(X), function(i) scale(.mat(X[[i]]), center = object$xmeans[[i]], scale = object$xscales[[i]]))
  if(object$blockscaling==TRUE){X <- lapply(1:length(X), function(i) X[[i]]/object$Xnorms[i])}
  
  Xconc <- do.call("cbind",X)
  
  T <- Xconc %*% object$P[, seq_len(nlv), drop = FALSE]
  colnames(T) <- paste("pc", seq_len(dim(T)[2]), sep = "")
  T
}
