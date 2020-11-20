
#' Internal control function for cca1d
#' 
#' list of parameters for controling cca1d fitting
#' 
#' @param max.itr number of maximaum iteration for optimization
#' 
#' @return a list of controling parameter.
cca1d.control <- function(max.itr=100, tol=1e-6, centering=TRUE, verbose=FALSE){
  return(list(max.itr=max.itr, tol=tol, centering=centering, verbose=verbose))
}


#' Common Component Analysis (CCA)
#' 
#' WIP
#' 
#' @param control a list of internal parameters controlling the model fitting
cca1d <- function(X, rank, control=list()){
  control <- do.call("cca1d.control", control)
  
  cat("#[info] Start : CCA for vectorized data \n")
  n <- dim(X)[1]
  p <- dim(X)[2]
  num.group <- dim(X)[3]
  
  X.org <- X
  X.mean <- array(0, c(p,num.group))
  if(control$centering){
    cat("#[info] Centering... \n")

    for(g in 1:num.group){
      X.mean[,g] <- apply(X[,,g], 2, mean)
      for(i in 1:n){
        X[i,,g] <- X[i,,g] - X.mean[,g]
      }
    }
  }
  
  # calculate covariance matrix
  S <- array(0, c(p, p, num.group))
  for(g in 1:num.group){
    S[,,g] <- (t(X[,,g]) %*% X[,,g])/n
  }
  
  # initialization ----------------------------------------------------------
  M <- matrix(0, nrow=p, ncol=p)
  for(g in 1:num.group){
    M <- M + S[,,g]
  }  
  V.full <- eigen(M, symmetric=TRUE)$vectors
  V <- V.full[, 1:rank, drop=FALSE]  
  
  ae.histroy <- numeric()
  # Local Optimization ------------------------------------------------------
  for(itr in 1:control$max.itr){
    # V step
    M <- matrix(0, nrow=p, ncol=p)
    for(g in 1:num.group){
      M <- M + tcrossprod(S[,,g] %*% V)
    }  
    V.full <- eigen(M, symmetric=TRUE)$vectors
    V <- V.full[, 1:rank, drop=FALSE] 
    
    # latent covariance matrix
    Lambda <- array(0, c(rank, rank, num.group))
    for(g in 1:num.group){
      Lambda[,,g] <- t(V) %*% S[,,g] %*% V
    } 
    
    # calculate approximate error
    ae.histroy[itr] <- calc.AE.cca1d(S, V, Lambda)
    if(control$verbose){
      cat("#[info] Approximate Error (AE): ", ae.histroy[itr], "\n")
    }
    if(itr != 1 && (abs(ae.histroy[itr] - ae.histroy[itr-1])/ae.histroy[itr-1]) < control$tol ){
      break
    }
  }
  
  X.est <- array(0, c(n, p, num.group))
  for (g in 1:num.group) {
    X.est[,,g] <- sweep(X[,,g] %*% V %*% t(V), 2, X.mean[,g], FUN="+")
  }
  
  cat("#[info] Finish : CCA for vectorized data \n")
  return(list(X=X, X.est=X.est, X.mean=X.mean, V = V, Lambda=Lambda, ae.histroy=ae.histroy))
}

cca1d.plot <- function(cca){
  plot(cca$ae.histroy, type="l")
}

calc.AE.cca1d <- function(S, V, Lambda){
  numerator <- 0
  denominator <- 0
  num.group <- dim(S)[3]
  for(g in 1:num.group){
    numerator <- numerator + norm( S[,,g] - V %*% Lambda[,,g] %*% t(V), type="F")^2
    denominator <- denominator + norm(S[,,g], type="F")^2
  }
  AE <- numerator/denominator
  return(AE)
}



