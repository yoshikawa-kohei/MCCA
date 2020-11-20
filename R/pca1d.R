#' Internal control function for cca1d
#' 
#' list of parameters for controling cca1d fitting
#' 
#' @param max.itr number of maximaum iteration for optimization
#' 
#' @return a list of controling parameter.
pca1d.control <- function(max.itr=100, centering=TRUE, verbose=FALSE){
  return(list(max.itr=max.itr, centering=centering, verbose=verbose))
}

#' Principal Component Analysis (PCA)
#' 
#' WIP
#' 
#' @param control a list of internal parameters controlling the model fitting
#' 
#' @return 
pca1d <- function(X, rank, control=list()){
  control <- do.call("pca1d.control", control)
  cat("#[info] PCA for vectorized data \n")
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  X.mean <- apply(X, 2, mean)
  if(control$centering){
    cat("#[info] Centering... \n")
      X.mean <- apply(X, 2, mean)
      for(i in 1:n){
        X[i, ] <- X[i, ] - X.mean
      }
    }

  # calculate covariance matrix
  S <- (t(X) %*% X)/n

# Optimization ------------------------------------------------------------
  V.full <- eigen(S, symmetric=TRUE)
  V <- V.full$vectors[, 1:rank, drop=FALSE]  
  
  # latent covariance matrix
  Lambda <- diag(V.full$values[1:rank, drop=FALSE], nrow=rank, ncol=rank)
  
  # calculate approximate error
  AE <- calc.AE.pca1d(S, V, Lambda)
  if(control$verbose){
    cat("#[info] Approximate Error (AE): ", AE, "\n")
  }
  
  X.est <- sweep(X %*% V %*% t(V), 2, X.mean,  FUN="+")

  cat("#[info] Finish : PCA for vectorized data \n")
  return(list(X.est=X.est, X.mean=X.mean, V = V, Lambda=Lambda, ae.histroy=AE))
}

#'Calculate approximate error of PCA
#' 
#' WIP
#' 
#' @param S
#' @param V
#' @param Lambda
#' 
#' @return approximate error
calc.AE.pca1d <- function(S, V, Lambda){
  AE <-  (norm( S - V %*% Lambda %*% t(V), type="F")^2)/(norm(S, type="F")^2)
  return(AE)
}
