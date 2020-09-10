generate.sim.data <- function(n, p1, p2, num.group){
  # covariance matrixes
  Sigma.col <- 0.5^abs(outer(1:p2, 1:p2,FUN="-"))
  Sigma.row <- 0.5^abs(outer(1:p1, 1:p1,FUN="-"))
  
  # Common matrix
  #Z <- rmatrixnorm(matrix(0,p1,p2), 100*diag(1:p1), 100*diag(1:p2))
  Z <- rmatrixnorm(matrix(0, p1, p2), Sigma.row, Sigma.col)
  
  # Group Common matrix
  W <- list()
  for(g in 1:num.group){
    #W[[g]] <- rmatrixnorm(matrix(g, p1,p2), diag(p1), diag(p2))
    W[[g]] <- rmatrixnorm(matrix(g, p1, p2), diag(p1), diag(p2))
  }
  
  # generate data
  X <- array(0, c(p1, p2, n, num.group))
  for(g in 1:num.group){
    for(i in 1:n){
      # X[row, colmun, sample]
      X[,,i,g] <- W[[g]] + Z + rmatrixnorm(matrix(0, p1, p2), Sigma.row, Sigma.col)
    }
  }
  
  # store data to list
  data = list()
  for(g in 1:num.group){
    data[[g]] <- as.tensor(X[,,,g])
  }
  return(data)
}
