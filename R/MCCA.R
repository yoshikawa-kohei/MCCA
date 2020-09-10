
#' Internal control function for cca
#' 
#' list of parameters for controling cca1d fitting
#' 
#' @param max.itr number of maximaum iteration for optimization
#' @param tol hoge
#' @param centering hoge
#' @param verbose hoge
#'
#' @return a list of controling parameter.
mcca.control <- function(max.itr=100, tol=1e-6, centering=TRUE, verbose=FALSE){
  return(list(max.itr=max.itr, tol=tol, centering=centering, verbose=verbose))
}

#' Multilinear Common Component Analysis (MCCA)
#' 
#' @param tnsr tensor object (dim1, dim2, ..., num of sample, num of group)
#' @param ranks hoge
#' @param control a list of internal parameters controlling the model fitting
mcca <- function(tnsr, ranks = NULL, control=list()){
  if (is.null(ranks)) {
    stop("ranks must be specified")
  }
  if (length(tnsr)<2) {
    stop("number of groups must be more than 1")
  }
  for (g in 1:length(tnsr)) {
    stopifnot(is(tnsr[[g]], "Tensor"))
  }
  if (length(ranks) != (tnsr[[1]]@num_modes - 1)) {
    stop("dimension of ranks must be same as (the mode of tensor -1)")
  }
  if (sum(ranks > tnsr[[1]]@modes[1:length(ranks)]) != 0) {
    stop("ranks must be smaller than the corresponding mode")
  }
  if (sum(ranks <= 0) != 0) {
    stop("ranks must be positive")
  }


  # get the parameters controling function of mcca to "control"
  message("#[info] Start : Multilinear Common Component Analysis (MCCA)")
  control <- do.call("mcca.control", control)
  
  num.groups <- length(tnsr)
  num.max.mode <- tnsr[[1]]@num_modes
  
  tnsr.org <- tnsr
  tnsr.mean <- vector("list", num.groups)
  if (control$centering) {
    message("#[info] Centering tensor by group...")
    tnsr.mean <- lapply(tnsr, function(x){return(modeMean(x, m=x@num_modes, drop=TRUE))})
    tnsr <- mapply(function(tnsr1, tnsr2){return(as.tensor(sweep(tnsr1@data, (1:(num.max.mode-1)), tnsr2@data, "-")))}, tnsr, tnsr.mean)
  }
  
  # calculate covariance matrix
  # S.all = [S(group1) , S(group2), ..., S(kronecker product)]
  # S(group1) = [S(mode1), S(mode2), ...]
  S.all <- vector("list", num.groups)
  for (g in 1:num.groups) {
    S.all[[g]] <- calc.tnsr.covmat(tnsr[[g]])
  }

  # store approximate error to "ae.history"
  ae.history <- numeric()
  # loss.history <-  vector("list", num.max.mode-1)
  # initialization ----------------------------------------------------------
  message("#[info] Initializing...")
  # calc weight
  weight <- array(0, c(num.groups, num.max.mode-1))
  for (g in 1:num.groups) {
    for (mode in 1:(num.max.mode-1)) {
      weight[g, mode] <- sum(eigen(tcrossprod(S.all[[g]]$S[[mode]]), symmetric=TRUE)$values[1:ranks[mode]])
    }
  }
  
  # calc weight cov mat
  V <- vector("list", num.max.mode-1)
  
  #weighted.covmat <- vector("list", num.max.mode-1)
  for (mode in 1:(num.max.mode-1)) {
    weighted.covmat <- matrix(0, nrow=tnsr[[1]]@modes[mode], ncol=tnsr[[1]]@modes[mode])
    for (g in 1:num.groups) {
      weighted.covmat <- weighted.covmat + prod(weight[g,-mode]) * tcrossprod(S.all[[g]]$S[[mode]])
    }
    V.eig <- eigen(weighted.covmat, symmetric=TRUE)$vectors
    V[[mode]] <- V.eig[, 1:ranks[mode], drop=FALSE]
  }

  # Lambda step
  Lambda <- vector("list", num.groups)
  for (g in 1:num.groups) {
    Lambda[[g]] <- vector("list", num.max.mode-1)
    for (mode in 1:(num.max.mode-1)) {
      Lambda[[g]][[mode]] <- t(V[[mode]]) %*% S.all[[g]]$S[[mode]] %*% V[[mode]]
    }
  }
  
  # check the convergence property
  # calculate alpha
  alpha <- numeric(num.max.mode-1)
  if(TRUE){
  for (mode in 1:(num.max.mode-1)) {
    weighted.covmat <- matrix(0, nrow=tnsr[[1]]@modes[mode], ncol=tnsr[[1]]@modes[mode])
    for (g in 1:num.groups) {
      weighted.covmat <- weighted.covmat + prod(weight[g,-mode]) * tcrossprod(S.all[[g]]$S[[mode]])
    }
    alpha[mode] <- tr(t(V[[mode]]) %*% weighted.covmat %*% V[[mode]])/tr(weighted.covmat)
  }
  print(alpha)
  }
  
  # Local Optimization ------------------------------------------------------
  pb <- progress_bar$new(total = control$max.itr, format = "#[info] Searching optimal parameters... [:bar] in :elapsed", clear = FALSE, width=100)
  for (itr in 1:control$max.itr) {
    pb$tick()
    # calc weight
    weight <- array(0, c(num.groups, num.max.mode-1))
    for (g in 1:num.groups) {
      for (mode in 1:(num.max.mode-1)) {
        weight[g, mode] <- tr(S.all[[g]]$S[[mode]] %*% V[[mode]] %*%  t(V[[mode]]) %*% S.all[[g]]$S[[mode]])
      }
    }
    
    # update V
    for (mode in 1:(num.max.mode-1)) {
      weighted.covmat <- matrix(0, nrow=tnsr[[1]]@modes[mode], ncol=tnsr[[1]]@modes[mode])
      for (g in 1:num.groups) {
        weighted.covmat <- weighted.covmat + prod(weight[g,-mode]) * tcrossprod(S.all[[g]]$S[[mode]] %*% V[[mode]])
      }
      V.eig <- eigen(weighted.covmat, symmetric=TRUE)$vectors
      V[[mode]] <- V.eig[, 1:ranks[mode], drop=FALSE]
    }

    # Lambda update
    for (g in 1:num.groups) {
      for (mode in 1:(num.max.mode-1)) {
        Lambda[[g]][[mode]] <- t(V[[mode]]) %*% S.all[[g]]$S[[mode]] %*% V[[mode]]
      }
    }
    
    # calculate loss
    ae.history[itr] <- calc.ae.mcca(S.all, V, Lambda)
    
    
    if (control$verbose) {
      message("#[info] MCCA : Approximate Error (AE): ", ae.history[itr])
    }
    
    if (itr != 1 && (abs(ae.history[itr]- ae.history[itr-1])/ae.history[itr-1]) < control$tol  ) {
      pb$tick(control$max.itr)
      break
    }
  }
  

  
  # tnsr.core is the dimensionally reduced tensor
  # tnsr.est is low-rank approximate tensor
  tnsr.core <- lapply(tnsr, function(tnsr){return(ttl(tnsr, lapply(V, t), ms=(1:(num.max.mode-1))))})
  tnsr.est <- lapply(tnsr.core, function(tnsr){return(ttl(tnsr, V, ms=(1:(num.max.mode-1))))})
  
  if(control$centering){
    tnsr.reconst <- mapply(function(tnsr1, tnsr2){return(as.tensor(sweep(tnsr1@data, (1:(num.max.mode-1)), tnsr2@data, "+")))}, tnsr.est, tnsr.mean)
  }else{
    tnsr.reconst <- tnsr.est
  }
  tnsr.reconst.combined <- as.tensor(abind(lapply(tnsr.reconst, function(x){return(x@data)}), along=num.max.mode))

  # calculate original tnsor norm
  tnsr.combined <- as.tensor(abind(lapply(tnsr.org, function(x){return(x@data)}), along=num.max.mode))
  
  tnsr.reconst.combined.fnorm <- fnorm(tnsr.combined - tnsr.reconst.combined)
  tnsr.fnorm <- fnorm(tnsr.combined)
  
  reconst.rate <- ((tnsr.reconst.combined.fnorm^2)/(tnsr.fnorm^2))*100
  
  message("#[info] Finish!")
  # return values
  return(list(core=tnsr.core, est=tnsr.est, tnsr.mean=tnsr.mean, tnsr.reconst=tnsr.reconst, V=V, Lambda=Lambda, ae.history=ae.history, alpha=alpha, reconst.rate=reconst.rate))
  
}


# auxiliary function ------------------------------------------------------

#' matrix trace
#' 
#' @param X hoge
tr <- function(X){
  return(sum(diag(X)))
}

#' sym mat
#' 
#' @param X hoge
sym <- function(X){
  return((X+t(X))/2)
}

#' Calculate Cov
#' 
#' @param tnsr hoge
calc.tnsr.covmat <- function(tnsr){
  
  S <- vector("list", tnsr@num_modes-1)
  for (mode in 1:(tnsr@num_modes-1)) {
    unf.tnsr <- k_unfold(tnsr, mode)
    S[[mode]] <- (unf.tnsr@data %*% t(unf.tnsr@data))/prod(tnsr@modes[-mode])
  }
  
  S.kp <- S[[1]]
  for (mode in 2:(tnsr@num_modes-1)) {
    S.kp <- kronecker(S.kp, S[[mode]])
  }
  
  return(list(S=S, S.kp=S.kp))
}


#' Prot
#' 
#' @param cca hoge
mcca.plot <- function(cca){
  df <- data.frame(
    itr=1:length(cca$ae.history),
    loss=cca$ae.history)
  ggplot(df, aes(x=itr, y=loss)) + geom_line() + theme_bw()
  # plot(cca$ae.history, type="l")
}


#' Calculate Approximate Error 
#' 
#' @param S hoge
calc.ae.mcca <- function(S, V, Lambda){
  numerator <- 0
  denominator <- 0
  
  V.kp <- V[[1]]
  for (mode in 2:length(V)) {
    V.kp <- kronecker(V.kp, V[[mode]])
  }
  
  Lambda.kp <- vector("list", length(Lambda))
  for (g in 1:length(Lambda)) {
    Lambda.kp[[g]] <- Lambda[[g]][[1]]
  }
  for (g in 1:length(Lambda)) {
    for (mode in 2:length(V)) {
      Lambda.kp[[g]] <- kronecker(Lambda.kp[[g]], Lambda[[g]][[mode]])
    }
  }
  
  for(g in 1:length(Lambda)){
    numerator <- numerator + norm(S[[g]]$S.kp - V.kp  %*% Lambda.kp[[g]] %*% t(V.kp), type="F")^2
    denominator <- denominator + norm(S[[g]]$S.kp, type="F")^2
  }
  ARE <- numerator/denominator
  return(ARE)
}
