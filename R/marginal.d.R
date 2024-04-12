
#' @title  Marginal predictive density 
#' @description It calculates the marginal density (Lewis and Raftery, 1997) from the DIFM sample using R.
#'
#' @param data  The dataset
#' @param model.attributes  Model attributes generated from \code{difm.model.attributes}.
#' @param hyp.parm  Hyperparameters generated from \code{difm.hyp.parm}.
#' @param Gibbs  Result of Gibbs sampler from DIFM function.
#' @param burnin  Burn-in period. If not specified, one tenths of the iterations will be the burn-in period.
#' @param verbose  Print out the iteration process.
#'
#' @importFrom LaplacesDemon rinvwishart rinvgamma
#'
#' @return  Metropolis-Laplace estimator of the Marginal density 
#' @export

marginal.d <- function(data, model.attributes, hyp.parm, Gibbs, burnin = NA, verbose = TRUE){
  n <- nrow(data)
  r <- ncol(data)
  k <- model.attributes$L
  n.iter <- dim(Gibbs$B)[1]
  if(is.na(burnin)){burnin <- round(n.iter/10)}
  times <- seq(burnin + 1, n.iter)
  GG <- model.attributes$GG
  npars <- r*k - sum(seq(1, k)) + r + 2*k^2 + k + k
  parmat <- cbind(rep(0, n.iter - burnin ))
  for(l in 1:k){parmat <- cbind(parmat, Gibbs$B[times, -seq(1,l), l])}
  for(l in seq(1, 2*k)){
    for(l2 in seq(l, 2*k)){
      parmat <- cbind(parmat, Gibbs$W[times, l, l2])
    }
  }
  parmat <- parmat[,-1] # Remove the first column
  parmat <- cbind(parmat, Gibbs$sigma2[times,],Gibbs$tau[times,])
  parcov <- cov(parmat)
  eigenparcov <- eigen(parcov)$values
  detparcov <- sum(log(eigenparcov))
  tot.result.mat <- matrix(NA, n.iter - burnin, n)
  intgrlike.mat <- matrix(NA, n.iter - burnin, n)
  
  for(g in (burnin + 1):n.iter){
    tot.result <- rep(npars/2*log(2*pi)/n + 0.5*detparcov/n,n)
    W0 <-Gibbs$W[g,,]
    V0 <-diag(Gibbs$sigma2[g,])
    B0 <- Gibbs$B[g,,]
    FF0 <- B0%*%model.attributes$FF
    a = matrix(0, nrow = n, ncol = 2*k)
    m = matrix(0, nrow = n, ncol = 2*k)
    A = array(0, dim = c(n, 2*k, r))
    f = matrix(0, nrow = n, ncol = r)
    e = matrix(0, nrow = n, ncol = r)
    
    # Define auxiliary arrays (T by p by arrays):
    R = array(0, dim = c(n, 2*k, 2*k))
    C = array(0, dim = c(n, 2*k, 2*k))
    Q = array(0,dim = c(n, r, r))
    
    # Kalman filter:
    a[1,] <- as.vector(model.attributes$GG %*% hyp.parm$m0)
    R[1,,] <- as.matrix(model.attributes$GG %*% hyp.parm$C0 %*% t(GG) + W0)
    f[1,] <- FF0 %*% a[1,]
    Q[1,,] <- as.matrix(FF0 %*% R[1,,] %*% t(FF0) + V0)
    A[1,,] <- R[1,,] %*% t(FF0) %*% solve(Q[1,,])
    e[1,] <- data[1,] - f[1,]
    m[1,] <- a[1,] + A[1,,] %*% e[1,]
    C[1,,] <- R[1,,] - A[1,,] %*% Q[1,,] %*% t(A[1,,])
    
    tot.result[1] <- tot.result[1] -(r/2)*log(2*pi) - 0.5*determinant(Q[1,,])$modulus[1] -0.5*t(matrix(e[1,]))%*%solve(Q[1,,])%*%matrix(e[1,])
    
    for(t in 2:n){
      a[t,] <- as.vector(GG %*% m[t-1,])
      R[t,,] <- GG %*% C[t-1,,] %*% t(GG) + W0
      f[t,] <- FF0 %*% a[t,]
      Q[t,,] <- FF0 %*% R[t,,] %*% t(FF0) + V0
      A[t,,] <- R[t,,] %*% t(FF0) %*% solve(Q[t,,])
      e[t,] <- data[t,] - f[t,]
      m[t,] <- a[t,] + A[t,,] %*% e[t,]
      C[t,,] <- R[t,,] - A[t,,] %*% Q[t,,] %*% t(A[t,,])
      
      tot.result[t] <- tot.result[t] -(r/2)*log(2*pi) - 0.5*determinant(Q[t,,])$modulus[1] -0.5*t(matrix(e[t,]))%*%solve(Q[t,,])%*%matrix(e[t,])
    }
    
    mar1 <- sum(dinvgamma(Gibbs$sigma2[g,], hyp.parm$n.sigma/2, hyp.parm$n.s2.sigma/2, log=TRUE))
    mar2 <- dinvwishart(Gibbs$W[g,,], hyp.parm$n.w + n - 1, hyp.parm$Psi, log=TRUE)
    mar3 <- sum(dinvgamma(Gibbs$tau[g,], hyp.parm$n.tau/2, hyp.parm$n.s2.tau/2, log=TRUE))
    mar4 <- 0

    for(l in 1:k){
      locone <- 1:l
      loctwo <- (l+1):r
      mufix <- numeric(l)
      mufix[l] <- 1
      b.b <- cbind(numeric(r-l))
      b.b <- b.b - hyp.parm$Hplus[loctwo,locone] %*% solve(hyp.parm$Hplus[locone,locone]) %*% mufix
      b.Hplus <- hyp.parm$Hplus[loctwo,loctwo] - hyp.parm$Hplus[loctwo,locone] %*% solve(hyp.parm$Hplus[locone,locone]) %*% hyp.parm$Hplus[locone, loctwo]
      b.Hplus.eigen <- eigen(b.Hplus)$values[-(r - l)]
      log.bmar <- -(r - 1 - l)/2*log(2*pi) - 0.5*sum(log(b.Hplus.eigen)) - 0.5*t(Gibbs$B[g, (l+1):r, l] - b.b) %*% matrix.Moore.Penrose(b.Hplus) %*% (Gibbs$B[g, (l+1):r, l] - b.b)
      mar4 <- mar4 + log.bmar
    }
    
    tot.result.mat[g-burnin,] <- tot.result + rep(mar1/n,n) + rep(mar2/n,n) + rep(mar3/n,n) + rep(mar4/n,n)
    intgrlike.mat[g-burnin,] <- tot.result
    if(verbose & g %% 100 ==0){
      cat(k,"factors done",g - burnin,"th step\n")}
  }
  Maximum <- max(apply(tot.result.mat, 1, sum))
  result.list <- list(tot.result.mat, intgrlike.mat, Maximum, detparcov)
  names(result.list) <- c("Total", "IL", "Maximum", "detparcov")
  return(result.list)
}

