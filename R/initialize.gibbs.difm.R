
#' @title  Initialize the first Gibbs sample 
#' @description  It computes the spatial covariance and precision matrix of the neighboring subregions using Intrinsice Autoregressive Conditional (ICAR) process.
#' @details The off-digonal values are -1 when two subregions are neighbors. Otherwise, we assign 0. The diagonal values are the sum of the values of its own row.
#'
#' @param model.attributes  Model attributes from \code{difm.model.attributes}
#'
#' @return  A list of the initialized parameters.
#' @noRd

initialize.gibbs.difm <- function(model.attributes){
  
  n <- model.attributes$N
  r <- model.attributes$R
  k <- model.attributes$L
  GG <- model.attributes$GG
  
  tau.current <- rep(.01, k) # Strength of spatial relationship
  
  b.current <- matrix(0, nrow = r, ncol = k)
  for(l in 1:k) b.current[l,l] <- 1
  
  x.current <- matrix(0,nrow = n,ncol = k)
  x.current[,1] <- rnorm(n)
  
  theta.length <- nrow(GG)
  theta.current <- matrix(0, nrow = n, ncol = theta.length)
  theta.current[,1] <- x.current[,1]
  
  sigma2.current <- rep(1, r)
  
  W.current <- diag(2*k)
  
  Gibbs.current <- list(b.current, sigma2.current, x.current, theta.current, tau.current, W.current)
  names(Gibbs.current) <- c("B", "sigma2", "X", "theta", "tau", "W")
  
  return(Gibbs.current)
  
}
