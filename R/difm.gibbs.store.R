
#' @title  List of Gibbs sampler of each parameter
#' @description  A Gibbs sampler list of parameters that are generated each iteration. The first sample is generated through \code{initialize.gibbs.difm}.
#'
#' @param model.attributes  Model attributes from \code{difm.model.attributes}
#' @param n.save  Number of total iterations to save the results
#'
#' @return  A list of Gibbs sample of factor loadings B, sigma squared, spatial strenth parameter tau, common factors X, temporal components theta, and evolution covariance W.
#' @noRd

difm.gibbs.store <- function(model.attributes, n.save){
  
  b.Gibbs <- array(NA, dim = c(n.save, model.attributes$R, model.attributes$L))
  sigma2.Gibbs <- matrix(NA, nrow = n.save, ncol = model.attributes$R)
  tau.Gibbs <- matrix(NA, nrow = n.save, ncol = model.attributes$L)
  x.Gibbs <- array(NA, dim = c(n.save, model.attributes$N, model.attributes$L))
  theta.Gibbs <- array(NA, dim = c(n.save, model.attributes$N, ncol(model.attributes$GG)))
  W.Gibbs <- array(NA, dim = c(n.save, ncol(model.attributes$GG), ncol(model.attributes$GG)))
  
  current.Gibbs <- initialize.gibbs.difm(model.attributes)
  
  b.Gibbs[1,,] <- current.Gibbs$B
  sigma2.Gibbs[1,] <- current.Gibbs$sigma2
  x.Gibbs[1,,] <- current.Gibbs$X
  theta.Gibbs[1,,] <- current.Gibbs$theta
  tau.Gibbs[1,] <- current.Gibbs$tau
  W.Gibbs[1,,] <- current.Gibbs$W
  
  Gibbs <- list(b.Gibbs, sigma2.Gibbs, x.Gibbs, theta.Gibbs, tau.Gibbs, W.Gibbs)
  names(Gibbs) <- c("B", "sigma2", "X", "theta", "tau", "W")
  return(Gibbs)
  
}
