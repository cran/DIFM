
#' @title  Hyperparameters for DIFM
#' @description  Sets the hyperparameters to generate Gibbs sampler of DIFM
#'
#' @param model.attributes  Model attributes from \code{difm.model.attributes}
#' @param n.tau  Shape parameter for tau
#' @param n.s2.tau  Rate parameter for tau
#' @param n.sigma  Shape parameter for sigma squared
#' @param n.s2.sigma  Rate parameter for sigma squared
#' @param Hlist Neighborhood matrix
#' @param Psi.size  The magnitude of covariance for the evolution matrix
#'
#' @return  A list of hyperparameters of tau, W, sigma, and theta.
#' @export

difm.hyp.parm <- function(model.attributes, n.tau = 2.2, n.s2.tau = .1, n.sigma = 2.2, n.s2.sigma = .1, Hlist, Psi.size = .01){
  
  k <- model.attributes$L
  polyorder <- nrow(model.attributes$GG) / k
  
  n.tau <- n.tau
  n.s2.tau <- n.s2.tau
  n.w <- polyorder*k + 2
  n.sigma <-  n.sigma
  n.s2.sigma <- n.s2.sigma
  m0 <- rep(0, polyorder*k)
  C0 <- 10^4 * diag(polyorder*k)
  Psi <- diag(polyorder*k) * Psi.size
  
  hyp.parm <- list(n.tau, n.s2.tau, n.w, n.sigma, n.s2.sigma, m0, C0, Psi, Hlist$H, Hlist$Hplus)
  names(hyp.parm) <- c("n.tau", "n.s2.tau", "n.w", "n.sigma", "n.s2.sigma", "m0", "C0", "Psi", "H", "Hplus")
  return(hyp.parm)
  
}
