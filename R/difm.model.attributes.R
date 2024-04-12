
#' @title  Initialize model attributes for DIFM
#' @description  It initialize the basic parameters and model attributes for DIFM
#'
#' @param data  The dataset
#' @param n.factors  Number of factors to run DIFM
#' @param n.iter Number of iterations
#' @param G0  The basic evolution matrix for one factor
#'
#' @return  A list of number of timepoints, subregions, factors, matrix of evolution matrix, and matrix to extract common factors.
#' @export

difm.model.attributes <- function(data, n.iter, n.factors, G0){
  
  N <- nrow(data)
  R <- ncol(data)
  L <- n.factors
  Gpower <- ncol(G0)
  GG <- matrix(0, n.factors * Gpower, n.factors * Gpower)
  FF <- matrix(0, n.factors, n.factors * Gpower)
  for(i in 1:n.factors){
    GG[(1 + Gpower*(i-1)):(Gpower + Gpower*(i-1)), (1 + Gpower*(i-1)):(Gpower + Gpower*(i-1))] <- G0
    FF[i,(1 + Gpower*(i-1))] <- 1
  }
  
  model.attributes <- list(N, R, L, n.iter, GG, FF)
  names(model.attributes) <- c("N", "R", "L", "n.iter", "GG", "FF")
  
  return(model.attributes)
  
}
