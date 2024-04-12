
#' @title  Spatial dependence matrix of the factor loadings
#' @description  It computes the spatial covariance and precision matrix of the neighboring subregions using Intrinsice Autoregressive Conditional (ICAR) process.
#' @details The off-digonal values are -1 when two subregions are neighbors. Otherwise, we assign 0. The diagonal values are the sum of the values of its own row.
#'
#' @param areapoly  The polygon of the areas. We can obtain this through \code{readOGR} function from \code{sp} matrix. 
#' @param permutation  Permutation order of the subregions
#'
#' @importFrom spdep poly2nb nb2mat 
#'
#' @return  A list of two matrices: Precision matrix H and the covariance matrix obtained through Moore-Penrose inverse of H.
#' @export
buildH <- function(areapoly, permutation = NA){
  
  if(is.na(permutation[1])){
    permutation <- seq(1, nrow(areapoly))
  }
  areaNB <- poly2nb(areapoly)
  areaNBmat <- nb2mat(areaNB,style="B")
  areaH <- areaNBmat[permutation,permutation]
  areaH[areaH == 1]<--1
  R <- ncol(areaH)
  for(i in 1:R){areaH[i,i] <- -sum(areaH[i,-i])}
  Hplus<-matrix.Moore.Penrose(areaH)
  
  Hlist <- list(areaH, Hplus)
  names(Hlist) <- c("H", "Hplus")
  return(Hlist)
  
}
