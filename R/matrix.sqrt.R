#' @noRd

matrix.sqrt <- function(H){
  # Computes square root of nonnegative definite symmetric matrix using spectral decomposition
  
  if(nrow(H)==1) {H.sqrt = matrix(sqrt(H),nrow=1,ncol=1)} else
  {
    H.eigen = eigen(H)
    H.eigen.values = H.eigen$values    
    H.eigen.values[abs(H.eigen$values) < 10^(-10)] = 0
    H.sqrt = H.eigen$vectors %*% diag(sqrt(H.eigen.values)) %*% t(H.eigen$vectors)
  }  
  
  H.sqrt
}
# compile function to "make it go faster"
#matrix.sqrt = cmpfun(aux)
