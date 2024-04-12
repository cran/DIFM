#' @title  Matrix Moore Penrose Matrix
#' @description  It calcualtes the inverse of Moore-Penrose form
#'
#' @param H  The matrix that needs Moore-Penrose inveres
#'
#' @return  Moore-Penrose inveres of \code{H}
#' @noRd

matrix.Moore.Penrose <- function(H){
  # Computes Moore Penrose generalized inverse of symmetric matrix using spectral decomposition
  
  H.eigen = eigen(H)
  inverse.values = rep(0,nrow(H))
  inverse.values[abs(H.eigen$values) > 10^(-7)] = 1 / H.eigen$values[abs(H.eigen$values) > 10^(-7)]
  H.MP = H.eigen$vectors %*% diag(inverse.values) %*% t(H.eigen$vectors) 
  
  return(H.MP)
}
# compile function to "make it go faster"
#matrix.Moore.Penrose = cmpfun(aux1)
