#' @title  Order of permutation by the largest absolute value in each eigenvector
#' @description It finds the vector of permutation to permute data by its largest absolute value in each eigenvector. It sets the order by specified number of factors, and the rest is ordered as they were.
#'
#' @param data  The dataset
#' @param n.factors  Number of factors
#'
#' @return  The numeric vector of permutation
#' @export

permutation.order <- function(data, n.factors){

  r <- dim(data)[2]
  eigvec <- eigen(cor(data))$vectors
  rownames(eigvec) <- 1:r
  colnames(eigvec) <- 1:r
  colseq <- seq(1, r)
  win.factors <- numeric(n.factors)
  win.factors[1] <- colseq[rank(abs(eigvec[,1])) == r]
  for(i in 2:n.factors){
    choose_again <- TRUE
    j <- 0
    while(choose_again){
      win.factors[i] <- colseq[rank(abs(eigvec[,i])) == (r-j)]
      choose_again <- any(win.factors[i] == win.factors[1:(i-1)])
      j <- j + 1
    }
  }
  rest.vars <- colseq[-win.factors]
  permutation <- c(win.factors, rest.vars)
  return(permutation)
}
