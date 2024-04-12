
#' @title  Credible interval plot of factor loadings
#' @description  The functions builds a column-wise plots of factor loadings. The parameters fixed at 1 are displayed with red dashed vertical lines.
#'
#' @param Gibbs  Result of Gibbs sampler from DIFM function 
#' @param true.val  True values of factor loadings. If not available, NA.  
#' @param burnin  Number of burn-in. If not specified, it uses the first tenths as burn-in period.
#' @param permutation  Permutation of variables. If not specified, no permutation.
#' @param main.bool  Add title of the plots.
#' @param layout.dim  Dimension of panel layout for multiple factor loadings. If not specificed, factor loadings plots are layout in one column.
#'
#' @importFrom stats quantile
#' @importFrom graphics layout lines points abline
#' @importFrom grDevices rgb
#'
#' @return  Factor loadings credible interval plots
#'
#' @export

plot_B.CI <- function(Gibbs, true.val = NA, burnin = NA, permutation = NA, main.bool = TRUE, layout.dim = NA){
  if(is.na(burnin)){burnin <- round(dim(Gibbs$B)[1] / 10)}
  n.iter <- dim(Gibbs$B)[1]
  turns <- seq(burnin + 1, n.iter)
  r <- dim(Gibbs$B)[2]
  k <- dim(Gibbs$B)[3]
  countOrder <- c("First", "Second", "Third", "Fourth", "Fifth" ,"Sixth", "Seventh", "Eighth", "Nineth", "Tenth")
  repermute <- 1:r
  if( !is.na(permutation[1]) ){repermute <- order(permutation)}
  if( is.na(layout.dim[1]) ){layout.out <- cbind(1:k)}else{
    layout.out <- matrix(layout.dim[1]*layout.dim[2], layout.dim[1], layout.dim[2], byrow = TRUE)
  }
  
  B.Gibbs.CI <- array(NA, dim = c(r,k,3))
  B.Gibbs.CI[,,1] <- apply(Gibbs$B[turns,repermute,], c(2,3), quantile, p = 0.025)
  B.Gibbs.CI[,,2] <- apply(Gibbs$B[turns,repermute,], c(2,3), mean)
  B.Gibbs.CI[,,3] <- apply(Gibbs$B[turns,repermute,], c(2,3), quantile, p = 0.975)
  
  layout(layout.out)
  for(i in 1:k){
    yrange <- c(min(B.Gibbs.CI[,i,1]), max(B.Gibbs.CI[,i,3]))
    if(main.bool){maint = paste(countOrder[i], "factor loadings")}
    if(!main.bool){maint = ""}
    plot(B.Gibbs.CI[,i,2], pch = 19, ylim = yrange, ylab = "", xlab = "Variables", main = maint)
    for(j in 1:r){lines(c(j,j), c(B.Gibbs.CI[j,i,1], B.Gibbs.CI[j,i,3]))}
    if(!is.na(true.val[1])){
      points(1:r, true.val[,i], pch = 17, cex = 1.2, col = rgb(0,0,1,0.5))
    }
    abline(h = 0, col = "gray", lty = 2)
    abline(v = permutation[i], col = "red", lty = 2)
  }
}
