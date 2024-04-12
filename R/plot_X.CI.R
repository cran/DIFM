
#' @title  Credible interval plot of common factors
#' @description  The functions builds the plot of 95\% confidence intervals of the common realizations, X. The black solid lines are the posterior mean and the dased lines are the 95\% confidence intervals.
#'
#' @param Gibbs  Result of Gibbs sampler from DIFM function.
#' @param burnin  Number of burn-in. If not specified, it uses the first tenths as burn-in period.
#' @param main.bool  Add title of the plots.
#' @param layout.dim Dimension of panel layout for multiple common factors. If not specificed, common factor plots are layout in one column.
#'
#' @importFrom stats quantile
#' @importFrom graphics layout lines points abline
#' @importFrom grDevices rgb
#'
#' @return  Credible interval plots of common factors
#'
#' @export

plot_X.CI <- function(Gibbs, burnin = NA, main.bool = FALSE, layout.dim = NA){
  
  if(is.na(burnin)){burnin <- round(dim(Gibbs$X)[1] / 10)}
  n.iter <- dim(Gibbs$X)[1]
  turns <- seq(burnin + 1, n.iter)
  n <- dim(Gibbs$X)[2]
  k <- dim(Gibbs$X)[3]
  r <- dim(Gibbs$B)[2]
  countOrder <- c("First", "Second", "Third", "Fourth", "Fifth" ,"Sixth", "Seventh", "Eighth", "Nineth", "Tenth")
  repermute <- 1:r
  tseq <- 1:n
  if( is.na(layout.dim[1]) ){layout.out <- cbind(1:k)}else{
    layout.out <- matrix(layout.dim[1]*layout.dim[2], layout.dim[1], layout.dim[2], byrow = TRUE)
  }
  
  X.Gibbs.CI <- array(NA, dim = c(n,k,3))
  X.Gibbs.CI[,,1] <- apply(Gibbs$X[turns,,], c(2,3), quantile, p = 0.025)
  X.Gibbs.CI[,,2] <- apply(Gibbs$X[turns,,], c(2,3), mean)
  X.Gibbs.CI[,,3] <- apply(Gibbs$X[turns,,], c(2,3), quantile, p = 0.975)
  
  layout(layout.out)
  for(i in 1:k){
    yrange <- c(min(X.Gibbs.CI[,i,1]), max(X.Gibbs.CI[,i,3]))
    if(main.bool){maint = paste(countOrder[i], "common factors")}
    if(!main.bool){maint = ""}
    plot(tseq, X.Gibbs.CI[,i,2], pch = 19, ylim = yrange, ylab = "", xlab = "Time", type = "l", lwd = 2, main = maint)
    lines(tseq, X.Gibbs.CI[,i,1], lwd=2, lty=2)
    lines(tseq, X.Gibbs.CI[,i,3], lwd=2, lty=2)
  }
  
}
