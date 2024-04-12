
#' @title Credible interval plot of factor loadings variance
#' @description  It returns a credible interval plot of factor loadings covariance, tau. The lines are 95% intervals, while the circles are posterior mean.
#'
#' @param Gibbs  Result of Gibbs sampler from DIFM function.
#' @param burnin  Number of burn-in. If not specified, it uses the first tenths as burn-in period.
#' @param true.val  True values of tau. If not available, NA.
#' @param main.bool  Add title of the plots.
#'
#' @return  Credible interval plot of tau
#' @export

plot_tau.CI <- function(Gibbs, burnin = NA, true.val = NA, main.bool = TRUE){
  n.iter <- dim(Gibbs$tau)[1]
  n.factors <- dim(Gibbs$tau)[2]
  if(is.na(burnin)){burnin <- round(n.iter / 10)}
  turns <- seq(burnin + 1, n.iter)
  
  tau.CI <- matrix(NA, 3, n.factors)
  tau.CI[1,] <- apply(Gibbs$tau[turns,], 2, quantile, 0.025)
  tau.CI[2,] <- apply(Gibbs$tau[turns,], 2, mean)
  tau.CI[3,] <- apply(Gibbs$tau[turns,], 2, quantile, 0.975)
  
  yrange <- c(min(tau.CI[1,]), max(tau.CI[3,]))
  if(main.bool){maint = expression(paste(tau, " Confidence Interval"))}
  if(!main.bool){maint = ""}
  plot(tau.CI[2,], main = maint, pch = 19, ylim = yrange, xlab = "Variables", ylab = "", xaxt = "n")
  axis(1, at = 1:n.factors, labels = 1:n.factors)
  for(i in 1:n.factors){lines(rep(i,2), c(tau.CI[1,i], tau.CI[3,i]))}
  if(!is.na(true.val[1])){points(1:n.factors, true.val, pch = 17, cex = 1.2, col = rgb(0,0,1,0.5))}
}
