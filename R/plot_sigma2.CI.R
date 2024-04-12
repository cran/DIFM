
#' @title A credible interval plot of posterior of sigma squared
#' @description  It returns a credible interval plot of idiosyncratic variance, sigma squared. The lines are 95% intervals, while the circles are posterior mean.
#'
#' @param Gibbs  Result of Gibbs sampler from DIFM function.
#' @param burnin  Number of burn-in. If not specified, it uses the first tenths as burn-in period.
#' @param permutation  Permutation of variables. If not specified, no permutation.
#' @param main.bool  Add title of the plots.
#'
#' @return  A credible interval plot of sigma squared
#' @export

plot_sigma2.CI <- function(Gibbs, burnin = NA, permutation = NA, main.bool = TRUE){
  n.iter <- dim(Gibbs$sigma2)[1]
  r <- dim(Gibbs$sigma2)[2]
  if(is.na(burnin)){burnin <- round(n.iter / 10)}
  turns <- seq(burnin + 1, n.iter)
  repermute <- 1:r
  if( !is.na(permutation[1]) ){repermute <- order(permutation)}
  
  sigma2.CI <- matrix(NA, 3, r)
  sigma2.CI[1,] <- apply(Gibbs$sigma2[turns,repermute], 2, quantile, 0.025)
  sigma2.CI[2,] <- apply(Gibbs$sigma2[turns,repermute], 2, mean)
  sigma2.CI[3,] <- apply(Gibbs$sigma2[turns,repermute], 2, quantile, 0.975)
  
  yrange <- c(min(sigma2.CI[1,]), max(sigma2.CI[3,]))
  if(main.bool){maint = expression(paste(sigma^2, " Confidence Interval"))}
  if(!main.bool){maint = ""}
  plot(sigma2.CI[2,], main = maint, pch = 19, ylim = yrange, xlab = "Variables", ylab = "")
  for(i in 1:r){lines(rep(i,2), c(sigma2.CI[1,i], sigma2.CI[3,i]))}
}
