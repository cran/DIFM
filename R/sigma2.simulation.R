
#' @noRd

sigma2.simulation <- function(model.attributes, hyp.parm, data, parm){
  
  sigma2.current <- parm$sigma2
  
  for(j in 1:model.attributes$R){
    n.sigma2.aux <- hyp.parm$n.sigma + model.attributes$N
    n.s2.sigma2.aux <- hyp.parm$n.s2.sigma 
    for (i in 1:model.attributes$N) n.s2.sigma2.aux <- n.s2.sigma2.aux + (data[i,j] - parm$B[j,] %*% parm$X[i,])^2
    sigma2.current[j] = 1 / rgamma(1, n.sigma2.aux/2, n.s2.sigma2.aux/2)
  }
  
  return(sigma2.current)
  
}
