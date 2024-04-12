
#' @noRd

tau.simulation <- function(model.attributes, hyp.parm, parm){
  
  tau.current <- parm$tau
  
  for(l in 1:model.attributes$L){
    n.tau.aux <- hyp.parm$n.tau + model.attributes$N - l
    locone <- 1:l
    loctwo <- (l+1):model.attributes$R
    bfix <- numeric(l)
    bfix[l] <- 1
    b.tau <- hyp.parm$H[loctwo,locone] %*% solve(hyp.parm$H[locone,locone]) %*% bfix
    b.tau <- cbind(parm$B[-(1:l),l] - b.tau)
    StateH.tau <- hyp.parm$H[loctwo,loctwo] - hyp.parm$H[loctwo,locone]%*%solve(hyp.parm$H[locone,locone])%*%hyp.parm$H[locone,loctwo]
    n.s2.tau.aux <- hyp.parm$n.s2.tau + t(b.tau)%*%StateH.tau%*%b.tau
    tau.current[l] <- 1 / rgamma(1, (n.tau.aux)/2, n.s2.tau.aux/2)
  }
  
  return(tau.current)
  
}

