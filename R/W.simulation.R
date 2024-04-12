
#' @noRd

W.simulation <- function(model.attributes, hyp.parm, parm){
  
  W.current <- parm$W
  
  n.w.aux <- hyp.parm$n.w + model.attributes$N - 1
  n.s2.w.aux <- matrix(0, nrow(model.attributes$GG), nrow(model.attributes$GG))
  for(tt in 2:model.attributes$N){ 
    n.s2.w.aux <- n.s2.w.aux + crossprod(t(parm$theta[tt,]- model.attributes$GG%*%parm$theta[tt-1,]))
  }
  W.current <- rinvwishart(n.w.aux, hyp.parm$Psi + n.s2.w.aux)
  
  return(W.current)
  
}
