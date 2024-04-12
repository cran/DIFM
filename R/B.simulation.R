
#' @noRd

B.simulation <- function(model.attributes, hyp.parm, data, parm){
  
  b.current <- parm$B
  
  Bvar <- list()
  Bvar.inv <- list()
  kseq <- c()
  mufix <- numeric(sum(1:model.attributes$L))
  mufix[cumsum(1:model.attributes$L)] <- 1
  for(l1 in 1:model.attributes$L){
    for(l2 in 1:l1){kseq <- c(kseq, model.attributes$R*(l1-1) + l2)}}
  for(l in 1:model.attributes$L){
    Bvar[[l]] <- hyp.parm$Hplus * parm$tau[l]
    Bvar.inv[[l]] <- hyp.parm$H / parm$tau[l]
  }
  Bvar.inv <- bdiag(Bvar.inv)
  V.inv <- diag(1 / parm$sigma2)
  Sig.b.inv <- Bvar.inv
  mu.b <- rep(0, model.attributes$R*model.attributes$L)
  for(i in 1:model.attributes$N){
    x.star <- kronecker(t(parm$X[i,]),diag(model.attributes$R))
    Sig.b.inv <- Sig.b.inv + t(x.star) %*% V.inv %*% x.star
    mu.b <- mu.b + t(x.star) %*% V.inv %*% cbind(as.numeric(data[i,]))
  }
  Sig.b <- solve(Sig.b.inv)
  mu.b <- Sig.b %*% mu.b
  mu.b.star <- mu.b[-kseq] + Sig.b[-kseq,kseq] %*% solve(Sig.b[kseq,kseq]) %*% cbind(mufix - mu.b[kseq])
  if(model.attributes$L > 1){
    Sig.b.star<- Sig.b[-kseq,-kseq] - Sig.b[-kseq,kseq] %*% matrix.Moore.Penrose(Sig.b[kseq,kseq]) %*% Sig.b[kseq,-kseq]
  }
  if(model.attributes$L == 1){
    Sig.b.star<- Sig.b[-kseq,-kseq] - cbind(Sig.b[-kseq,kseq]) %*% rbind(Sig.b[kseq,-kseq]) / Sig.b[kseq,kseq]
  }
  
  ### START FROM HERE!!
  
  aux <- mu.b.star + matrix.sqrt(Sig.b.star) %*% rnorm(model.attributes$R*model.attributes$L - sum(1:model.attributes$L))
  for(l in 1:model.attributes$L){
    b.current[(l+1):model.attributes$R,l] = aux[(model.attributes$R*(l-1) + 1 - sum(0:(l-1))):(model.attributes$R*l - sum(1:l))]
  }
  
  return(b.current)
  
}
