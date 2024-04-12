
#' @noRd

X.simulation <- function(model.attributes, parm){
  
  x.current <- parm$X
  for(tt in 1:model.attributes$N){
    x.current[tt,] = model.attributes$FF %*% parm$theta[tt,]
  }
  
  return(x.current)
  
}
