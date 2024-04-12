
#' @title  Run Dynamic ICAR Factors Model (DIFM), with C++ codes
#' @description This function runs Dynamic ICAR factors Model (DIFM), simulated from C++ codes
#'
#' @param model.attributes  Model attributes from \code{difm.model.attributes}
#' @param hyp.parm  Hyperparameters from \code{difm.hyp.parm}
#' @param data  The dataset
#' @param every  Save \code{every} iterations to final result
#' @param verbose Print out the iteration process
#'
#' @return  The Gibbs sampler of DIFM
#'
#' @export
#' 
DIFMcpp <- function(model.attributes, hyp.parm, data, every = 1, verbose = TRUE){
  
  g.save <- 1
  n.save <- floor(model.attributes$n.iter / every)
  data <- as.matrix(data)
  
  Gibbs <- difm.gibbs.store(model.attributes, n.save)
  parm <- NULL
  parm$B <- Gibbs$B[1,,]
  if(model.attributes$L == 1){parm$B <- cbind(parm$B)}
  parm$sigma2 <- Gibbs$sigma2[1,]
  parm$tau <- Gibbs$tau[1,]
  parm$W <- Gibbs$W[1,,]
  parm$X <- Gibbs$X[1,,]
  if(model.attributes$L == 1){parm$X <- cbind(parm$X)}
  parm$theta <- Gibbs$theta[1,,]
  
  for(g in 1:model.attributes$n.iter){
    parm$theta <- theta_simulation(model.attributes, hyp.parm, data, parm)
    parm$X <- X_simulation(model.attributes, parm)
    parm$B <- B_simulation(model.attributes, hyp.parm, data, parm)
    parm$sigma2 <- sigma2_simulation(model.attributes, hyp.parm, data, parm)
    parm$tau <- tau_simulation(model.attributes, hyp.parm, parm)
    parm$W  <- W_simulation(model.attributes, hyp.parm, parm)
    
    if(g == 1 & every == 1){
      Gibbs$theta[g.save,,] <- parm$theta
      Gibbs$X[g.save,,] <- parm$X
      Gibbs$B[g.save,,] <- parm$B
      Gibbs$sigma2[g.save,] <- parm$sigma2
      Gibbs$tau[g.save,] <- parm$tau
      Gibbs$W[g.save,,] <- parm$W
      # g.save <- g.save + 1
    }
    
    if(g %% every == 0){
      Gibbs$theta[g.save,,] <- parm$theta
      Gibbs$X[g.save,,] <- parm$X
      Gibbs$B[g.save,,] <- parm$B
      Gibbs$sigma2[g.save,] <- parm$sigma2
      Gibbs$tau[g.save,] <- parm$tau
      Gibbs$W[g.save,,] <- parm$W
      g.save <- g.save + 1
    }
    
    if(verbose & (g %% 100 == 0)){
      cat("Factors =", model.attributes$L, ", current step:", g, "\n")
    }
    
  }
  
  return(Gibbs)
  
}
