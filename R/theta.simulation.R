
#' @noRd

theta.simulation <- function(model.attributes, hyp.parm, data, parm){
  # Simulate the latent process
  
  p <- ncol(model.attributes$GG)
  FF0 <- parm$B %*% model.attributes$FF
  
  # Define auxiliary matrices (T by p matrices):
  a <- matrix(0, nrow = model.attributes$N, ncol = p)
  m <- matrix(0, nrow = model.attributes$N, ncol = p)
  A <- array(0, dim = c(model.attributes$N, p, model.attributes$R))
  f <- matrix(0, nrow = model.attributes$N, ncol = model.attributes$R)
  e <- matrix(0, nrow = model.attributes$N, ncol = model.attributes$R)
  
  # Define auxiliary arrays (T by p by arrays):
  R <- array(0, dim = c(model.attributes$N, p, p))
  C <- array(0, dim = c(model.attributes$N,p,p))
  Q <- array(0, dim = c(model.attributes$N, model.attributes$R, model.attributes$R))
  
  # This matrix will store the simulated theta_t, t=1,...,T:
  theta <- matrix(0, nrow = model.attributes$N, ncol = p)
  
  # Kalman filter:
  a[1,] <- as.vector(model.attributes$GG %*% hyp.parm$m0)
  R[1,,] <- as.matrix(model.attributes$GG %*% hyp.parm$C0 %*% t(model.attributes$GG) + parm$W)
  f[1,] <- FF0 %*% a[1,]
  Q[1,,] <- as.matrix(FF0 %*% R[1,,] %*% t(FF0) + diag(parm$sigma2))
  A[1,,] <- R[1,,] %*% t(FF0) %*% solve(Q[1,,])
  e[1,] <- data[1,] - f[1,]
  # m[1,] <- a[1,] + A[1,,] %*% e[1,]
  # C[1,,] <- R[1,,] - A[1,,] %*% Q[1,,] %*% t(A[1,,])
  
  C[1,,] <- solve(solve(R[1,,]) + t(FF0) %*% diag(1/parm$sigma2) %*% FF0 )
  m[1,] <- C[1,,] %*% (solve(R[1,,]) %*% a[1,] + t(FF0) %*% diag(1/parm$sigma2) %*% data[1,])
  
  for(tt in 2:model.attributes$N)
  {
    a[tt,] <- as.vector(model.attributes$GG %*% m[tt-1,])
    R[tt,,] <- model.attributes$GG %*% C[tt-1,,] %*% t(model.attributes$GG) + parm$W
    f[tt,] <- FF0 %*% a[tt,]
    Q[tt,,] <- FF0 %*% R[tt,,] %*% t(FF0) + diag(parm$sigma2)
    A[tt,,] <- R[tt,,] %*% t(FF0) %*% solve(Q[tt,,])
    e[tt,] <- data[tt,] - f[tt,]
    # m[t,] <- a[t,] + A[t,,] %*% e[t,]
    # C[t,,] <- R[t,,] - A[t,,] %*% Q[t,,] %*% t(A[t,,])
    
    C[tt,,] <- solve(solve(R[tt,,]) + t(FF0) %*% diag(1/parm$sigma2) %*% FF0 )
    m[tt,] <- C[tt,,] %*% (solve(R[tt,,]) %*% a[tt,] + t(FF0) %*% diag(1/parm$sigma2) %*% data[tt,])
  }
  
  # Now, the backward sampler:
  C.sqrt = matrix.sqrt(C[model.attributes$N,,])
  theta[model.attributes$N,] <- m[model.attributes$N,] + C.sqrt %*% rnorm(p, mean = 0, sd = 1)
  
  for(tt in (model.attributes$N - 1):1)
  {
    B <- C[tt,,] %*% t(model.attributes$GG) %*% matrix.Moore.Penrose(R[tt+1,,])
    h <- m[tt,] + B %*% (theta[tt+1,] - a[tt+1,])
    H <- C[tt,,] - B %*% R[tt+1,,] %*% t(B)   
    
    H.sqrt = matrix.sqrt(H)
    
    #    H <- solve(solve(C[t,,]) + t(GG) %*% solve(W) %*% GG ) 
    #    h <- H %*% (solve(C[t,,]) %*% m[t,] + t(GG) %*% solve(W) %*% theta[t+1,])
    theta[tt,] <- h + H.sqrt %*% rnorm(p, mean = 0, sd = 1)
  }
  
  return(theta)
} 
