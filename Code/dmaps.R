library(DTMCPack)
library(Matrix)

#### print out diffusion maps from undirected grph
dmap <- function(trans, times){
  # require : DTMCPack, Matrix
  # # input
  # trans is a transition probability
  # times is a set of diffusion times you want to print out
  ################################################
  maps <- list()
  
  # stationary distribution (probability)
  pi <- statdistr(trans)
  pi.mat1 <- as.matrix(Diagonal(length(pi), pi))
  pi.mat2 <- as.matrix(Diagonal(length(pi), pi^(-1/2)))
  
  # symmetric kernel
  Q <- sqrt(pi.mat1) %*% trans %*% pi.mat2
  
  lambda <- eigen(Q)$values
  psi <- eigen(Q)$vectors
  
  phi <- pi.mat2 %*% psi
  
  Lambda <- as.matrix(Diagonal(length(lambda), lambda))

  # each row of maps is a diffusion maps of each vertex
  for(t in 1:length(times)){
    maps[[t]] <- phi %*% Lambda^(times[t])  
  }

  return(maps)
}