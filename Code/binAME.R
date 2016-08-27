library(igraph)

binAME <- function(popn, link.g){
	# this function returns a graph objects 
	# which is generated from additive and multiplicative model
	# of uniformly distributed latent factor
	###########################################################
	# # input
	# link.g : measurable function of g : [0,1]^2 -> [0,1]
	# popn : number of 
	############################################################ 
	
	U <- runif(popn, 0, 1)
	X <- rbinom(popn, 1, U)
	
  	A <- matrix(0, popn, popn) 
 
  	for(i in 1: (popn-1)){
    	for(j in (i+1):popn){
      		p <- link.g(U[i],  U[j])
      		A[i,j] <- rbinom(1,1,p)
      		A[j,i] <- A[i,j]
    	}
  	}
  		
  	G  <- graph.adjacency(A)
	V(G)$X <- X
	
	return(G)
	
}
