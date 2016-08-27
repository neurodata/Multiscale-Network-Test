###### network test
source("dmap.R")
source("disToRanks.R")
source("MGCScaleVerify.R")
source("MGCPermutationTest.R")
source("LocalCorr.R")
source("connect_graph.R")

NetworkTest <- function(G, X, option, dstep, n.perm){
  # Auther : Youjin Lee
  # This function shows prints out the pvalues of testing
  # independence between network G and nodal attributes X 
  # through testing diffusion maps U and X
  ################################################
  # # Input 
  # G : igraph object - right now restrict it to connected, undirected graph
  # X : matrix of any dimensional nodal attributes
  # option : 1,2,3 for MGC, dCov, HHG
  # dstep :  diffusion time you are going to use in testing
  # n.perm : the number of permutations
  #################################################
  # require : HHG, DTMCPack, Matrix, igraph, energy
  #################################################
  # # Output
  # pvalues : a vector of pvalues from selected testing method
  #           In case of MGC, pvalue maps are printed aout
  #           In case of hhg, pvalue is based on likelihood ratio scores statistic
  #################################################
  
  if(!is.connected(G)){
    G <- connect_graph(G)
  }
  
  A <- as.matrix(get.adjacency(G))
  P <- A / rowSums(A)
  
  U.list <- dmap(P, dstep) 
  
  pvalue.map <- list()
  pvalues <- c()
  
  for(s in 1:length(dstep)){
    
    U <- U.list[[s]]
    Dx = as.matrix(dist((U)), diag = TRUE, upper = TRUE) 
    Dy = as.matrix(dist((X)), diag = TRUE, upper = TRUE)
    
    if(option == 1){
      mgc <- MGCPermutationTest(Dx, Dy, n.perm, option = 1)
      pvalue.map[[s]] <- mgc
    }
    
    if(option == 2){
      dco <- dcov.test(Dx, Dy, index = 1.0, R = n.perm)
      pvalues[s] <- dco$p.value
    }
    
    
    if(option == 3){
      hhg <- hhg.test(Dx, Dy, nr.perm = n.perm)
      pvalues[s] <- hhg$perm.pval.hhg.sl
    }
  }
  
  if(option==1){
    return(pvalue.map)
  }else{
    return(pvalues)
  }
}
