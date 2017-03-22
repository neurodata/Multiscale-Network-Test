library(HHG)
library(DTMCPack)
library(Matrix)
library(energy)
library(ecodist)
library(ADMTools)

source(diffmaps.R)

NetworkTest.q = function(G, X, option, diffusion, dstep, n.perm, q){
  # This function shows prints out the pvalues of testing
  # independence between network G and nodal attributes X 
  # through testing diffusion maps U and X at given t=dstep and q = q.
  ################################################
  # # Input 
  # G : igraph object - right now restrict it to connected, undirected graph
  # X : matrix of any dimensional nodal attributes
  # option : 1,2,3 for MGC, dCov, HHG applied to Euclidean distance
  # diffusion : if TRUE, distance metric for G is diffusion distance across dstep 
  # dstep :  diffusion time you are going to use in testing
  # n.perm : the number of permutations
  # q : dimension of diffusion maps
  #################################################
  # require : HHG, DTMCPack, Matrix, energy, ecodist, ADMTools
  #################################################
  # # Output
  # pvalues : a vector of pvalues from selected testing method
  #           In case of MGC, pvalue maps are printed aout
  #           In case of hhg, pvalue is based on likelihood ratio scores statistic
  #################################################
  
  if(!is.connected(G)){
    G = connect_graph(G)
  }
  
  A = as.matrix(get.adjacency(G))
  P = A / pmax(rowSums(A), 1)
  
  if(diffusion == FALSE){ 
  	Dx = as.matrix(dist((A)), diag = TRUE, upper = TRUE) 
  	Dy = as.matrix(dist((X)), diag = TRUE, upper = TRUE)
  	
  	if(option == 1){
  		result = MGCPermutationTest(Dx, Dy, n.perm, option = 'mcor')
  		return(result)
  	}else if(option == 2){
  		pvalues = dcov.test(Dx, Dy, index = 1.0, R = n.perm)$p.value
  		return(pvalues)
  	}else if(option == 3){
  		pvalues = hhg.test(Dx, Dy, nr.perm = n.perm)$perm.pval.hhg.sl
  		return(pvalues)
  	}
  }
  
  U.list = dmap.q(P, dstep, q) 
  
  all.result = list()
  pvalues = c()
  
  for(s in 1:length(dstep)){
    
    U = U.list[[s]]
    Dx = as.matrix(dist(U), diag = TRUE, upper = TRUE) 
    Dy = as.matrix(dist(X), diag = TRUE, upper = TRUE)
    
    if(option == 1){
      mgc = MGCPermutationTest(Dx, Dy, n.perm, option = 'mcor')
      all.result[[s]] = mgc
    }
    
    if(option == 2){
      dco = dcor.ttest(Dx, Dy, distance = TRUE)
      pvalues[s] = dco$p.value
    }
    
    
    if(option == 3){
      hhg = hhg.test(Dx, Dy, nr.perm = n.perm)
      pvalues[s] = hhg$perm.pval.hhg.sl
    }
  }
  
  if(option == 1){
    return(all.result)
  }else{
    return(pvalues)
  }
}