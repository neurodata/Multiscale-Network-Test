varyA = function(popn, p, q, tau){
  n.group = 2
  ## 1. generate an outcome variable
  outcome = rbern(popn, 0.5)
  theta = runif(popn, 1-tau, 1+tau)
  
  ## 2. divide sample into two blocks
  multip = matrix(0, nrow = popn, ncol = n.group)
  for(i in 1:popn){
    if(outcome[i] <= 0.5){
      multip[i,] = c(0.6, 0.4)
    }else{
      multip[i,] = c(0.4, 0.6)
    }
  }

  group = rep(0, popn)
  for(i in 1:popn){
    group[i] = which(rmultinom(1, 1, multip[i,]) == 1)
  }
  
  ## 3. generate an edge variable
  adj = matrix(0, popn, popn)
  for(i in 1:popn){
  	for(j in i:popn){
  		if(i == j){
  			adj[i,j] = 0
  		}else{
  			adj[i,j] = rbinom(1, 1, theta[i]*theta[j]*q)
  		}
  		adj[j,i] = adj[i,j]
  	}
  }
  
  prob = matrix(0, popn, popn)
  group1.index = which(group == 1)
  group2.index = which(group == 2)
  
  # for group1
  for(i in 1:length(group1.index)){
    for(j in i:length(group1.index)){
      id1 = group1.index[i]
      id2 = group1.index[j]
      if (id1 == id2){
        prob[id1, id2] = 0.0
      }else{
        prob[id1, id2] = p*theta[id1]*theta[id2]
      }  
      adj[id1, id2] = rbinom(1,1, prob[id1, id2])
      adj[id2, id1] = adj[id1, id2]
    }
  }

  # for group2
  for(i in 1:length(group2.index)){
    for(j in i:length(group2.index)){
      id1 = group2.index[i]
      id2 = group2.index[j]
      if (id1 == id2){
        prob[id1, id2] = 0.0
      }else{
        prob[id1, id2] = p*theta[id1]*theta[id2]
      }  
      adj[id1, id2] = rbinom(1,1,prob[id1, id2])
      adj[id2, id1] = adj[id1, id2]
    }
  }

  G = graph.adjacency(adj, "undirected")
  G = connect_graph(G)
  V(G)$outcome = outcome
  V(G)$group = group

  return(G)
  
}
