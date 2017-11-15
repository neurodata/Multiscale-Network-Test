library(amen)

FH_factor = function(A, k.range){
  # this function test independence between network topology and nodal attributes
  # via Fosdick & Hoff (2015) 's method of testing where A is binary
  ################################################
  # # input
  # A : binary adjacency matrix of network 
  # X : any dimensional matrix of nodal attributes
  # k.range : range of dimension of multiplicative factor model
  # n.iter : the number of iteration beyond burn-in
  ################################################
  # require : amen
  #################################################
  # # output
  # factors : estimated network factors
  ################################################
  factors = list()
  
  for(k in 1:length(k.range)){
    
    fit = ame(A, burn = 1, nscan = n.iter, odens = 1,
               model = "bin", print = FALSE, rvar = TRUE, cvar = TRUE,
               R = k.range[k] , plot = FALSE, symmetric = TRUE) 
    factors[[k]] = cbind(fit$APM, fit$U)
    
    }
  
  return(factors)
  
}

FH_test = function(A, X, k.range, n.iter){
  # this function test independence between network topology and nodal attributes
  # via Fosdick & Hoff (2015) 's method of testing where A is binary
  ################################################
  # # input
  # A : binary adjacency matrix of network 
  # X : any dimensional matrix of nodal attributes
  # k.range : range of dimension of multiplicative factor model
  # n.iter : the number of iteration beyond burn-in
  ################################################
  # require : amen
  #################################################
  # # output
  # pvalues : pvalues of each model of different dimension of factors
  #           pavlues are based on Wilks' statistic 
  ################################################
  wilks.pvalues = c()
  
  for(k in 1:length(k.range)){
    
    fit = ame(A, burn = 1, nscan = n.iter, odens = 1,
               model = "bin", print = FALSE, rvar = TRUE, cvar = TRUE,
               R = k.range[k] , plot = FALSE, symmetric = TRUE) 
    factors = cbind(fit$APM, fit$U)
    
    # regress using Multivariate ANOVE
    anov = lm(factors ~ X)
    if(ncol(factors) == 1){
    	wilks.pvalues[k] = anova(anov)[1,5]
    }else{    
    	wilks.pvalues[k] = summary(manova(anov), test = "Wilks")$stats[1,6]
    }
  }
  
  return(wilks.pvalues)
}