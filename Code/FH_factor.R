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
  # pvalues : pvalues of each model of different dimension of factors
  #           pavlues are based on Wilks' statistic 
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