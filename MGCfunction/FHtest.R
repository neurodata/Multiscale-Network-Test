#' This function prints out network testing results applied to binary adjacency matrix based \href{https://cran.r-project.org/web/packages/amen/amen.pdf}{AME model}.
#'
#' @param A [n x n] binary adjacency matrix;
#' @param X  corresponds to the nodal attributes:
#' \describe{
#'    \item{a [nxd] matrix}
#'    \item{a vector of length [n]}
#' }
#' @param k.range is a range of dimension of multiplicative network factors;
#' @param n.iter is the number of iteration beyond burn-in.
#' @importFrom amen ame
#' @return p-values of each model having different dimension of network factors;
#' @return p-values are based on Wilks' statistic.
#' @author Youjin Lee
#' @export
#'

FH_test = function(A, X, k.range, n.iter){

  wilks.pvalues = c()
  
  for(k in 1:length(k.range)){
    
    fit = ame(A, burn = 1, nscan = n.iter, odens = 1,
               model = "bin", print = FALSE, rvar = TRUE, cvar = TRUE,
               R = k.range[k] , plot = FALSE, symmetric = TRUE) 
    factors = cbind(fit$APM, fit$U)
    
    # regression using Multivariate ANOVA
    anov = lm(factors ~ X)
    if(ncol(factors) == 1){
    	wilks.pvalues[k] = anova(anov)[1,5]
    }else{    
    	wilks.pvalues[k] = summary(manova(anov), test = "Wilks")$stats[1,6]
    }
  }
  
  return(wilks.pvalues)
  
}


#' This function prints out the network factors applied to a binary adjacency matrix based \href{https://cran.r-project.org/web/packages/amen/amen.pdf}{AME model}.
#' As the network factors depend on the pre-specified dimension of multiplicative factors, this returns a set of factors with size of \code{k.range}.
#'
#' @param A [n x n] binary adjacency matrix;
#' @param k.range is a range of dimension of multiplicative network factors;
#' @param n.iter is the number of iteration beyond burn-in.
#' @importFrom amen ame
#' @return a list of network factors;
#' @author Youjin Lee
#' @export
#'

FH_factor = function(A, k.range){

  factors = list()
  
  for(k in 1:length(k.range)){
    
    fit = ame(A, burn = 1, nscan = n.iter, odens = 1,
               model = "bin", print = FALSE, rvar = TRUE, cvar = TRUE,
               R = k.range[k] , plot = FALSE, symmetric = TRUE) 
    factors[[k]] <- cbind(fit$APM, fit$U) 
    }
  return(factors) 
}