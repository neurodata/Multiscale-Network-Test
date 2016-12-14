#############################
## FH test with multivariate-covariates ####
FH_test2 = function(A, X, k.range, n.iter){
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
  pvalues <- c()
  
  for(k in 1:length(k.range)){
    
    fit <- ame(A, burn = 1, nscan = n.iter, odens = 1,
               model = "bin", print = FALSE, rvar = TRUE, cvar = TRUE,
               R = k.range[k] , plot = FALSE, symmetric = TRUE) 
    factors <- cbind(fit$APM, fit$U)
    
    # regress using Multivariate ANOVE
    anov <-lm(factors ~ as.factor(X))
    
    
    if(ncol(factors) == 1){
      pvalues[k] <- Anova(anov)[1,4]
    }else{    
      outtests <- car:::print.Anova.mlm
      
      # allow the function to return the results and disable print
      body(outtests)[[16]] <- quote(invisible(tests))
      body(outtests)[[15]] <- NULL
      
      # Now run the regression
      mod <- lm(factors ~ as.factor(X))
      
      # Run the Anova over all tests  
      tab <- lapply(c("Pillai", "Wilks", "Roy"), 
                    function(i)  outtests(Anova(mod, test.statistic=i)))
      
      tab <- do.call(rbind, tab)
      row.names(tab) <- c("Pillai", "Wilks", "Roy")
      pvalues[k] =  tab[2,6]  
    }
    
    
  }
  
  return(pvalues)}

# collaboration; directed network
collab <- read.csv(file = "../Data/politics/climate0205-collab.csv", 
                   header = TRUE, row.names = 1, sep = ";") ## 34 x 34 matrix (not symmetric)

# type of organization; vector with five character types
types <- read.csv(file="../Data/politics/climate0205-type.csv", 
                  header = TRUE, row.names = 1, sep = ";")[,2] ## categorical 5 types


######## transform directed network into undirected network #########
undirected.G = graph.adjacency(as.matrix(collab), mode = "undirected")
Dy.types = matrix(0, 34, 34)
for(i in 1:34){
  for(j in 1:34){
    Dy.types[i,j] = ifelse(types[i] == types[j], 0, 1)
  }
}


#dstep = c(1:10) # diffusion times
n.perm = 500 # the number of permutations
#n.iter = n.perm # the number of iterations (FH test)
#k.range = c(0:10) # the dimension of multivariate factors in FH model

G.adj = as.matrix(get.adjacency(undirected.G))
#### n.random = number of random edges.

n.random = as.integer(34*c(1:10)*0.1)

# the number of edges 1651 / 

Dx.G1 = matrix(0, nrow = 20, ncol = 20); Dx.G2 =  matrix(0, nrow = 20, ncol = 20)
hhg.1 =  matrix(0, nrow = 20, ncol = 20); hhg.2 =  matrix(0, nrow = 20, ncol = 20)
mcorr.1 = matrix(0, nrow = 20, ncol = 20); mcorr.2 = matrix(0, nrow = 20, ncol = 20)
fh.G = matrix(0, nrow = 20, ncol = 20)
mgc.1 = matrix(0, nrow = 20, ncol = 20);  mgc.2 = matrix(0, nrow = 20, ncol = 20)

hhg.5 = matrix(0, nrow=20, ncol = 20); mcorr.5 = matrix(0, nrow = 20, ncol = 20)
mgc.5 = matrix(0, nrow = 20, ncol = 20)


for(nr in 1:10){
  for(ii in 1:20){
    set.seed(ii)
  ## permutation 
      permute.ind = sample(c(1:34), n.random[nr])
      notpermute.ind = setdiff(c(1:34), permute.ind)
      delete.G.adj = G.adj[c(permute.ind, notpermute.ind), c(permute.ind, notpermute.ind)]

      delete.G.adj[c(1: n.random[nr]) ,c(1: n.random[nr])] = 0
      for(i in 1: (n.random[nr]-1)){
        for(j in (i+1): n.random[nr]){
          delete.G.adj[i,j] = rbinom(1,1,0.1)
          delete.G.adj[j,i] = delete.G.adj[i,j]
        }
      }
  
      p.types = types[c(permute.ind[order(permute.ind)], notpermute.ind)]
      p.Dy.types = Dy.types[c(permute.ind[order(permute.ind)], notpermute.ind), c(permute.ind[order(permute.ind)], notpermute.ind)]
 
      tmp = ifelse(rowSums(delete.G.adj) == 0, 1, rowSums(delete.G.adj))
      delete.G.P = delete.G.adj / tmp
  
      delete.G.dmaps = dmap(delete.G.P, c(1,2,5))
  
      Dx.G1 = as.matrix(dist(delete.G.dmaps[[1]], diag = TRUE, upper = TRUE)) 
      mcorr.1[nr, ii] <- dcor.ttest(Dx.G1, p.Dy.types, distance = TRUE)$p.value
      hhg.1[nr, ii] <- hhg.test(Dx.G1, p.Dy.types, nr.perm = n.perm)$perm.pval.hhg.sl
      mgc.1[nr, ii] <- MGCPermutationTest(Dx.G1, p.Dy.types, rep = n.perm, option = 'mcor')[[1]]
      
      Dx.G2 = as.matrix(dist(delete.G.dmaps[[2]], diag = TRUE, upper = TRUE)) 
      mcorr.2[nr, ii] <- dcor.ttest(Dx.G1, p.Dy.types, distance = TRUE)$p.value
      hhg.2[nr, ii] <- hhg.test(Dx.G1, p.Dy.types, nr.perm = n.perm)$perm.pval.hhg.sl
      mgc.2[nr, ii] <- MGCPermutationTest(Dx.G1, p.Dy.types, rep = n.perm, option = 'mcor')[[1]]
  
      Dx.G5 = as.matrix(dist(delete.G.dmaps[[3]], diag = TRUE, upper = TRUE)) 
      mcorr.5[nr, ii] <- dcor.ttest(Dx.G2, p.Dy.types, distance = TRUE)$p.value
      hhg.5[nr, ii] <- hhg.test(Dx.G2, p.Dy.types, nr.perm = n.perm)$perm.pval.hhg.sl
      mgc.5[nr, ii] <- MGCPermutationTest(Dx.G2, p.Dy.types, rep = n.perm, option = 'mcor')[[1]]
      print(mgc.1[nr, ii])
    
      #### Latent space model
      fh.G[nr, ii] <- min(FH_test2(delete.G.adj, p.types, k.range = c(0:5), n.iter = 300))
  
      print(fh.G[nr,ii])
  }
}


mgc.p = c()
mcorr.p = c()
hhg.p = c()
fh.p = c()
mgc.p = rowMeans(mgc.2)
mcorr.p = rowMeans(mcorr.2)
hhg.p = rowMeans(hhg.2)
fh.p = rowMeans(fh.G)

pdf("figure/tmppvalue.pdf")
par(mfrow = c(1,1), cex.lab = 3, cex.axis = 3,
    mar = c(8,8,3,3), tcl = 0.5)
plot(seq(10,90,10), mgc.p[1:9], col = "red", 
     lty = 1, lwd = 4, ylab = "p-value (avg over 10)",
     ylim = c(0, 0.6), type = "l", mgp = c(6,2,0),
     xlab = "contamination rate (%)")
lines(seq(10,90,10), mcorr.p[1:9], col = "dodgerblue", 
      lty = 2, lwd = 4, yaxis = NULL, type = "l")
lines(seq(10,90,10), hhg.p[1:9], col = "hotpink", 
      lty = 3, lwd = 4, yaxis = NULL, type = "l")
lines(seq(10,90,10),fh.p[1:9], col = "darkgreen", 
      lty = 4, lwd = 4, yaxis = NULL,  type = "l")
legend("topleft",
       c("MGC", "mCorr", "HHG", "FH"), seg.len = 3,
       col = c("red", "dodgerblue", "hotpink", "darkgreen"),
       lty = c(1,2,3,4), lwd = 4, bty = 'n',  xpd = NA, cex = 2)
dev.off()

