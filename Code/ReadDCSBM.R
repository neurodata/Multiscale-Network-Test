library(gplots)
library(igraph)
library(MASS)
library(HHG)
library(lattice)
library(DTMCPack)
library(Matrix)
library(energy)
library(Rlab)
library(amen)
library(doParallel)
library(ecodist)
library(SDMTools)

source("notsimple3block.R")
source("diffmaps.R")
source("FH_test.R")
source("getElbows.R")
source("NetworkTest.R")
source("optimalT.R")
###
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

### setting
tau = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)

popn = 200
n.perm = 500; n.iter = n.perm
M = 100; alpha = 0.05
dcSBM = list(); results = list()


for(N in 1:6){
    # for each tau
    for(i in 1:M){
      set.seed(i)
      G  = varyA(popn, p = 0.2, q = 0.05, tau = tau[N])
      if(!is.connected(G)) G = connect_graph(G)
      A = as.matrix(get.adjacency(G))
      X = V(G)$outcome
      
      mgc.result =  NetworkTest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
      adj.result =  NetworkTest.diffusion.stat(G, X, option = 1, diffusion = FALSE, t.range = c(0:10), n.perm = n.perm)
      # estimate the matrix M based on SVD
      SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
      # FH
      fh.k = min( max(getElbows(abs(SVD.A$d), n = 2, plot = FALSE)), popn-10)
      fh.result = FH_test(A, X, k.range = fh.k, n.iter = n.iter)
      factors = FH_factor(A, k.range = fh.k)
      Dx = as.matrix(dist(factors[[1]]), diag = TRUE, upper = TRUE)
      Dy = as.matrix(dist(X), diag = TRUE, upper = TRUE)
      fmgc.result = MGCPermutationTest(Dx, Dy, rep = n.perm, option = 'mcor')
  
      results[[i]] = list(mgc.result, adj.result, fmgc.result, fh.result))  
    }
    dcSBM[[N]] = results
}

## save the pvalues
save(dcSBM, file = "../Data/dcSBM.RData")

multi_tau0 = dcSBM[[1]]; multi_tau2 = dcSBM[[2]]
multi_tau4 = dcSBM[[3]]; multi_tau6 = dcSBM[[4]]
multi_tau8 = dcSBM[[5]]; multi_tau10 = dcSBM[[6]]


df.tau0 = c(); adj.tau0 = c(); lf.tau0 = c(); fh.tau0 = c()
for(i in 1:length(multi_tau0)){
  df.tau0[i] = print.stat.optimal(multi_tau0[[i]][[1]], 4)$pvalue
  adj.tau0[i] = multi_tau0[[i]][[2]]$pMGC
  lf.tau0[i] = multi_tau0[[i]][[3]]$pMGC
  fh.tau0[i] = multi_tau0[[i]][[4]]
}

df.tau2 = c(); adj.tau2 = c(); lf.tau2 = c(); fh.tau2 = c()
for(i in 1:length(multi_tau2)){
  df.tau2[i] = print.stat.optimal(multi_tau2[[i]][[1]], 4)$pvalue
  adj.tau2[i] = multi_tau2[[i]][[2]]$pMGC
  lf.tau2[i] = multi_tau2[[i]][[3]]$pMGC
  fh.tau2[i] = multi_tau2[[i]][[4]]
}

df.tau4 = c(); adj.tau4 = c(); lf.tau4 = c(); fh.tau4 = c()
for(i in 1:length(multi_tau4)){
  df.tau4[i] = print.stat.optimal(multi_tau4[[i]][[1]], 4)$pvalue
  adj.tau4[i] = multi_tau4[[i]][[2]]$pMGC
  lf.tau4[i] = multi_tau4[[i]][[3]]$pMGC
  fh.tau4[i] = multi_tau4[[i]][[4]]
}

df.tau6 = c(); adj.tau6 = c(); lf.tau6 = c(); fh.tau6 = c()
for(i in 1:length(multi_tau6)){
  df.tau6[i] = print.stat.optimal(multi_tau6[[i]][[1]], 4)$pvalue
  adj.tau6[i] = multi_tau6[[i]][[2]]$pMGC
  lf.tau6[i] = multi_tau6[[i]][[3]]$pMGC
  fh.tau6[i] = multi_tau6[[i]][[4]]
}

df.tau8 = c(); adj.tau8 = c(); lf.tau8 = c(); fh.tau8 = c()
for(i in 1:length(multi_tau8)){
  df.tau8[i] = print.stat.optimal(multi_tau8[[i]][[1]], 4)$pvalue
  adj.tau8[i] = multi_tau8[[i]][[2]]$pMGC
  lf.tau8[i] = multi_tau8[[i]][[3]]$pMGC
  fh.tau8[i] = multi_tau8[[i]][[4]]
}

df.tau10 = c(); adj.tau10 = c(); lf.tau10 = c(); fh.tau10 = c()
for(i in 1:length(multi_tau10)){
  df.tau10[i] = print.stat.optimal(multi_tau10[[i]][[1]], 4)$pvalue
  adj.tau10[i] = multi_tau10[[i]][[2]]$pMGC
  lf.tau10[i] = multi_tau10[[i]][[3]]$pMGC
  fh.tau10[i] = multi_tau10[[i]][[4]]
}

df.power = colMeans(cbind(df.tau0, df.tau2, df.tau4, df.tau6, df.tau8, df.tau10) <= 0.05)
adj.power = colMeans(cbind(adj.tau0, adj.tau2, adj.tau4, adj.tau6, adj.tau8, adj.tau10) <= 0.05)
lf.power = colMeans(cbind(lf.tau0, lf.tau2, lf.tau4, lf.tau6, lf.tau8, lf.tau10) <= 0.05)
fh.power = colMeans(cbind(fh.tau0, fh.tau2, fh.tau4, fh.tau6, fh.tau8, fh.tau10) <= 0.05)

### Make plot
pdf("../Figure/multi_DCSBM.pdf", width = 15, height = 6)
par(mfrow = c(1,1), cex.lab = 3, cex.axis = 2,
    mar = c(6,8,3,17), tcl = 0.5)
plot(seq(0,1,0.2),  df.power, col = "red", 
     lty = 1, lwd = 4, ylab = "Power",
     ylim = c(0, 1.0), type = "l", mgp = c(5,2,0),
     xlab = expression(paste(tau)), yaxt = "n")
axis(side = 2, at = c(0.0, 0.25, 0.5,  0.75, 1.0), 
     labels = c(0.0, 0.25, 0.5,  0.75, 1.0), 
     tck = 0.05)
lines(seq(0,1,0.2),  adj.power, col = "mediumblue", 
      lty = 2, lwd = 4, yaxis = NULL, type = "l")
lines(seq(0,1,0.2),  lf.power, col = "gold4", 
      lty = 5, lwd = 4, yaxis = NULL,  type = "l")
lines(seq(0,1,0.2),  fh.power, col = "darkgreen", 
      lty = 4, lwd = 4, yaxis = NULL,  type = "l")
legend("topright", inset=c(-0.32, 0.5), 
       c(expression(MGC %.% DM), expression(MGC %.% AM),
         expression(paste( MGC %.% LF)), "FH Test"), seg.len = 3,
       col = c("red", "mediumblue", "gold4", "darkgreen"),
       lty = c(1,2,5,4), lwd = 4, bty = 'n',  xpd = NA, cex = 2)
dev.off()
