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

source("FH_factor.R")
source("FH_test.R")
source("getElbows.R")
source("NetworkTest.R")
source("diffmaps.R")
source("optimalT.R")

### setting
n.perm = 500
n.iter = n.perm
alpha = 0.05

### generate two graphs from RDPG but different link functions
GraphRDPG = function(popn){

  # generate multi-dimensional latent variable U
  # and one-dimensional latent variable W
  U = matrix(0, nrow = popn, ncol = 5)
  for(j in 1:5){
    U[,j] = runif(popn, 0, 1)
  }
  A1 = matrix(0, popn, popn); A2 = matrix(0, popn, popn) 
 
    for(i in 1: (popn-1)){
      for(j in (i+1):popn){
          p1 = sum(U[i,]*U[j,])/5
          A1[i,j] = rbinom(1,1,p1)
          A1[j,i] = A1[i,j]

          p2 =  sum((1-U[i,1])^2*(1-U[j,1])^2)/5
          A2[i,j] = rbinom(1,1,p2)
          A2[j,i] = A2[i,j]

      }
    }

    G1 = graph.adjacency(A1, "undirected")
    G2 = graph.adjacency(A2, "undirected")
   

  return(list(G1, G2)) 
}


#### implement
popn = seq(50, 150, 20) # the number of nodes
results = list(); twoGraphs = list();

for(N in 1:length(popn)){
  
  for(i in 1:100){
    set.seed(i)
    generateG = GraphRDPG(popn)

    G1 = generateG[[1]]
    G2 = generateG[[2]]
    
    G1 = connect_graph(G1)
    G2 = connect_graph(G2)
 
    mgc.result =  TwoGraphs.diffusion.stat(G1, G2, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
    mcorr.result =  TwoGraphs.diffusion.stat(G1, G2, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
    hhg.result =  TwoGraphs.diffusion.stat(G1, G2, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
        
    results[[i]] = list(mgc.result, mcorr.result, hhg.result) 
  }
    
    twoGraphs[[N]] = results

}

save(twoGraphs, file = "../Data/twoGraphs.RData")
############################################
multi_graph50 = twoGraphs[[1]]; multi_graph70 = twoGraphs[[2]]
multi_graph90 = twoGraphs[[3]]; multi_graph110 = twoGraphs[[4]]
multi_graph130 = twoGraphs[[5]]; multi_graph150 = twoGraphs[[6]]

pval.mgc = c(); pval.mcorr = c(); pval.hhg = c()
for(i in 1:length(multi_graph50)){
  pval.mgc[i] = print.stat.optimal(multi_graph50[[i]][[1]], 4)$pvalue
  pval.mcorr[i] = print.stat.optimal(multi_graph50[[i]][[2]], 4)$pvalue
  pval.hhg[i] = print.stat.optimal(multi_graph50[[i]][[3]], 4)$pvalue
}
power50 = colMeans(cbind(pval.mgc, pval.mcorr, pval.hhg) <= 0.05)

for(i in 1:length(multi_graph70)){
  pval.mgc[i] = print.stat.optimal(multi_graph70[[i]][[1]], 4)$pvalue
  pval.mcorr[i] = print.stat.optimal(multi_graph70[[i]][[2]], 4)$pvalue
  pval.hhg[i] = print.stat.optimal(multi_graph70[[i]][[3]], 4)$pvalue
}
power70 = colMeans(cbind(pval.mgc, pval.mcorr, pval.hhg) <= 0.05)

for(i in 1:length(multi_graph90)){
  pval.mgc[i] = print.stat.optimal(multi_graph90[[i]][[1]], 4)$pvalue
  pval.mcorr[i] = print.stat.optimal(multi_graph90[[i]][[2]], 4)$pvalue
  pval.hhg[i] = print.stat.optimal(multi_graph90[[i]][[3]], 4)$pvalue
}
power90 = colMeans(cbind(pval.mgc, pval.mcorr, pval.hhg) <= 0.05)

for(i in 1:length(multi_graph110)){
  pval.mgc[i] = print.stat.optimal(multi_graph110[[i]][[1]], 4)$pvalue
  pval.mcorr[i] = print.stat.optimal(multi_graph110[[i]][[2]], 4)$pvalue
  pval.hhg[i] = print.stat.optimal(multi_graph110[[i]][[3]], 4)$pvalue
}
power110 = colMeans(cbind(pval.mgc, pval.mcorr, pval.hhg) <= 0.05)

for(i in 1:length(multi_graph130)){
  pval.mgc[i] = print.stat.optimal(multi_graph130[[i]][[1]], 4)$pvalue
  pval.mcorr[i] = print.stat.optimal(multi_graph130[[i]][[2]], 4)$pvalue
  pval.hhg[i] = print.stat.optimal(multi_graph130[[i]][[3]], 4)$pvalue
}
power130 = colMeans(cbind(pval.mgc, pval.mcorr, pval.hhg) <= 0.05)

for(i in 1:length(multi_graph150)){
  pval.mgc[i] = print.stat.optimal(multi_graph150[[i]][[1]], 4)$pvalue
  pval.mcorr[i] = print.stat.optimal(multi_graph150[[i]][[2]], 4)$pvalue
  pval.hhg[i] = print.stat.optimal(multi_graph150[[i]][[3]], 4)$pvalue
}
power150 = colMeans(cbind(pval.mgc, pval.mcorr, pval.hhg) <= 0.05)


mat = rbind(power50, power70, power90, power110, power130, power150)
mgc.power = mat[,1]; mcorr.power = mat[,2];  hhg.power = mat[,3]

## Make Figure
pdf("../Figure/multi_graphs.pdf", width = 13, height = 6)
par(mfrow = c(1,1), cex.lab = 3, cex.axis = 2,
    mar = c(5,8,3,16), tcl = 0.5)
plot(c(50, 70, 90, 110, 130, 150), mgc.power, col = "red",
     lty = 1, lwd = 5, ylab = "Power",
     ylim = c(0,1), type = "l", mgp = c(4,2,0),
     xlab = "number of nodes", yaxt = 'n', xaxt = 'n')
axis(side = 2, at = seq(0, 1, 0.2),
     labels = seq(0, 1, 0.2),
     tck = 0.05)
axis(side = 1, at = seq(50, 150, 20),
     labels = seq(50, 150, 20),
     tck = 0.05)
lines(c(50, 70, 90, 110, 130, 150), mcorr.power, col = "dodgerblue",
      lty =2, lwd = 5,  type = "l")
lines(c(50, 70, 90, 110, 130, 150), hhg.power, col = "lightsalmon4",
      lty =5, lwd = 5,  type = "l")
legend("topright", inset=c(-0.3, 0.5),
       c(expression(MGC), expression(dCorr), expression(HHG)),
       col = c("red", "dodgerblue", "lightsalmon4"), seg.len = 3,
       lty = c(1,2,5), lwd = 4, bty = 'n', cex = 2, xpd = NA)
dev.off()
