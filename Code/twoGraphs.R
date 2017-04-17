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

source("diffmaps.R")
source("dmaps.R")
source("elbowmap.R")
source("FH_factor.R")
source("FH_test.R")
source("getElbows.R")
source("NetworkTest.R")


### setting
dstep = 3
n.perm = 500
n.iter = n.perm

M = 500
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
popn = seq(40, 120, 10) # the number of nodes
pval = matrix(0, 500, length(popn)); twoGraphs = list()

for(N in 1:length(popn)){
  for(i in 1:500){
    set.seed(i)
    generateG = GraphRDPG(popn)
  
    G1 = generateG[[1]]
    G2 = generateG[[2]]
  
    # connect the disconnected components
    G1 = connect_graph(G1)
    G2 = connect_graph(G2)
  
    A1 = as.matrix(get.adjacency(G1))
    A2 = as.matrix(get.adjacency(G2))
  
    P1 = A1  / pmax(rowSums(A1), 1)
    P2 = A2 / pmax(rowSums(A2), 1)
  
    tmp = try( max(getElbows(print.lambda(P1, times = 3)[[1]], plot = FALSE, n = 3)), silent = TRUE)
    if(class(tmp) == "try-error"){
      diffusion.q1 = 1 # when the elbow method does not work well
    }else{
      diffusion.q1  =  min( max(getElbows(print.lambda(P1, times = 3)[[1]], plot = FALSE, n = 3)), popn/2)
    }
  
  
    tmp = try( max(getElbows(print.lambda(P2, times = 3)[[1]], plot = FALSE, n = 3)), silent = TRUE)
    if(class(tmp) == "try-error"){
      diffusion.q2 = 1 # when the elbow method does not work well
    }else{
      diffusion.q2  =  min( max(getElbows(print.lambda(P2, times = 3)[[1]], plot = FALSE, n = 3)), popn/2)
   }
  
    dmap1 = dmap.q(P1, times = 3, q = diffusion.q1)[[1]]
    dmap2 = dmap.q(P2, times = 3, q = diffusion.q2)[[1]] 
  
    Dx = as.matrix(dist(dmap1), diag = TRUE, upper = TRUE)
    Dy = as.matrix(dist(dmap2), diag = TRUE, upper = TRUE)
  
    mgc.result = MGCPermutationTest(Dx, Dy, n.perm, option = 'mcor')[[1]]
    mcorr.result = dcor.ttest(Dx, Dy, distance = TRUE)$p.value
    hhg.result = hhg.test(Dx, Dy, nr.perm = n.perm)$perm.pval.hhg.sl
  
  
    pvals[i,] = c(mgc.result, mcorr.result, hhg.result)
  }
    twoGraphs[[N]] = pvals # matrix with p-values from 500 replicates.
}


save(twoGraphs, file = "../Data/twoGraphs.RData")


for(N in 1:9){
  
  mgc.result = twoGraphs[[N]][,1]
  mcorr.result = twoGraphs[[N]][,2]
  hhg.result = twoGraphs[[N]][,3]
  
  assign( paste("mgc.power", popn[N], sep=""),  rep(0,1))   
  assign( paste("mcorr.power", popn[N], sep=""),  rep(0,1))
  assign( paste("hhg.power", popn[N], sep=""),  rep(0,1))
 
  
  for(i in 1:M){
    assign(paste("mgc.power", popn[N], sep="") , eval(parse(text = paste("mgc.power", popn[N], sep=""))) + (mgc.result[i] <= alpha) / M)
    assign(paste("mcorr.power", popn[N], sep="") , eval(parse(text = paste("mcorr.power", popn[N], sep=""))) + (mcorr.result[i] <= alpha) / M)
    assign(paste("hhg.power", popn[N], sep="") , eval(parse(text = paste("hhg.power", popn[N], sep=""))) + (hhg.result[i] <= alpha) / M)
  }
  
}

mgc.power = c(mgc.power40, mgc.power50, mgc.power60, mgc.power70, mgc.power80,
              mgc.power90, mgc.power100, mgc.power110, mgc.power120)
mcorr.power = c(mcorr.power40, mcorr.power50, mcorr.power60, mcorr.power70, mcorr.power80,
                mcorr.power90, mcorr.power100, mcorr.power110, mcorr.power120)
hhg.power =  c(hhg.power40, hhg.power50, hhg.power60, hhg.power70, hhg.power80,
               hhg.power90, hhg.power100, hhg.power110, hhg.power120)

## Make Figure
pdf("../Figure/Graphs.pdf", width = 12, height = 6)
par(mfrow = c(1,1), cex.lab = 5, cex.axis = 3,
    mar = c(8,10,3,20), tcl = 0.5)
plot(seq(40, 120, 10), mgc.power, col = "red", 
     lty = 1, lwd = 5, ylab = "Power",
     ylim = c(0,1.1), type = "l", mgp = c(6,2,0),
     xlab = "number of nodes", yaxt = 'n', frame.plot = FALSE)
axis(side = 2, at = c(0.0, 0.25, 0.5,  0.75, 1.0), 
     labels = c(0.0, 0.25, 0.5,  0.75, 1.0), 
     tck = 0.05)
lines(seq(40, 120, 10), mcorr.power, col = "dodgerblue", 
      lty =2, lwd = 5,  type = "l")
lines(seq(40, 120, 10), hhg.power, col = "lightsalmon4", 
      lty =3, lwd = 5,  type = "l")
legend("topright", inset=c(-0.6, 0.5), 
       c(expression(MGC), expression(mCorr), expression(HHG)),
       col = c("red", "dodgerblue", "lightsalmon4"), seg.len = 3,
       lty = c(1,2,3), lwd = 4, bty = 'n', cex = 2.5, xpd = NA)
dev.off()


