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
  U = runif(popn, 0, 1)
  
  A1 = matrix(0, popn, popn); A2 = matrix(0, popn, popn) 
  
  for(i in 1: (popn-1)){
    for(j in (i+1):popn){
      p1 = (U[i]*U[j])
      A1[i,j] = rbinom(1,1,p1)
      A1[j,i] = A1[i,j]
      
      p2 = ((U[i])^2*(U[j])^2)
      A2[i,j] = rbinom(1,1,p2)
      A2[j,i] = A2[i,j]
      
    }
  }
  
  G1 = graph.adjacency(A1, "undirected")
  V(G1)$U = U
  
  G2 = graph.adjacency(A2, "undirected")
  V(G2)$U = U^2
  
  return(list(G1, G2)) 
}


#### implement
popn = c(seq(10, 30, 2), 40) # the number of nodes
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
  
    mgc.result = MGCPermutationTest(Dy, Dy, n.perm, option = 'mcor')[[1]]
    mcorr.result = dcor.ttest(Dx, Dy, distance = TRUE)$p.value
    hhg.result = hhg.test(Dx, Dy, nr.perm = n.perm)$perm.pval.hhg.sl
  
    Dx = as.matrix(dist(A1), diag = TRUE, upper = TRUE)
    Dy = as.matrix(dist(A2), diag = TRUE, upper = TRUE)
  
    adj.mgc.result = MGCPermutationTest(Dy, Dy, n.perm, option = 'mcor')[[1]]
    adj.mcorr.result = dcor.ttest(Dx, Dy, distance = TRUE)$p.value
    adj.hhg.result = hhg.test(Dx, Dy, nr.perm = n.perm)$perm.pval.hhg.sl 
  
    pvals[i,] = c(mgc.result, mcorr.result, hhg.result, 
              adj.mgc.result, adj.mcorr.result, adj.hhg.result)
  }
    twoGraphs[[N]] = pvals
    twoGraphs[[N]][is.na(twoGraphs[[N]])] = 1 #  we rarely have NA due limited precision
}


save(twoGraphs, file = "../Data/twoGraphs.RData")

for(N in 1:12){
  
  mgc.result = twoGraphs[[N]][,1]
  mcorr.result = twoGraphs[[N]][,2]
  hhg.result = twoGraphs[[N]][,3]
  
  adj.mgc.result = twoGraphs[[N]][,4]
  adj.mcorr.result = twoGraphs[[N]][,5]
  adj.hhg.result = twoGraphs[[N]][,6]
  
  assign( paste("mgc.power", popn[N], sep=""),  rep(0,1))   
  assign( paste("mcorr.power", popn[N], sep=""),  rep(0,1))
  assign( paste("hhg.power", popn[N], sep=""),  rep(0,1))
  assign( paste("adj.mgc.power", popn[N], sep=""),  rep(0,1))   
  assign( paste("adj.mcorr.power", popn[N], sep=""),  rep(0,1))
  assign( paste("adj.hhg.power", popn[N], sep=""),  rep(0,1))
  
  
  for(i in 1:M){
    assign(paste("mgc.power", popn[N], sep="") , eval(parse(text = paste("mgc.power", popn[N], sep=""))) + (mgc.result[i] <= alpha) / M)
    assign(paste("mcorr.power", popn[N], sep="") , eval(parse(text = paste("mcorr.power", popn[N], sep=""))) + (mcorr.result[i] <= alpha) / M)
    assign(paste("hhg.power", popn[N], sep="") , eval(parse(text = paste("hhg.power", popn[N], sep=""))) + (hhg.result[i] <= alpha) / M)
    
    assign(paste("adj.mgc.power", popn[N], sep="") , eval(parse(text = paste("adj.gc.power", popn[N], sep=""))) + (adj.mgc.result[i] <= alpha) / M)
    assign(paste("adj.mcorr.power", popn[N], sep="") , eval(parse(text = paste("adj.mcorr.power", popn[N], sep=""))) + (adj.mcorr.result[i] <= alpha) / M)
    assign(paste("adj.hhg.power", popn[N], sep="") , eval(parse(text = paste("adj.hhg.power", popn[N], sep=""))) + (adj.hhg.result[i] <= alpha) / M)
  }
  
}

mgc.power = c(mgc.power10, mgc.power12, mgc.power14, mgc.power16, mgc.power18, mgc.power20,
              mgc.power22, mgc.power24, mgc.power26, mgc.power28, mgc.power30, mgc.power40)
mcorr.power = c(mcorr.power10, mcorr.power12, mcorr.power14, mcorr.power16, mcorr.power18, mcorr.power20,
                mcorr.power22, mcorr.power24, mcorr.power26, mcorr.power28, mcorr.power30, mcorr.power40)
hhg.power = c(hhg.power10, hhg.power12, hhg.power14, hhg.power16, hhg.power18, hhg.power20,
              hhg.power22, hhg.power24, hhg.power26, hhg.power28, hhg.power30, hhg.power40)

adj.mgc.power = c(adj.mgc.power10, adj.mgc.power12, adj.mgc.power14, adj.mgc.power16, adj.mgc.power18, adj.mgc.power20,
                  adj.mgc.power22, adj.mgc.power24, adj.mgc.power26, adj.mgc.power28, adj.mgc.power30, adj.mgc.power40)
adj.mcorr.power = c(adj.mcorr.power10, adj.mcorr.power12, adj.mcorr.power14, adj.mcorr.power16, adj.mcorr.power18, adj.mcorr.power20,
                    adj.mcorr.power22, adj.mcorr.power24, adj.mcorr.power26, adj.mcorr.power28, adj.mcorr.power30, adj.mcorr.power40)
adj.hhg.power = c(adj.hhg.power10, adj.hhg.power12, adj.hhg.power14, adj.hhg.power16, adj.hhg.power18, adj.hhg.power20,
                  adj.hhg.power22, adj.hhg.power24, adj.hhg.power26, adj.hhg.power28, adj.hhg.power30, adj.hhg.power40)


## Make Figure
pdf("../Figure/twoGraphs1.pdf", width = 12, height = 6)
par(mfrow = c(1,1), cex.lab = 5, cex.axis = 3,
    mar = c(8,10,3,20), tcl = 0.5)
plot(c(seq(10, 30, 2),40), mgc.power, col = "red", 
     lty = 1, lwd = 5, ylab = "Power",
     ylim = c(0,1.1), type = "l", mgp = c(6,2,0),
     xlab = "number of nodes", yaxt = 'n', frame.plot = FALSE)
axis(side = 2, at = c(0.0, 0.25, 0.5,  0.75, 1.0), 
     labels = c(0.0, 0.25, 0.5,  0.75, 1.0), 
     tck = 0.05)
lines(c(seq(10, 30, 2),40), mcorr.power, col = "dodgerblue", 
      lty =2, lwd = 5,  type = "l")
lines(c(seq(10, 30, 2),40), hhg.power, col = "lightsalmon4", 
      lty =3, lwd = 5,  type = "l")
legend("topright", inset=c(-0.6, 0.5), 
       c(expression(MGC), expression(mCorr), expression(HHG)),
       col = c("red", "dodgerblue", "lightsalmon4"), seg.len = 3,
       lty = c(1,2,3), lwd = 4, bty = 'n', cex = 2.5, xpd = NA)
dev.off()

