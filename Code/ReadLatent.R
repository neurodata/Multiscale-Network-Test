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
n.sample = seq(20, 100, 10)

simpleRDPG = function(popn){
    U = runif(popn, 0, 1)
    
    A = matrix(0, popn, popn)
    
    for(i in 1: (popn-1)){
        for(j in (i+1):popn){
            p = (U[i]*U[j])
            A[i,j] = rbinom(1,1,p)
            A[j,i] = A[i,j]
        }
    }
    
    G = graph.adjacency(A, "undirected")
    V(G)$U = U
    
    return(G) 
}

###############################
RDPG = list(); mgc.results = c(); fh.results = c(); fmgc.results = c()

for(N in 1:9){
    # for each sample size
    for(i in 1:M){
        set.seed(i)
        
      G = simpleRDPG(n.sample[N])
      G = connect_graph(G)
      U = V(G)$U # latent vector
      X = rnorm(popn, U, 0.5) # nodal attributes
      A = as.matrix(get.adjacency(G))
      P = A  / pmax(rowSums(A), 1)
      diffusion.q  =  min( max(getElbows(print.lambda(P, times = 3)[[1]], plot = FALSE, n = 3)), n.sample[N]-1)
      mgc.results[i] =  NetworkTest.q(G, X, option = 1, diffusion = TRUE, dstep = 3, n.perm = n.perm, q = diffusion.q)[[1]][[1]]
          
  
      SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
      fh.k = min( max(getElbows(SVD.A$d, n = 3, plot = FALSE)), n.sample[N]-1)
      fh.results[i] = FH_test(A, X, k.range = fh.k, n.iter = n.iter)
      

      factors = FH_factor(A, k.range = fh.k)
      Dx = as.matrix(dist(factors[[1]]), diag = TRUE, upper = TRUE)
      Dy = as.matrix(dist(X), diag = TRUE, upper = TRUE)
      fmgc.results[i] = MGCPermutationTest(Dx, Dy, rep = n.perm, option = 'mcor')[[1]]

    }
    
    RDPG[[N]] = cbind(mgc.results, fh.results, fmgc.results)
}


save(RDPG, file = "../Data/RDPG.RData")

for(N in 1:9){
  
  Smgc.results = RDPG[[N]][,1]
  fh.results = RDPG[[N]][,2]
  Fmgc.results = RDPG[[N]][,3]
  
  assign( paste("Smgc.power", n.sample[N], sep=""),  rep(0,1))   
  assign( paste("fh.power", n.sample[N], sep=""),  rep(0,1))
  assign( paste("Fmgc.power", n.sample[N], sep=""),  rep(0,1))
  
  for(i in 1:M){
    assign(paste("Smgc.power", n.sample[N], sep="") , eval(parse(text = paste("Smgc.power", n.sample[N], sep=""))) + (Smgc.results[i] <= alpha) / M)
    assign(paste("fh.power", n.sample[N], sep="") , eval(parse(text = paste("fh.power", n.sample[N], sep=""))) + (fh.results[i] <= alpha) / M)
    assign(paste("Fmgc.power", n.sample[N], sep="") , eval(parse(text = paste("Fmgc.power", n.sample[N], sep=""))) + (Fmgc.results[i] <= alpha) / M)
  }
}

## Summarize power
Smgc.power = c(Smgc.power20, Smgc.power30, Smgc.power40, Smgc.power50,
               Smgc.power60, Smgc.power70, Smgc.power80, Smgc.power90, Smgc.power100)
fh.power = c(fh.power20, fh.power30, fh.power40, fh.power50,
             fh.power60, fh.power70, fh.power80, fh.power90, fh.power100)
Fmgc.power = c(Fmgc.power20, Fmgc.power30, Fmgc.power40, Fmgc.power50,
                Fmgc.power60, Fmgc.power70, Fmgc.power80, Fmgc.power90, Fmgc.power100)

## Make Figure
pdf("../Figure/elbow3_RDPG.pdf", width = 15, height = 6)
par(mfrow = c(1,1), cex.lab = 4, cex.axis = 3,
    mar = c(8,10,3,20), tcl = 0.5)
plot(seq(20, 100, 10), Smgc.power, col = "red", 
     lty = 1, lwd = 5, ylab = "Power",
     ylim = c(0,1), type = "l", mgp = c(6,2,0),
     xlab = "Sample Size", yaxt = 'n')
axis(side = 2, at = c(0.0, 0.25, 0.5,  0.75, 1.0), 
     labels = c(0.0, 0.25, 0.5,  0.75, 1.0), 
     tck = 0.05)
lines(seq(20, 100, 10), Fmgc.power, col = "midnightblue", 
      lty =2, lwd = 5,  type = "l")
lines(seq(20, 100, 10), fh.power, col = "darkgreen", 
      lty =3, lwd = 5,  type = "l")
legend("topright", inset=c(-0.4, 0.5), 
       c(expression(MGC %.% DM), expression(MGC %.% LF), "FH Test"),
       col = c("red", "midnightblue", "darkgreen"), seg.len = 3,
       lty = c(1,2,3), lwd = 4, bty = 'n', cex = 2, xpd = NA)
abline(v = 0.2, lwd = 2, col = "black")
dev.off()
