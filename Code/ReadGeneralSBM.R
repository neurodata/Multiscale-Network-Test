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
source("dmaps.R")
source("elbowmap.R")
source("FH_factor.R")
source("FH_test.R")
source("getElbows.R")
source("NetworkTest.R")

popn = 100; n.group = 3; M = 500
p = 0.5; q = 0.2
r = seq(0.05, 0.6, 0.05)

dstep <- 3
n.perm <- 500
n.iter <- n.perm

mgc.results = c(); mcorr.results = c(); hhg.results = c(); fh.results = c()
GeneralSBM = list()

for(N in 1:12){
  # different block probabilities
  for(i in 1:M){
    set.seed(i)

    G = notsimple3block(popn, p, q, r[N])
    X = V(G)$outcome
    A = as.matrix(get.adjacency(G))
    P =  A / pmax(rowSums(A),1)
    diffusion.q  =  min( max(getElbows(print.lambda(P, times = 3)[[1]], plot = FALSE, n = 3, threshold = 0)), popn-1)

    mgc.results[i] =  NetworkTest.q(G, X, option = 1, diffusion = TRUE, dstep = 3, n.perm = n.perm, q = diffusion.q)[[1]][[1]]
    mcorr.results[i] =  NetworkTest.q(G, X, option = 2, diffusion = TRUE, dstep = 3, n.perm = n.perm, q = diffusion.q)
    hhg.results[i] =  NetworkTest.q(G, X, option = 3, diffusion = TRUE, dstep = 3, n.perm = n.perm, q = diffusion.q)
        
    # estimate the matrix M based on SVD
    SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
    # FH
    fh.k = min( max(getElbows(SVD.A$d, n = 3, plot = FALSE)), popn-1)
    fh.results[i] = FH_test(A, X, k.range = fh.k, n.iter = n.iter)
  }

  GeneralSBM[[N]] = cbind(mgc.results, mcorr.results, hhg.results, fh.results)

}

save(GeneralSBM, file = "../Data/GeneralSBM.RData")


theta = seq(5, 60, 5)
alpha = 0.05
M = 500
theta[1] = "05"

for(N in 1:12){

  mgc.results = GeneralSBM[[N]][,1]
  mcorr.results = GeneralSBM[[N]][,2]
  hhg.results = GeneralSBM[[N]][,3]
  fh.results = GeneralSBM[[N]][,4]
  
  assign( paste("mgc.power", theta[N], sep=""),  rep(0,1))   
  assign( paste("mcorr.power", theta[N], sep=""),  rep(0,1)) 
  assign( paste("hhg.power", theta[N], sep=""),  rep(0,1))
  assign( paste("fh.power", theta[N], sep=""),  rep(0,1))
  
  # calculate power
  for(i in 1:M){
    assign(paste("mgc.power", theta[N], sep="") , eval(parse(text = paste("mgc.power", theta[N], sep=""))) + (mgc.results[i] <= alpha) / M)
    assign(paste("mcorr.power", theta[N], sep="") , eval(parse(text = paste("mcorr.power", theta[N], sep=""))) + (mcorr.results[i] <= alpha) / M)
    assign(paste("hhg.power", theta[N],  sep="") , eval(parse(text = paste("hhg.power", theta[N], sep=""))) + (hhg.results[i] <= alpha) / M)
    assign(paste("fh.power", theta[N], sep="") , eval(parse(text = paste("fh.power", theta[N], sep=""))) + (fh.results[i] <= alpha) / M)
  }

}

### Summarize power
mgc.power = c(mgc.power05, mgc.power10, mgc.power15, mgc.power20,
               mgc.power25, mgc.power30, mgc.power35, mgc.power40,
               mgc.power45, mgc.power50, mgc.power55, mgc.power60)
mcorr.power = c(mcorr.power05, mcorr.power10, mcorr.power15, mcorr.power20,
               mcorr.power25, mcorr.power30, mcorr.power35, mcorr.power40,
               mcorr.power45, mcorr.power50, mcorr.power55, mcorr.power60)
hhg.power = c(hhg.power05, hhg.power.q10, hhg.power.q15, hhg.power.q20,
              hhg.power25, hhg.power.q30, hhg.power.q35, hhg.power.q40,
              hhg.power45, hhg.power.q50, hhg.power.q55, hhg.power.q60)
fh.power = c(fh.power05, fh.power10, fh.power15, fh.power20,
             fh.power25, fh.power30, fh.power35, fh.power40,
              fh.power45, fh.power50, fh.power55, fh.power60)

### Make plot
pdf("../Figure/monoelbow3_t3.pdf", width = 15, height = 6)
par(mfrow = c(1,1), cex.lab = 5, cex.axis = 3,
    mar = c(8,10,3,20), tcl = 0.5)
plot( seq(0.05, 0.60, 0.05),  mgc.power, col = "red", 
      lty = 1, lwd = 5, ylab = "Power",
      ylim = c(0, 1), type = "l", mgp = c(6,2,0),
      xlab = expression(theta), yaxt = 'n')
axis(side = 2, at = c(0.0, 0.25, 0.5,  0.75, 1.0), 
     labels = c(0.0, 0.25, 0.5,  0.75, 1.0), 
     tck = 0.05)
lines( seq(0.05, 0.60, 0.05),  mcorr.power, col = "dodgerblue", 
       lty = 2, lwd = 5,  type = "l")
lines( seq(0.05, 0.60, 0.05),  hhg.power, col = "lightsalmon4", 
       lty = 3, lwd = 5,  type = "l")
lines( seq(0.05, 0.60, 0.05),  fh.power, col = "darkgreen", 
       lty =4, lwd = 5,  type = "l")
legend("topright", inset=c(-0.4, 0.5), 
       c( expression(MGC %.% DM), expression(mCorr %.% DM), 
          expression(HHG %.% DM), "FH Test"), 
       col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"), seg.len = 3,
       lty = c(1,2,3,4), lwd = 4, bty = 'n', cex = 2, xpd = NA)
abline(v = 0.2, lwd = 2, col = "black")
dev.off()

