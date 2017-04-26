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
source("varyA.R")

### setting
tau = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)

popn = 200
dstep = 3
n.perm = 500
n.iter = n.perm

M = 500
rho = c(0,2,4,6,8,10)
alpha = 0.05

dcSBM = list()
mgc.result = c(); adj.result = c(); fh.result = c(); fmgc.results = c()

for(N in 1:6){
    # for each tau
    for(i in 1:M){
        set.seed(i)
        
        G  = varyA(popn, p = 0.2, q = 0.05, tau = tau[N])
       
        A = as.matrix(get.adjacency(G))
        P = A / pmax(rowSums(A), 1)
        X = V(G)$outcome
        
        # upper bound of dimension is set to 100.
        diffusion.q = min( max(getElbows(print.lambda(P, times = 3)[[1]], plot = FALSE, n = 3, threshold = 0)), 100)
        mgc.result[i] =  NetworkTest.q(G, X, option = 1, diffusion = TRUE, dstep = 3, n.perm = n.perm, q = diffusion.q)[[1]][[1]]
        adj.result[i] =  NetworkTest.q(G, X, option = 1, diffusion = FALSE, dstep = 3, n.perm = n.perm, q = diffusion.q)
        
        SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
        # upper bound of dimension is set to 100
        fh.k = min( max(getElbows(SVD.A$d, n = 3, plot = FALSE)), 100)
        
        fh.result[i] = FH_test(A, X, k.range = fh.k, n.iter = n.iter)
        factors = FH_factor(A, k.range = fh.k)
        Dx = as.matrix(dist(factors[[1]]), diag = TRUE, upper = TRUE)
        Dy = as.matrix(dist(X), diag = TRUE, upper = TRUE)
        fmgc.result[i] = MGCPermutationTest(Dx, Dy, rep = n.perm, option = 'mcor')[[1]]
    }
    dcSBM[[N]] = cbind(mgc.result, adj.result, fh.result, fmgc.result)
}

## save the pvalues
save(dcSBM, file = "../Data/dcSBM.RData")


for(N in 1:6){

  df.results = dcSBM[[N]][,1]
  adj.results = dcSBM[[N]][,2]
  fmgc.results = dcSBM[[N]][,3]
  fh.results = dcSBM[[N]][,4]
    
  assign( paste("df.power", rho[N], sep=""),  rep(0,1))   
  assign( paste("adj.power", rho[N], sep=""),  rep(0,1)) 
  assign( paste("lf.power", rho[N], sep=""),  rep(0,1)) 
  assign( paste("fh.power", rho[N], sep=""),  rep(0,1))
  
  for(i in 1:M){
    assign(paste("df.power", rho[N], sep="") , eval(parse(text = paste("df.power", rho[N], sep=""))) + (df.results[i] <= alpha) / M)
    assign(paste("adj.power", rho[N], sep="") , eval(parse(text = paste("adj.power", rho[N], sep=""))) + (adj.results[i] <= alpha) / M)
    assign(paste("lf.power", rho[N], sep="") , eval(parse(text = paste("lf.power", rho[N], sep=""))) + (fmgc.results[i] <= alpha) / M)
    assign(paste("fh.power", rho[N], sep="") , eval(parse(text = paste("fh.power", rho[N], sep=""))) + (fh.results[i] <= alpha) / M)
  }

}

df.power = c(df.power0, df.power2, df.power4, df.power6, df.power8, df.power10)
adj.power = c(adj.power0, adj.power2, adj.power4, adj.power6, adj.power8, adj.power10)
lf.power = c(lf.power0, lf.power2, lf.power4, lf.power6, lf.power8, lf.power10)
fh.power = c(fh.power0, fh.power2, fh.power4, fh.power6, fh.power8, fh.power10)


### Make plot
pdf("../Figure/elbow3_dcSBM.pdf", width = 15, height = 6)
par(mfrow = c(1,1), cex.lab = 5, cex.axis = 3,
    mar = c(8,10,3,20), tcl = 0.5)
plot(seq(0,1,0.2),  df.power, col = "red", 
     lty = 1, lwd = 4, ylab = "Power",
     ylim = c(0, 1.0), type = "l", mgp = c(6,2,0),
     xlab = expression(paste(tau)), yaxt = "n")
axis(side = 2, at = c(0.0, 0.25, 0.5,  0.75, 1.0), 
     labels = c(0.0, 0.25, 0.5,  0.75, 1.0), 
     tck = 0.05)
lines(seq(0,1,0.2),  adj.power, col = "dodgerblue", 
      lty = 2, lwd = 4, yaxis = NULL, type = "l")
lines(seq(0,1,0.2),  lf.power, col = "midnightblue", 
      lty = 3, lwd = 4, yaxis = NULL,  type = "l")
lines(seq(0,1,0.2),  fh.power, col = "darkgreen", 
      lty = 4, lwd = 4, yaxis = NULL,  type = "l")
legend("topright", inset=c(-0.4, 0.5), 
       c(expression(MGC %.% DM), expression(MGC %.% AM),
         expression(paste( MGC %.% LF)), "FH Test"), seg.len = 3,
       col = c("red", "dodgerblue", "midnightblue", "darkgreen"),
       lty = c(1,2,3,4), lwd = 4, bty = 'n',  xpd = NA, cex = 2)
dev.off()

