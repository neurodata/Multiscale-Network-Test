library(ColorPalette)
library(RColorBrewer)
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

RedPalette = colorRampPalette(brewer.pal(9,"Reds"))(20)

source("notsimple3block.R")
source("diffmaps.R")
source("dmaps.R")
source("elbowmap.R")
source("FH_factor.R")
source("FH_test.R")
source("getElbows.R")
source("NetworkTest.R")

mgc.DM.result = c(); mcorr.DM.result = c(); hhg.DM.result = c()
mgc.AM.result = c(); mcorr.AM.result = c(); hhg.AM.result = c()
mgc.LF.result = c(); mcorr.LF.result = c(); hhg.LF.result = c()
fh.result = c()

popn = 100 # size of node; sample size
n.group = 3 # number of blocks
M = 500 # number of iterations
 
p = 0.5; q = 0.2; r = 0.3

dstep = 3; n.perm = 500; n.iter = n.perm

for(i in 1:M){
  set.seed(i)
  G = notsimple3block(popn, p, q, r)
  X = V(G)$outcome
      
  A = as.matrix(get.adjacency(G))
      
  P = A  / pmax(rowSums(A),1)
  diffusion.q  =  min( max(getElbows(print.lambda(P, times = 3)[[1]], plot = FALSE, n = 3)), popn-1)

  mgc.DM.result[i] =  NetworkTest.q(G, X, option = 1, diffusion = TRUE, dstep = 3, n.perm = n.perm, q = diffusion.q)[[1]][[1]]
  mcorr.DM.result[i] =  NetworkTest.q(G, X, option = 2, diffusion = TRUE, dstep = 3, n.perm = n.perm, q = diffusion.q)
  hhg.DM.result[i] =  NetworkTest.q(G, X, option = 3, diffusion = TRUE, dstep = 3, n.perm = n.perm, q = diffusion.q)

  mgc.AM.result[i] =  NetworkTest.q(G, X, option = 1, diffusion = FALSE, dstep = 3, n.perm = n.perm, q = diffusion.q)[[1]][[1]]
  mcorr.AM.result[i] =  NetworkTest.q(G, X, option = 2, diffusion = FALSE, dstep = 3, n.perm = n.perm, q = diffusion.q)
  hhg.AM.result[i] =  NetworkTest.q(G, X, option = 3, diffusion = FALSE, dstep = 3, n.perm = n.perm, q = diffusion.q)

  SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
  fh.k = min( max(getElbows(SVD.A$d, n = 3, plot = FALSE)), popn-1)
  factors = FH_factor(A, k.range = fh.k) # latent position
  
  Dx = as.matrix(dist(factors[[1]]), diag = TRUE, upper = TRUE) # distance matrix of latent positions
  Dy = as.matrix(dist(X), diag = TRUE, upper = TRUE) # distance matrix of nodal atttributes
  
  mgc.LF.result[i] = MGCPermutationTest(Dx, Dy, rep = n.perm, option = 'mcor')[[1]]
  mcorr.LF.result[i] = dcor.ttest(Dx, Dy, distance = TRUE)$p.value
  hhg.LF.result[i] = hhg.test(Dx, Dy, nr.perm = n.perm)$perm.pval.hhg.sl

  fh.result[i] = FH_test(A, X, k.range = fh.k, n.iter = n.iter)
}


ThreeSBM = cbind(mgc.DM.result, mcorr.DM.result, hhg.DM.result,
                mgc.AM.result, mcorr.AM.result, hhg.AM.result,
                mgc.LF.result, mcorr.LF.result, hhg.LF.result, fh.result)
write.csv(ThreeSBM, file = "../Data/ThreeSBM.csv", row.names = FALSE)


M = 500 # number of iterations
alpha = 0.05 # type-I error

mgc.DM.power = 0; mcorr.DM.power = 0; hhg.DM.power = 0
mgc.AM.power = 0; mcorr.AM.power = 0; hhg.AM.power = 0
mgc.LF.power = 0; mcorr.LF.power = 0; hhg.LF.power = 0
fh.power = 0; 

for(i in 1:M){
  mgc.DM.power = mgc.DM.power + (mgc.DM.result[i] <= alpha) / M
  mcorr.DM.power = mcorr.DM.power + (mcorr.DM.result[i] <= alpha) / M
  hhg.DM.power = hhg.DM.power + (hhg.DM.result[i] <= alpha) / M
  
  mgc.AM.power = mgc.AM.power + (mgc.AM.result[i] <= alpha) / M
  mcorr.AM.power = mcorr.AM.power + (mcorr.AM.result[i] <= alpha) / M
  hhg.AM.power = hhg.AM.power + (hhg.AM.result[i] <= alpha) / M
  
  mgc.LF.power = mgc.LF.power + (mgc.LF.result[i] <= alpha) / M
  mcorr.LF.power = mcorr.LF.power + (mcorr.LF.result[i] <= alpha) / M
  hhg.LF.power = hhg.LF.power + (hhg.LF.result[i] <= alpha) / M
  
  fh.power = fh.power + (fh.result[i] <= alpha) / M
}


#################################################
mat = matrix(0, nrow = 3, ncol = 4)
mat[1,] <- c(mgc.DM.power, mcorr.DM.power, hhg.DM.power, 0)
mat[2,] <- c(mgc.AM.power, mcorr.AM.power, hhg.AM.power, 0)
mat[3,] <- c(mgc.LF.power, mcorr.LF.power, hhg.LF.power, fh.power)

png("figure/ThreeSBM_Elbow3.png")
par(cex.main=2, cex.lab = 2.5, cex.axis = 3, 
    mar=c(7,7,5,2), tcl = 0.5, oma = c(3, 5, 0, 2), xpd = NA)
heatmap.2(mat, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
          xlab = "",
          ylab = "", 
          main = "Empirical Power (n = 100)",
          labRow = c("DM","AM"," LF")
          , labCol = c("MGC", 
                       "mCorr", "HHG", "FH"), 
          keysize = 1.5,
          key.par=list(mar=c(3,0,0,4)),
          key.xlab = "",
          key.ylab = "",
          margins = c(5,5),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(0,3),c(2,1),c(0,4)), lhei=c(0.7, 2.5, 0.5), lwid=c(0.5,4),
          key.title = "", mgp = c(3,1,0),
          cexRow=2.5,
          cexCol=2.5, col = RedPalette,
          srtCol=0, adjCol = c(0.3,-20),
          adjRow = c(-2,0),
          key.xtickfun=function() {
            breaks <- parent.frame()$breaks
            mtext("0", 1, cex = 2, line = -1.1, adj = 0.1, outer = TRUE, xpd = NA)
            mtext("0.2", 1, cex = 2, line = -1.1, adj = 0.32, outer = TRUE, xpd = NA)
            mtext("0.4", 1, cex = 2, line = -1.1, adj = 0.62, outer = TRUE, xpd = NA)
            mtext("0.6", 1, cex = 2, line = -1.1, adj = 0.92, outer = TRUE, xpd = NA)
            return(list(
              at=parent.frame()$scale01(c(0, 0.2, 0.4, 0.6)),
              labels=FALSE
            ))
          },
          breaks=seq(0, 0.6, 0.6/20),
          offsetRow = -45, density.info = "none")
title(xlab = "Test statistics", line = 1.5, adj = 0)
title(ylab = "Metrics", line = 10, adj = 0.5)
dev.off()