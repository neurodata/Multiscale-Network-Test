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
source("FH_test.R")
source("getElbows.R")
source("NetworkTest.R")
source("optimalT.R")


popn = 100 # size of node; sample size
M = 100 # number of iterations
 
p = 0.5; q = 0.2; r = 0.4
n.perm = 500; n.iter = n.perm
adjthree = lfthree = dfthree = list()

for(i in 1:M){
  set.seed(i)
  G = notsimple3block(popn, p, q, r)
  X = V(G)$outcome
  A = as.matrix(get.adjacency(G))
      
  mgc.result =  NetworkTest.diffusion.stat.sym(G, X, option = 1, diffusion = FALSE, t.range = c(0:10), n.perm = n.perm)
  dcov.result =  NetworkTest.diffusion.stat.sym(G, X, option = 2, diffusion = FALSE, t.range = c(0:10), n.perm = n.perm)
  hhg.result =  NetworkTest.diffusion.stat.sym(G, X, option = 3, diffusion = FALSE, t.range = c(0:10), n.perm = n.perm)
  adjthree[[i]] = list(mgc.result, dcov.result, hhg.result)


  # estimate the matrix M based on SVD
  SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
  # FH
  fh.k = min( max(getElbows(SVD.A$d, n = 2, plot = FALSE, threshold = 0)), popn-1)
  fh.result = FH_test(A, X, k.range = fh.k, n.iter = n.iter)
  factors = FH_factor(A, k.range = fh.k)
  Dx = as.matrix(dist(factors[[1]]), diag = TRUE, upper = TRUE)
  Dy = as.matrix(dist(X), diag = TRUE, upper = TRUE)
  mgc.result = MGCPermutationTest(Dx, Dy, rep = n.perm, option = 'mcor')
  dcov.result = dcor.ttest(Dx, Dy, distance = TRUE)$p.value
  hhg.result = hhg.test(Dx, Dy, nr.perm = n.perm)$perm.pval.hhg.sl
  lfthree[[i]] = list(mgc.result, dcov.result, hhg.result, fh.result)  

  mgc.result =  NetworkTest.diffusion.stat.sym(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  dcov.result =  NetworkTest.diffusion.stat.sym(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  hhg.result =  NetworkTest.diffusion.stat.sym(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  dfthree[[i]] = list(mgc.result, dcov.result, hhg.result)             
}

ThreeSBM = list(adjthree, lfthree, dfthree)

save(ThreeSBM, file = "../Data/ThreeSBM.RData")

###
multi_Adjthree_karl2 = ThreeSBM[[1]]
multi_LFthree_karl2 = ThreeSBM[[2]]
multi_mono40_karl2 = ThreeSBM[[3]]
###
mgc.DM.result = c(); dcov.DM.result = c(); hhg.DM.result = c()
mgc.AM.result = c(); dcov.AM.result = c(); hhg.AM.result = c()
mgc.LF.result = c(); dcov.LF.result = c(); hhg.LF.result = c()
fh.result = c(); fh.result2 = c()


mgc.DM.power = 0; dcov.DM.power = 0; hhg.DM.power = 0
mgc.AM.power = 0; dcov.AM.power = 0; hhg.AM.power = 0
mgc.LF.power = 0; dcov.LF.power = 0; hhg.LF.power = 0
fh.power = 0;
###
for(i in 1:M){
  mgc.DM.result[i] = print.stat.optimal(multi_mono40_karl2[[i]][[1]], 4)$pvalue
  dcov.DM.result[i] = print.stat.optimal(multi_mono40_karl2[[i]][[2]], 4)$pvalue
  hhg.DM.result[i] = print.stat.optimal(multi_mono40_karl2[[i]][[3]], 4)$pvalue

  mgc.AM.result[i] = multi_Adjthree_karl2[[i]][[1]]$pMGC
  dcov.AM.result[i] = multi_Adjthree_karl2[[i]][[2]]
  hhg.AM.result[i] = multi_Adjthree_karl2[[i]][[3]]
  
  mgc.LF.result[i] = multi_LFthree_karl2[[i]][[1]]$pMGC
  dcov.LF.result[i] = multi_LFthree_karl2[[i]][[2]]
  hhg.LF.result[i] = multi_LFthree_karl2[[i]][[3]]
  fh.result[i] = multi_LFthree_karl2[[i]][[4]]
  
  mgc.DM.power = mgc.DM.power + (mgc.DM.result[i] <= alpha) / M
  dcov.DM.power = dcov.DM.power + (dcov.DM.result[i] <= alpha) / M
  hhg.DM.power = hhg.DM.power + (hhg.DM.result[i] <= alpha) / M
  
  mgc.AM.power = mgc.AM.power + (mgc.AM.result[i] <= alpha) / M
  dcov.AM.power = dcov.AM.power + (dcov.AM.result[i] <= alpha) / M
  hhg.AM.power = hhg.AM.power + (hhg.AM.result[i] <= alpha) / M
  
  mgc.LF.power = mgc.LF.power + (mgc.LF.result[i] <= alpha) / M
  dcov.LF.power = dcov.LF.power + (dcov.LF.result[i] <= alpha) / M
  hhg.LF.power = hhg.LF.power + (hhg.LF.result[i] <= alpha) / M
  
  fh.power = fh.power + (fh.result[i] <= alpha) / M
}

###
mat = matrix(0, nrow = 3, ncol = 4)
rownames(mat) = c("DM", "AM", "LF")
colnames(mat) = c("MGC", "mCorr", "HHG", "FH")
mat[1,] = c(mgc.DM.power, dcov.DM.power, hhg.DM.power,NA)
mat[2,] = c(mgc.AM.power, dcov.AM.power, hhg.AM.power,NA)
mat[3,] = c(mgc.LF.power, dcov.LF.power, hhg.LF.power, fh.power)
###
pdf("../Figure/multi_ThreeSBM.pdf")
par(cex.main=2, cex.lab = 2.5, cex.axis = 3, 
    mar=c(7,7,5,2), tcl = 0.5, oma = c(3, 5, 0, 2), xpd = NA)
heatmap.2(mat, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
          xlab = "",
          ylab = "", 
          main = "Empirical Power (n = 100)",
          labRow = c("DM","AM"," LF")
          , labCol = c("MGC", 
                       "dCorr", "HHG", "FH"), 
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
          adjRow = c(-1,0), cellnote= formatC(mat, digits=2, format = "f"),
          notecol = "black", notecex = 2,
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
          breaks=seq(0, 0.60, 0.60/20),
          offsetRow = -45, density.info = "none")
title(xlab = "Test statistics", line = 1.5, adj = 0)
title(ylab = "Metrics", line = 10, adj = 0.5)
dev.off()
