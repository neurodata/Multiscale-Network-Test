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

popn = 100; M = 100
p = 0.5; q = 0.2
r = seq(0.05, 0.6, 0.05)

n.perm = 500; n.iter = n.perm

GeneralSBM = list(); results = list()

for(N in 1:12){
  # different block probabilities
  for(i in 1:M){
    set.seed(i)

    G = notsimple3block(popn, p, q, r)
    X = V(G)$outcome

    mgc.result =  Networktest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
    dcov.result =  Networktest.diffusion.stat(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
    hhg.result =  Networktest.diffusion.stat(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
        
   
    # estimate the matrix M based on SVD
    A = as.matrix(get.adjacency(G))
    SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
    # FH
    fh.k = min( max(getElbows(SVD.A$d, n = 2, plot = FALSE)), popn-1)
    fh.result = FH_test(A, X, k.range = fh.k, n.iter = n.iter)

    results[[i]] = list(mgc.result, dcov.result, hhg.result, fh.result))      
  }

  GeneralSBM[[N]] = results

}

save(GeneralSBM, file = "../Data/GeneralSBM.RData")


#########################################################
multi_mono05_karl2 = GeneralSBM[[1]]; multi_mono10_karl2 = GeneralSBM[[2]];  
multi_mono15_karl2 = GeneralSBM[[3]]; multi_mono20_karl2 = GeneralSBM[[4]]; 
multi_mono25_karl2 = GeneralSBM[[5]]; multi_mono30_karl2 = GeneralSBM[[6]]; 
multi_mono35_karl2 = GeneralSBM[[7]]; multi_mono40_karl2 = GeneralSBM[[8]];  
multi_mono45_karl2 = GeneralSBM[[9]]; multi_mono50_karl2 = GeneralSBM[[10]]; 
multi_mono55_karl2 = GeneralSBM[[11]]; multi_mono60_karl2 = GeneralSBM[[12]]

mono.mgc = matrix(0, nrow = 100, ncol = 12)
mono.mcorr = matrix(0, nrow = 100, ncol = 12)
mono.hhg = matrix(0, nrow = 100, ncol = 12)
mono.fh = matrix(0, nrow = 100, ncol = 12)
for(i in 1:100){
  mono.mgc[i,1] = print.stat.optimal(multi_mono05_karl2[[i]][[1]], 4)$pvalue; mono.mgc[i,2] = print.stat.optimal(multi_mono10_karl2[[i]][[1]], 4)$pvalue
  mono.mgc[i,3] = print.stat.optimal(multi_mono15_karl2[[i]][[1]], 4)$pvalue; mono.mgc[i,4] = print.stat.optimal(multi_mono20_karl2[[i]][[1]], 4)$pvalue
  mono.mgc[i,5] = print.stat.optimal(multi_mono25_karl2[[i]][[1]], 4)$pvalue; mono.mgc[i,6] = print.stat.optimal(multi_mono30_karl2[[i]][[1]], 4)$pvalue
  mono.mgc[i,7] = print.stat.optimal(multi_mono35_karl2[[i]][[1]], 4)$pvalue; mono.mgc[i,8] = print.stat.optimal(multi_mono40_karl2[[i]][[1]], 4)$pvalue
  mono.mgc[i,9] = print.stat.optimal(multi_mono45_karl2[[i]][[1]], 4)$pvalue; mono.mgc[i,10] = print.stat.optimal(multi_mono50_karl2[[i]][[1]], 4)$pvalue
  mono.mgc[i,11] = print.stat.optimal(multi_mono55_karl2[[i]][[1]], 4)$pvalue; mono.mgc[i,12] = print.stat.optimal(multi_mono60_karl2[[i]][[1]], 4)$pvalue
  
  mono.mcorr[i,1] = print.stat.optimal(multi_mono05_karl2[[i]][[2]], 4)$pvalue; mono.mcorr[i,2] = print.stat.optimal(multi_mono10_karl2[[i]][[2]], 4)$pvalue
  mono.mcorr[i,3] = print.stat.optimal(multi_mono15_karl2[[i]][[2]], 4)$pvalue; mono.mcorr[i,4] = print.stat.optimal(multi_mono20_karl2[[i]][[2]], 4)$pvalue
  mono.mcorr[i,5] = print.stat.optimal(multi_mono25_karl2[[i]][[2]], 4)$pvalue; mono.mcorr[i,6] = print.stat.optimal(multi_mono30_karl2[[i]][[2]], 4)$pvalue
  mono.mcorr[i,7] = print.stat.optimal(multi_mono35_karl2[[i]][[2]], 4)$pvalue; mono.mcorr[i,8] = print.stat.optimal(multi_mono40_karl2[[i]][[2]], 4)$pvalue
  mono.mcorr[i,9] = print.stat.optimal(multi_mono45_karl2[[i]][[2]], 4)$pvalue; mono.mcorr[i,10] = print.stat.optimal(multi_mono50_karl2[[i]][[2]], 4)$pvalue
  mono.mcorr[i,11] = print.stat.optimal(multi_mono55_karl2[[i]][[2]], 4)$pvalue; mono.mcorr[i,12] = print.stat.optimal(multi_mono60_karl2[[i]][[2]], 4)$pvalue
  
  mono.hhg[i,1] = print.stat.optimal(multi_mono05_karl2[[i]][[3]], 4)$pvalue; mono.hhg[i,2] = print.stat.optimal(multi_mono10_karl2[[i]][[3]], 4)$pvalue
  mono.hhg[i,3] = print.stat.optimal(multi_mono15_karl2[[i]][[3]], 4)$pvalue; mono.hhg[i,4] = print.stat.optimal(multi_mono20_karl2[[i]][[3]], 4)$pvalue
  mono.hhg[i,5] = print.stat.optimal(multi_mono25_karl2[[i]][[3]], 4)$pvalue; mono.hhg[i,6] = print.stat.optimal(multi_mono30_karl2[[i]][[3]], 4)$pvalue
  mono.hhg[i,7] = print.stat.optimal(multi_mono35_karl2[[i]][[3]], 4)$pvalue; mono.hhg[i,8] = print.stat.optimal(multi_mono40_karl2[[i]][[3]], 4)$pvalue
  mono.hhg[i,9] = print.stat.optimal(multi_mono45_karl2[[i]][[3]], 4)$pvalue; mono.hhg[i,10] = print.stat.optimal(multi_mono50_karl2[[i]][[3]], 4)$pvalue
  mono.hhg[i,11] = print.stat.optimal(multi_mono55_karl2[[i]][[3]], 4)$pvalue; mono.hhg[i,12] = print.stat.optimal(multi_mono60_karl2[[i]][[3]], 4)$pvalue
  
  mono.fh[i,1] = multi_mono05_karl2[[i]][[4]]; mono.fh[i,2] = multi_mono10_karl2[[i]][[4]]
  mono.fh[i,3] = multi_mono15_karl2[[i]][[4]]; mono.fh[i,4] = multi_mono20_karl2[[i]][[4]]
  mono.fh[i,5] = multi_mono25_karl2[[i]][[4]]; mono.fh[i,6] = multi_mono30_karl2[[i]][[4]]
  mono.fh[i,7] = multi_mono35_karl2[[i]][[4]]; mono.fh[i,8] = multi_mono40_karl2[[i]][[4]]
  mono.fh[i,9] = multi_mono45_karl2[[i]][[4]]; mono.fh[i,10] = multi_mono50_karl2[[i]][[4]]
  mono.fh[i,11] = multi_mono55_karl2[[i]][[4]]; mono.fh[i,12] = multi_mono60_karl2[[i]][[4]]
  
}

mgc.power = colMeans(mono.mgc <= 0.05); mcorr.power = colMeans(mono.mcorr <= 0.05)
hhg.power = colMeans(mono.hhg <= 0.05); fh.power = colMeans(mono.fh <= 0.05)

pdf("../Figure/multi_monoton.pdf", width = 15, height = 6)
par(mfrow = c(1,1), cex.lab = 3, cex.axis = 2,
    mar = c(6,8,3,17), tcl = 0.5)
plot(seq(0.05, 0.60, 0.05),  mgc.power, col = "red", 
      lty = 1, lwd = 5, ylab = "Power",
      ylim = c(0, 0.8), type = "l", mgp = c(5,2,0),
      xlab = expression(theta), yaxt = 'n')
axis(side = 2, at = seq(0.0, 0.8, 0.2), 
     labels = seq(0.0, 0.8, 0.2), 
     tck = 0.05)
lines( seq(0.05, 0.60, 0.05),  mcorr.power, col = "dodgerblue", 
       lty = 2, lwd = 5,  type = "l")
lines( seq(0.05, 0.60, 0.05),  hhg.power, col = "lightsalmon4", 
       lty = 5, lwd = 5,  type = "l")
lines( seq(0.05, 0.60, 0.05),  fh.power, col = "darkgreen", 
       lty =4, lwd = 5,  type = "l")
legend("topright", inset=c(-0.32, 0.5), 
       c( expression(MGC %.% DM), expression(dCorr %.% DM), 
          expression(HHG %.% DM), "FH Test"), 
       col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"), seg.len = 3,
       lty = c(1,2,5,4), lwd = 4, bty = 'n', cex = 2, xpd = NA)
abline(v = 0.2, lwd = 2, col = "black")
dev.off()