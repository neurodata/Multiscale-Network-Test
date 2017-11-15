source(optimalT.R)

########################################################
load("../Data/multi_circle50_karl2.RData")
mgc.pval = c(); mcorr.pval = c(); hhg.pval = c(); fh.pval = c()
for(i in 1:length(multi_circle50_karl2)){
  mgc.pval[i] = print.stat.optimal(multi_circle50_karl2[[i]][[1]], 4)$pvalue
  mcorr.pval[i] = print.stat.optimal(multi_circle50_karl2[[i]][[2]], 4)$pvalue
  hhg.pval[i] = print.stat.optimal(multi_circle50_karl2[[i]][[3]], 4)$pvalue
  fh.pval[i] = multi_circle50_karl2[[i]][[4]]
}
circle.mgc = mean(mgc.pval <= 0.05); circle.mcorr = mean(mcorr.pval <= 0.05)
circle.hhg = mean(hhg.pval <= 0.05); circle.fh = mean(fh.pval <= 0.05)

# power of each time point t
mgc.pval = matrix(0, nrow = 100, ncol = 11)
mcorr.pval = matrix(0, nrow = 100, ncol = 11)
hhg.pval = matrix(0, nrow = 100, ncol = 11)
for(i in 1:length(multi_circle50_karl2)){
  for(t in 1:11){
    mgc.pval[i,t] = mean(multi_circle50_karl2[[i]][[1]][[1]][t] <= multi_circle50_karl2[[i]][[1]][[2]][,t])
    mcorr.pval[i,t] = mean(multi_circle50_karl2[[i]][[2]][[1]][t] <= multi_circle50_karl2[[i]][[2]][[2]][,t])
    hhg.pval[i,t] = mean(multi_circle50_karl2[[i]][[3]][[1]][t] <= multi_circle50_karl2[[i]][[3]][[2]][,t])
  }
}

mgc.power = colMeans(mgc.pval <= 0.05)
mcorr.power = colMeans(mcorr.pval <= 0.05)
hhg.power = colMeans(hhg.pval <= 0.05)

# optimal t 
mgc.alt.t = c(); mgc.alt.no0.t = c(); mgc.alt.de0.t = c()
mcorr.alt.t = c(); mcorr.alt.no0.t = c(); mcorr.alt.de0.t = c()
hhg.alt.t = c(); hhg.alt.no0.t = c(); hhg.alt.de0.t = c()

for(i in 1:length(multi_circle50_karl2)){
  mgc.alt.t[i]  = print.stat.optimal(multi_circle50_karl2[[i]][[1]], 4)$alt.t
  mcorr.alt.t[i]  = print.stat.optimal(multi_circle50_karl2[[i]][[2]], 4)$alt.t
  hhg.alt.t[i]  = print.stat.optimal(multi_circle50_karl2[[i]][[3]], 4)$alt.t
}

#### Distribution of chosen t ########
alter.t = matrix(0, nrow = 11, ncol = 3)
alter.t[c(1:11) %in% names(table(mgc.alt.t)),1] = table(mgc.alt.t) / 100
alter.t[c(1:11) %in% names(table(mcorr.alt.t)),2] = table(mcorr.alt.t) / 100
alter.t[c(1:11) %in% names(table(hhg.alt.t)),3] = table(hhg.alt.t) / 100
colnames(alter.t) = c("MGC", "mCorr", "HHG")
rownames(alter.t) = c(0:10)

pdf("../Figure/circle_optimal.pdf", width = 15, height = 8)
par(mfrow = c(1,1), mar = c(3,3,3,10), tcl = 0.5,
    xpd = TRUE, cex.axis = 2, cex.main = 2.5, oma = c(3,5,3,5))
barplot(t(alter.t), beside = TRUE, 
        col = c("red", "dodgerblue", "lightsalmon4"),
        ylim = c(0, 1), xpd=TRUE, cex.lab = 2, mgp = c(4,1,0),
        xlab = "", 
        ylab = "", main = "Circle")
par(new=TRUE)
plot(seq(1, 22, 2), mgc.power,  xlab="", ylab="", ylim=c(0,1.0), 
     axes=FALSE, type="b", col="red", lty = 2, pch = 20, cex = 2)
## a little farther out (line=4) to make room for labels
lines(seq(1, 22, 2), mcorr.power, type = "b", lty = 2, pch = 20, cex = 2, col = "dodgerblue")
lines(seq(1, 22, 2), hhg.power, type = "b", lty = 2, pch = 20, cex = 2, col = "lightsalmon4")

lines(seq(1, 22, 2), rep(circle.mgc,11), type = "l", cex = 2, col = "red", lwd = 4)
lines(seq(1, 22, 2), rep(circle.mcorr,11), type = "l", cex = 2, col = "dodgerblue", lwd = 4)
lines(seq(1, 22, 2), rep(circle.hhg,11), type = "l", cex = 2, col = "lightsalmon4", lwd = 4)
mtext("Power", side=4, cex=2.5, line=3) 
axis(4, ylim=c(0, 1.0), at = seq(0, 1.0, 0.1))
mtext(paste("Proportion of selecting t as an optimal", sep = "        "), side= 2, cex= 2.5,  xpd=TRUE, 
      line = 3.5,
      adj = 0.5)
mtext("Markov Time t", side= 1, cex= 2.5,  xpd=TRUE, 
      line = 3.5)
legend("topright", inset=c(-0.2, 0) ,  c("MGC", "dCorr", "HHG"),
       col = c("red", "dodgerblue", "lightsalmon4"), pch = 20, cex = 2,
       bty ='n', xpd = TRUE)
dev.off()

########################################################
load("../Data/multi_mono50_karl2.RData")
mgc.pval = c(); mcorr.pval = c(); hhg.pval = c(); fh.pval = c()
for(i in 1:length(multi_mono50_karl2)){
  mgc.pval[i] = print.stat.optimal(multi_mono50_karl2[[i]][[1]], 4)$pvalue
  mcorr.pval[i] = print.stat.optimal(multi_mono50_karl2[[i]][[2]], 4)$pvalue
  hhg.pval[i] = print.stat.optimal(multi_mono50_karl2[[i]][[3]], 4)$pvalue
  fh.pval[i] = multi_mono50_karl2[[i]][[4]]
}
circle.mgc = mean(mgc.pval <= 0.05); circle.mcorr = mean(mcorr.pval <= 0.05)
circle.hhg = mean(hhg.pval <= 0.05); circle.fh = mean(fh.pval <= 0.05)

# power of each time point t
mgc.pval = matrix(0, nrow = 100, ncol = 11)
mcorr.pval = matrix(0, nrow = 100, ncol = 11)
hhg.pval = matrix(0, nrow = 100, ncol = 11)
for(i in 1:length(multi_mono50_karl2)){
  for(t in 1:11){
    mgc.pval[i,t] = mean(multi_mono50_karl2[[i]][[1]][[1]][t] <= multi_mono50_karl2[[i]][[1]][[2]][,t])
    mcorr.pval[i,t] = mean(multi_mono50_karl2[[i]][[2]][[1]][t] <= multi_mono50_karl2[[i]][[2]][[2]][,t])
    hhg.pval[i,t] = mean(multi_mono50_karl2[[i]][[3]][[1]][t] <= multi_mono50_karl2[[i]][[3]][[2]][,t])
  }
}

mgc.power = colMeans(mgc.pval <= 0.05)
mcorr.power = colMeans(mcorr.pval <= 0.05)
hhg.power = colMeans(hhg.pval <= 0.05)

# optimal t 
mgc.alt.t = c(); mcorr.alt.t = c(); hhg.alt.t = c();

for(i in 1:length(multi_mono50_karl2)){
  mgc.alt.t[i]  = print.stat.optimal(multi_mono50_karl2[[i]][[1]], 4)$alt.t
  mcorr.alt.t[i]  = print.stat.optimal(multi_mono50_karl2[[i]][[2]], 4)$alt.t
  hhg.alt.t[i]  = print.stat.optimal(multi_mono50_karl2[[i]][[3]], 4)$alt.t
}

#### Distribution of chosen t ########
alter.t = matrix(0, nrow = 11, ncol = 3)
alter.t[c(1:11) %in% names(table(mgc.alt.t)),1] = table(mgc.alt.t) / 100
alter.t[c(1:11) %in% names(table(mcorr.alt.t)),2] = table(mcorr.alt.t) / 100
alter.t[c(1:11) %in% names(table(hhg.alt.t)),3] = table(hhg.alt.t) / 100
colnames(alter.t) = c("MGC", "mCorr", "HHG")
rownames(alter.t) = c(0:10)

pdf("../Figure/monoton50_optimal.pdf", width = 15, height = 8)
par(mfrow = c(1,1), mar = c(3,3,3,10), tcl = 0.5,
    xpd = TRUE, cex.axis = 2, cex.main = 2.5, oma = c(3,5,3,5))
barplot(t(alter.t), beside = TRUE, 
        col = c("red", "dodgerblue", "lightsalmon4"),
        ylim = c(0, 1), xpd=TRUE, cex.lab = 2, mgp = c(4,1,0),
        xlab = "", 
        ylab = "", main = expression(paste("SBM (", theta, "=0.50)", sep ="")))
par(new=TRUE)
plot(seq(1, 44, 4), mgc.power,  xlab="", ylab="", ylim=c(0,1.0), 
     axes=FALSE, type="b", col="red", lty = 2, pch = 20, cex = 2)
## a little farther out (line=4) to make room for labels
lines(seq(1, 44, 4), mcorr.power, type = "b", lty = 2, pch = 20, cex = 2, col = "dodgerblue")
lines(seq(1, 44, 4), hhg.power, type = "b", lty = 2, pch = 20, cex = 2, col = "lightsalmon4")

lines(seq(1, 44, 4), rep(circle.mgc,11), type = "l", cex = 2, col = "red", lwd = 4)
lines(seq(1, 44, 4), rep(circle.mcorr,11), type = "l", cex = 2, col = "dodgerblue", lwd = 4)
lines(seq(1, 44, 4), rep(circle.hhg,11), type = "l", cex = 2, col = "lightsalmon4", lwd = 4)
mtext("Power", side=4, cex=2.5, line=3) 
axis(4, ylim=c(0, 1.0), at = seq(0, 1.0, 0.1))
mtext(paste("Proportion of selecting t as an optimal", sep = "        "), side= 2, cex= 2.5,  xpd=TRUE, 
      line = 3.5,
      adj = 0.5)
mtext("Markov Time t", side= 1, cex= 2.5,  xpd=TRUE, 
      line = 3.5)
legend("topright", inset=c(-0.2, 0) ,  c("MGC", "dCorr", "HHG"),
       col = c("red", "dodgerblue", "lightsalmon4"), pch = 20, cex = 2,
       bty ='n', xpd = TRUE)
dev.off()

