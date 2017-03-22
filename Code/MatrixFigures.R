### R code to generate heatmap figures ###
library(ColorPalette)
library(RColorBrewer)
library(gplots)
library(DTMCPack)

load("../Data/SampleGraph.RData")
source("diffmaps.R")

GnPalette = colorRampPalette(brewer.pal(9,"Greens"))(20)

# G : sample graph
### align the row and column by block membership
ind1 = which(V(G)$group == 1) ; n1 = length(ind1)
ind2 = which(V(G)$group == 2) ; n2 = length(ind2)
ind3 = which(V(G)$group == 3) ; n3 = length(ind3)

### get an adjacency matrix
A = get.adjacency(G)
A = as.matrix(A)
A.mat = A[c(ind1, ind2, ind3), c(ind1, ind2, ind3)]

### get a transition matrix 
P = A.mat / pmax(rowSums(A.mat),1)
########### figures of probability matrix/Adjacency matrix/Transition matrix ######
pdf("../Figure/Amat.pdf")
par(cex.main=3, cex.lab = 2.5, cex.axis = 3, 
   mar=c(7,2,5,2), tcl = 0.5, oma = c(3, 5, 0, 2), xpd = NA)
heatmap.2(A.mat, dendrogram='none', 
          Rowv=FALSE, Colv=FALSE,trace='none', 
          xlab = "", 
          ylab = "", 
          main = "A",
          labRow = ""
          , labCol = "",
          keysize = 1,
          key.par=list(mar=c(3,0,0,4)),
          key.xlab = "",
          key.ylab = "",
          margins = c(5,5),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(0,3),c(2,1),c(0,4)), lhei=c(0.5, 2.5, 0.5), lwid=c(0.5,4),
          key.title = "",
          cexCol=2.5, col = GnPalette,
          srtCol=3,
          key.xtickfun=function() {
            breaks <- parent.frame()$breaks
            mtext("0", 1, cex = 3, line = -0.6, adj = 0.1, outer = TRUE, xpd = NA)
            mtext("0.5", 1, cex = 3, line = -0.6, adj = 0.5, outer = TRUE, xpd = NA)
            mtext("1", 1, cex = 3, line = -0.6, adj = 0.95, outer = TRUE, xpd = NA)
            return(list(
              at=parent.frame()$scale01(c(0, 0.5, 1)),
              labels=FALSE
            ))
          },
          breaks=seq(0, 1, 1/20),
          offsetRow = -45, density.info = "none")
#mtext("A", 3, cex = 4, line = -3.2, adj = 0.4, outer = TRUE, xpd = NA)
dev.off()


### get a distance matrix
X = V(G)$outcome[c(ind1, ind2, ind3)]
Dy = as.matrix(dist((X)), diag = TRUE, upper = TRUE)

U.list = dmap.q(P, c(1:50), q = 10)
## diffusion maps at t = 1
Dx.1 = as.matrix(dist((U.list[[1]])), diag = TRUE, upper = TRUE) 
Dx.2 = as.matrix(dist((U.list[[2]])), diag = TRUE, upper = TRUE) 
Dx.3 = as.matrix(dist((U.list[[3]])), diag = TRUE, upper = TRUE) 
Dx.5 = as.matrix(dist((U.list[[5]])), diag = TRUE, upper = TRUE)
Dx.10 = as.matrix(dist((U.list[[10]])), diag = TRUE, upper = TRUE)
Dx.20 = as.matrix(dist((U.list[[20]])), diag = TRUE, upper = TRUE)
Dx.50 = as.matrix(dist((U.list[[50]])), diag = TRUE, upper = TRUE)

pdf("../Figure/Q10T2.pdf")
par(cex.main=3, cex.lab = 2.5, cex.axis = 3, 
    mar=c(7,2,5,2), tcl = 0.5, oma = c(3, 5, 0, 2), xpd = NA)
heatmap.2(Dx.2, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
          xlab = "",
          ylab = "", 
          main = "t=2   q=10",
          labRow = ""
          , labCol = "",
          keysize = 1.5,
          key.par=list(mar=c(3,0,0,4)),
          key.xlab = "",
          key.ylab = "",
          margins = c(5,5),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(0,3),c(2,1),c(0,4)), lhei=c(0.5, 2.5, 0.5), lwid=c(0.5,4),
          key.title = "",
          cexRow=2.5,
          cexCol=2.5, col = GnPalette,
          srtCol=3,
          key.xtickfun=function() {
            breaks <- parent.frame()$breaks
            mtext("0", 1, cex = 3, line = -0.6, adj = 0.1, outer = TRUE, xpd = NA)
            mtext("0.36", 1, cex = 3, line = -0.6, adj = 0.5, outer = TRUE, xpd = NA)
            mtext("0.72", 1, cex = 3, line = -0.6, adj = 0.95, outer = TRUE, xpd = NA)
            return(list(
              at=parent.frame()$scale01(c(0, 0.5, 1)),
              labels=FALSE
            ))
          },
          breaks=seq(0, 0.72, 0.72/20),
          offsetRow = -45, density.info = "none")
dev.off()


###
pdf("../Figure/Q10T3.pdf")
par(cex.main=3, cex.lab = 2.5, cex.axis = 3, 
    mar=c(7,2,5,2), tcl = 0.5, oma = c(3, 5, 0, 2), xpd = NA)
heatmap.2(Dx.3, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
          xlab = "",
          ylab = "", 
          main = "t=3   q=10",
          labRow = ""
          , labCol = "",
          keysize = 1.5,
          key.par=list(mar=c(3,0,0,4)),
          key.xlab = "",
          key.ylab = "",
          margins = c(5,5),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(0,3),c(2,1),c(0,4)), lhei=c(0.5, 2.5, 0.5), lwid=c(0.5,4),
          key.title = "",
          cexRow=2.5,
          cexCol=2.5, col = GnPalette,
          srtCol=3,
          key.xtickfun=function() {
            breaks <- parent.frame()$breaks
            mtext("0", 1, cex = 3, line = -0.6, adj = 0.1, outer = TRUE, xpd = NA)
            mtext("0.13", 1, cex = 3, line = -0.6, adj = 0.5, outer = TRUE, xpd = NA)
            mtext("0.26", 1, cex = 3, line = -0.6, adj = 0.95, outer = TRUE, xpd = NA)
            return(list(
              at=parent.frame()$scale01(c(0, 0.5, 1)),
              labels=FALSE
            ))
          },
          breaks=seq(0, 0.26, 0.26/20),
          offsetRow = -45, density.info = "none")
dev.off()

###
pdf("../Figure/Q10T5.pdf")
par(cex.main=3, cex.lab = 2.5, cex.axis = 3, 
    mar=c(7,2,5,2), tcl = 0.5, oma = c(3, 5, 0, 2), xpd = NA)
heatmap.2(Dx.5, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
          xlab = "",
          ylab = "", 
          main = "t=5   q=10",
          labRow = ""
          , labCol = "",
          keysize = 1.5,
          key.par=list(mar=c(3,0,0,4)),
          key.xlab = "",
          key.ylab = "",
          margins = c(5,5),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(0,3),c(2,1),c(0,4)), lhei=c(0.5, 2.5, 0.5), lwid=c(0.5,4),
          key.title = "",
          cexRow=2.5,
          cexCol=2.5, col = GnPalette,
          srtCol=3,
          key.xtickfun=function() {
            breaks <- parent.frame()$breaks
            mtext("0", 1, cex = 3, line = -0.6, adj = 0.1, outer = TRUE, xpd = NA)
            mtext("0.02", 1, cex = 3, line = -0.6, adj = 0.5, outer = TRUE, xpd = NA)
            mtext("0.04", 1, cex = 3, line = -0.6, adj = 0.95, outer = TRUE, xpd = NA)
            return(list(
              at=parent.frame()$scale01(c(0, 0.5, 1)),
              labels=FALSE
            ))
          },
          breaks=seq(0, 0.04, 0.04/20),
          offsetRow = -45, density.info = "none")
dev.off()


###
pdf("../Figure/Q10T10.pdf")
par(cex.main=3, cex.lab = 2.5, cex.axis = 3, 
    mar=c(7,2,5,2), tcl = 0.5, oma = c(3, 5, 0, 2), xpd = NA)
heatmap.2(Dx.10, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
          xlab = "",
          ylab = "", 
          main = "t=10   q=10",
          labRow = ""
          , labCol = "",
          keysize = 1.5,
          key.par=list(mar=c(3,0,0,4)),
          key.xlab = "",
          key.ylab = "",
          margins = c(5,5),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(0,3),c(2,1),c(0,4)), lhei=c(0.5, 2.5, 0.5), lwid=c(0.5,4),
          key.title = "",
          cexRow=2.5,
          cexCol=2.5, col = GnPalette,
          srtCol=3,
          key.xtickfun=function() {
            breaks <- parent.frame()$breaks
            mtext("0", 1, cex = 3, line = -0.4, adj = 0.1, outer = TRUE, xpd = NA)
            mtext(expression(2 %.% 10^-3), 1, cex = 3, line = -0.4, adj = 0.5, outer = TRUE, xpd = NA)
            mtext(expression(4 %.% 10^-3), 1, cex = 3, line = -0.4, adj = 0.95, outer = TRUE, xpd = NA)
            return(list(
              at=parent.frame()$scale01(c(0, 0.5, 1)),
              labels=FALSE
            ))
          },
          breaks=seq(0, 4*10^(-4), 4*10^(-4)/20),
          offsetRow = -45, density.info = "none")
dev.off()


###################################################

Q1 = dmap.q(P, 3, q = 1)[[1]]
Q2 = dmap.q(P, 3, q = 2)[[1]]
Q10 = dmap.q(P, 3, q = 10)[[1]]
Q100 = dmap.q(P, 3, q = 100)[[1]] 
## diffusion maps at t = 1
Dx.q1 = as.matrix(dist(Q1), diag = TRUE, upper = TRUE) 
Dx.q2 = as.matrix(dist(Q2), diag = TRUE, upper = TRUE) 
Dx.q10 = as.matrix(dist(Q10), diag = TRUE, upper = TRUE) 
Dx.q100 = as.matrix(dist(Q100), diag = TRUE, upper = TRUE) 

###
pdf("../Figure/Q1T3.pdf")
par(cex.main=3, cex.lab = 2.5, cex.axis = 3, 
    mar=c(7,2,5,2), tcl = 0.5, oma = c(3, 5, 0, 2), xpd = NA)
heatmap.2(Dx.q1, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
          xlab = "",
          ylab = "", 
          main = "t=3   q=1",
          labRow = ""
          , labCol = "",
          keysize = 1.5,
          key.par=list(mar=c(3,0,0,4)),
          key.xlab = "",
          key.ylab = "",
          margins = c(5,5),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(0,3),c(2,1),c(0,4)), lhei=c(0.5, 2.5, 0.5), lwid=c(0.5,4),
          key.title = "",
          cexRow=2.5,
          cexCol=2.5, col = GnPalette,
          srtCol=3,
          key.xtickfun=function() {
            breaks <- parent.frame()$breaks
            mtext("0", 1, cex = 3, line = -0.1, adj = 0.1, outer = TRUE, xpd = NA)
            mtext(expression(5^-14), 1, cex = 3, line = -0.1, adj = 0.55, outer = TRUE, xpd = NA)
            mtext(expression(10^-14), 1, cex = 3, line = -0.1, adj = 1.0, outer = TRUE, xpd = NA)
            return(list(
              at=parent.frame()$scale01(c(0, 0.5, 1)),
              labels=FALSE
            ))
          },
          breaks=seq(0, 10^(-14), 10^(-14)/20),
          offsetRow = -45, density.info = "none")
dev.off()

###
pdf("../Figure/Q2T3.pdf")
par(cex.main=3, cex.lab = 2.5, cex.axis = 3, 
    mar=c(7,2,5,2), tcl = 0.5, oma = c(3, 5, 0, 2), xpd = NA)
heatmap.2(Dx.q2, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
          xlab = "",
          ylab = "", 
          main = "t=3   q=2",
          labRow = ""
          , labCol = "",
          keysize = 1.5,
          key.par=list(mar=c(3,0,0,4)),
          key.xlab = "",
          key.ylab = "",
          margins = c(5,5),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(0,3),c(2,1),c(0,4)), lhei=c(0.5, 2.5, 0.5), lwid=c(0.5,4),
          key.title = "",
          cexRow=2.5,
          cexCol=2.5, col = GnPalette,
          srtCol=3,
          key.xtickfun=function() {
            breaks <- parent.frame()$breaks
            mtext("0", 1, cex = 3, line = -0.6, adj = 0.1, outer = TRUE, xpd = NA)
            mtext("0.14", 1, cex = 3, line = -0.6, adj = 0.5, outer = TRUE, xpd = NA)
            mtext("0.27", 1, cex = 3, line = -0.6, adj = 0.95, outer = TRUE, xpd = NA)
            return(list(
              at=parent.frame()$scale01(c(0, 0.5, 1)),
              labels=FALSE
            ))
          },
          breaks=seq(0, 0.27, 0.27/20),
          offsetRow = -45, density.info = "none")
dev.off()

###
pdf("../Figure/Q10T3.pdf")
par(cex.main=3, cex.lab = 2.5, cex.axis = 3, 
    mar=c(7,2,5,2), tcl = 0.5, oma = c(3, 5, 0, 2), xpd = NA)
heatmap.2(Dx.q10, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
          xlab = "",
          ylab = "", 
          main = "t=3   q=10",
          labRow = ""
          , labCol = "",
          keysize = 1.5,
          key.par=list(mar=c(3,0,0,4)),
          key.xlab = "",
          key.ylab = "",
          margins = c(5,5),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(0,3),c(2,1),c(0,4)), lhei=c(0.5, 2.5, 0.5), lwid=c(0.5,4),
          key.title = "",
          cexRow=2.5,
          cexCol=2.5, col = GnPalette,
          srtCol=3,
          key.xtickfun=function() {
            breaks <- parent.frame()$breaks
            mtext("0", 1, cex = 3, line = -0.6, adj = 0.1, outer = TRUE, xpd = NA)
            mtext("0.14", 1, cex = 3, line = -0.6, adj = 0.5, outer = TRUE, xpd = NA)
            mtext("0.27", 1, cex = 3, line = -0.6, adj = 0.95, outer = TRUE, xpd = NA)
            return(list(
              at=parent.frame()$scale01(c(0, 0.5, 1)),
              labels=FALSE
            ))
          },
          breaks=seq(0, 0.27, 0.27/20),
          offsetRow = -45, density.info = "none")
dev.off()


###
pdf("../Figure/Q100T3.pdf")
par(cex.main=3, cex.lab = 2.5, cex.axis = 3, 
    mar=c(7,2,5,2), tcl = 0.5, oma = c(3, 5, 0, 2), xpd = NA)
heatmap.2(Dx.q100, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
          xlab = "",
          ylab = "", 
          main = "t=3   q=100",
          labRow = ""
          , labCol = "",
          keysize = 1.5,
          key.par=list(mar=c(3,0,0,4)),
          key.xlab = "",
          key.ylab = "",
          margins = c(5,5),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(0,3),c(2,1),c(0,4)), lhei=c(0.5, 2.5, 0.5), lwid=c(0.5,4),
          key.title = "",
          cexRow=2.5,
          cexCol=2.5, col = GnPalette,
          srtCol=3,
          key.xtickfun=function() {
            breaks <- parent.frame()$breaks
            mtext("0", 1, cex = 3, line = -0.6, adj = 0.1, outer = TRUE, xpd = NA)
            mtext("0.14", 1, cex = 3, line = -0.6, adj = 0.5, outer = TRUE, xpd = NA)
            mtext("0.27", 1, cex = 3, line = -0.6, adj = 0.95, outer = TRUE, xpd = NA)
            return(list(
              at=parent.frame()$scale01(c(0, 0.5, 1)),
              labels=FALSE
            ))
          },
          breaks=seq(0, 0.27, 0.27/20),
          offsetRow = -45, density.info = "none")
dev.off()


################### Adjacency Spectral Embedding ######################
A = as.matrix(get.adjacency(G))
elbows = getElbows(svd(A)$d, n =5)

E1 = embed_adjacency_matrix(G, no = 1, weights = NULL)[[1]][c(ind1, ind2, ind3),] 
E3 = embed_adjacency_matrix(G, no = 3, weights = NULL)[[1]][c(ind1, ind2, ind3),]
E10 = embed_adjacency_matrix(G, no = 10, weights = NULL)[[1]][c(ind1, ind2, ind3),]
E99 = embed_adjacency_matrix(G, no = 99, weights = NULL)[[1]][c(ind1, ind2, ind3),]

E.1 = as.matrix(dist(E1), diag = TRUE, upper = TRUE)
E.2 = as.matrix(dist(E2), diag = TRUE, upper = TRUE) 
E.3 = as.matrix(dist(E3), diag = TRUE, upper = TRUE)
E.10 = as.matrix(dist(E10), diag = TRUE, upper = TRUE) 
E.99 = as.matrix(dist(E99), diag = TRUE, upper = TRUE) 
                                                               
pdf("../Figure/E1.pdf")
par(cex.main=3, cex.lab = 2.5, cex.axis = 3, 
    mar=c(7,2,5,2), tcl = 0.5, oma = c(3, 5, 0, 2), xpd = NA)
heatmap.2(E.1, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
          xlab = "",
          ylab = "", 
          main = "q = 1",
          labRow = ""
          , labCol = "",
          keysize = 1.5,
          key.par=list(mar=c(3,0,0,4)),
          key.xlab = "",
          key.ylab = "",
          margins = c(5,5),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(0,3),c(2,1),c(0,4)), lhei=c(0.5, 2.5, 0.5), lwid=c(0.5,4),
          key.title = "",
          cexRow=2.5,
          cexCol=2.5, col = GnPalette,
          srtCol=3,
          key.xtickfun=function() {
            breaks <- parent.frame()$breaks
            mtext("0", 1, cex = 3, line = -0.6, adj = 0.1, outer = TRUE, xpd = NA)
            mtext("0.21", 1, cex = 3, line = -0.6, adj = 0.5, outer = TRUE, xpd = NA)
            mtext("0.42", 1, cex = 3, line = -0.6, adj = 0.95, outer = TRUE, xpd = NA)
            return(list(
              at=parent.frame()$scale01(c(0, 0.5, 1)),
              labels=FALSE
            ))
          },
          breaks=seq(0, 0.42, 0.42/20),
          offsetRow = -45, density.info = "none")
dev.off()

###
pdf("../Figure/E3.pdf")
par(cex.main=3, cex.lab = 2.5, cex.axis = 3, 
    mar=c(7,2,5,2), tcl = 0.5, oma = c(3, 5, 0, 2), xpd = NA)
heatmap.2(E.3, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
          xlab = "",
          ylab = "", 
          main = "q = 3",
          labRow = ""
          , labCol = "",
          keysize = 1.5,
          key.par=list(mar=c(3,0,0,4)),
          key.xlab = "",
          key.ylab = "",
          margins = c(5,5),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(0,3),c(2,1),c(0,4)), lhei=c(0.5, 2.5, 0.5), lwid=c(0.5,4),
          key.title = "",
          cexRow=2.5,
          cexCol=2.5, col = GnPalette,
          srtCol=3,
          key.xtickfun=function() {
            breaks <- parent.frame()$breaks
            mtext("0", 1, cex = 3, line = -0.6, adj = 0.1, outer = TRUE, xpd = NA)
            mtext("0.72", 1, cex = 3, line = -0.6, adj = 0.5, outer = TRUE, xpd = NA)
            mtext("1.44", 1, cex = 3, line = -0.6, adj = 0.95, outer = TRUE, xpd = NA)
            return(list(
              at=parent.frame()$scale01(c(0, 0.5, 1)),
              labels=FALSE
            ))
          },
          breaks=seq(0, 1.44, 1.44/20),
          offsetRow = -45, density.info = "none")
dev.off()


###
pdf("../Figure/E10.pdf")
par(cex.main=3, cex.lab = 2.5, cex.axis = 3, 
    mar=c(7,2,5,2), tcl = 0.5, oma = c(3, 5, 0, 2), xpd = NA)
heatmap.2(E.10, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
          xlab = "",
          ylab = "", 
          main = "q = 10",
          labRow = ""
          , labCol = "",
          keysize = 1.5,
          key.par=list(mar=c(3,0,0,4)),
          key.xlab = "",
          key.ylab = "",
          margins = c(5,5),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(0,3),c(2,1),c(0,4)), lhei=c(0.5, 2.5, 0.5), lwid=c(0.5,4),
          key.title = "",
          cexRow=2.5,
          cexCol=2.5, col = GnPalette,
          srtCol=3,
          key.xtickfun=function() {
            breaks <- parent.frame()$breaks
            mtext("0", 1, cex = 3, line = -0.6, adj = 0.1, outer = TRUE, xpd = NA)
            mtext("1.05", 1, cex = 3, line = -0.6, adj = 0.5, outer = TRUE, xpd = NA)
            mtext("2.10", 1, cex = 3, line = -0.6, adj = 0.95, outer = TRUE, xpd = NA)
            return(list(
              at=parent.frame()$scale01(c(0, 0.5, 1)),
              labels=FALSE
            ))
          },
          breaks=seq(0, 2.10, 2.10/20),
          offsetRow = -45, density.info = "none")
dev.off()

###
pdf("../Figure/E99.pdf")
par(cex.main=3, cex.lab = 2.5, cex.axis = 3, 
    mar=c(7,2,5,2), tcl = 0.5, oma = c(3, 5, 0, 2), xpd = NA)
heatmap.2(E.99, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
          xlab = "",
          ylab = "", 
          main = "q = 99",
          labRow = ""
          , labCol = "",
          keysize = 1.5,
          key.par=list(mar=c(3,0,0,4)),
          key.xlab = "",
          key.ylab = "",
          margins = c(5,5),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(0,3),c(2,1),c(0,4)), lhei=c(0.5, 2.5, 0.5), lwid=c(0.5,4),
          key.title = "",
          cexRow=2.5,
          cexCol=2.5, col = GnPalette,
          srtCol=3,
          key.xtickfun=function() {
            breaks <- parent.frame()$breaks
            mtext("0", 1, cex = 3, line = -0.6, adj = 0.1, outer = TRUE, xpd = NA)
            mtext("1.56", 1, cex = 3, line = -0.6, adj = 0.5, outer = TRUE, xpd = NA)
            mtext("3.12", 1, cex = 3, line = -0.6, adj = 0.95, outer = TRUE, xpd = NA)
            return(list(
              at=parent.frame()$scale01(c(0, 0.5, 1)),
              labels=FALSE
            ))
          },
          breaks=seq(0, 3.12, 3.12/20),
          offsetRow = -45, density.info = "none")
dev.off()
