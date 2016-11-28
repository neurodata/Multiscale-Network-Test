library(igraph)
library(network)
library(MASS)
library(xtable)
library(parallel)
library(foreach)
library(HHG)
library(lattice)
library(DTMCPack)
library(Matrix)
library(energy)
library(Rlab)
library(VGAM)
library(amen)
library(doParallel)
library(ecodist)
library(SDMTools)

#source("R/DistCentering.R")
#source("R/DistRanks.R")
#source("R/FindLargestRectangles.R")
#source("R/MGCLocalCorr.R")
#source("R/MGCSampleStat.R")
#source("R/MGCPermutationTest.R")
#source("R/MGCSampleStat.R")
source("Code/dmaps.R")
source("Code/FH_test.R")

mortonXYZ <- function(pos) {
  x = c()
  y = c()
  z = c()
  idx <- DecToBin(pos)
  idx <- strsplit(idx, '')
  idx <- idx[[1]]
  idx <- idx[(length(idx)-24+1):length(idx)]
  length(idx)
  while (!is.na(idx[1])) {
    z <- c(z, idx[1])
    y <- c(y, idx[2])
    x <- c(x, idx[3])
    idx = idx[4:length(idx)]
  }
  triple <- c(BinToDec(x), BinToDec(y), BinToDec(z))
  return(triple)
}

DecToBin <- function(x) {
  paste(sapply(strsplit(paste(rev(intToBits(x))),""),`[[`,2),collapse="")
}

BinToDec <- function(x) {
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
}


## read Brain data
Brain <- read.graph("data/KKI2009_113_1_bg.graphml", format = "graphml")

# get only one largest connected component 
cl <- clusters(Brain)
brain1 <- induced_subgraph(Brain, which(cl[[1]] == which.max(cl[[2]] == 95)))


# return a triple of (X,Y,Z) location 
position1 <- matrix(0, nrow = length(V(brain1)), ncol = 3)
for(i in 1:length(V(brain1))){
  #print(i)
  position1[i,] <- mortonXYZ(V(brain1)$spatial_id[i])
}
dim(position1)

A <- as.matrix(get.adjacency(brain1)) # symmetric! yay!
X <- position1

dstep <- c(1:10)
n.perm <- 500
n.iter <- n.perm
k.range <- c(0:10)
  
A <- as.matrix(get.adjacency(brain1))
P <- A / rowSums(A)
brain1.dmaps <- dmap(P, c(1:10))

Dx <- list()
# Euclidean distance of nodal attributes
Dy = as.matrix(dist((X)), diag = TRUE, upper = TRUE) 

ldcorr <- c(); lmdcorr <- c(); hhgr <- c()

for(s in 1:length(dstep)){
  # Euclidean distance of diffusion maps
  Dx[[s]] = as.matrix(dist((brain1.dmaps[[s]])), diag = TRUE, upper = TRUE) 
  ldcorr[s] <- dcov.test(Dx[[s]], Dy, index = 1.0, R = n.perm)$p.value
  hhgr[s] <- hhg.test(Dx[[s]], Dy, nr.perm = 1000)$perm.pval.hhg.sl
}

fh.results <- FH_test(A, X, k.range = k.range, n.iter = n.iter)

######################################################
############ 3d graph ###############################

library(plot3D)
library(plotly)
x = position1[,1]
y = position1[,2]
z = position1[,3]
library(scatterplot3d)

pdf("Figure/intro.pdf")
par(mfrow = c(1,1), cex.lab = 2, cex.axis = 1.5, cex.main = 2,
    oma=c(1,1,1,1), mar = c(5,7,4,2) + 0.1, las = 2)
s <- scatterplot3d(x,y,z, main="", color = "skyblue",
                   tick.marks = FALSE, label.tick.marks = TRUE,
                   grid = TRUE, 
                   pch = 19, box = FALSE,
                   angle = 40,
                   scale.y = 1.3, xlab = "x",
                   ylab = "", zlab = "",
                   cex.symbols = 1,
                   cex.axis = 2,
                   cex.lab = 2, col.grid="blue",
                   col.axis = "blue",y.margin.add=0,
                   las = 1)
mtext( "y", side = 4, line = -4, outer = TRUE, cex = 2, srt = 0, padj = 1)
mtext( "z", side = 2, line = -2, outer = TRUE, cex = 2, srt = 90)
s$points3d(position1, col = gg.col(100), pch = 16, cex = 1.5)

## now draw a line between points 2 and 3
for(i in 1:length(V(brain1))){
  assign( paste("p", i, sep=""), s$xyz.convert(x[i],y[i],z[i]))
}

num.E = length(E(brain1))
edge.list = get.edgelist(brain1)
for(i in 1:num.E){
  tmp.p1 <- eval(parse(text = paste("p", edge.list[i,1], sep="")))
  tmp.p2 <- eval(parse(text = paste("p", edge.list[i,2], sep="")))
  segments(tmp.p1$x, tmp.p1$y, tmp.p2$x, tmp.p2$y, 
           lwd=1, col="grey80")
}
dev.off()


