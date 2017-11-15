source("diffmaps.R")
source("NetworkTest.R")
source("optimalT.R")

library(igraph)

### auxiliary functions to obtain physical coordinates of voxel.
mortonXYZ = function(pos) {
  x = c()
  y = c()
  z = c()
  idx = DecToBin(pos)
  idx = strsplit(idx, '')
  idx = idx[[1]]
  idx = idx[(length(idx)-24+1):length(idx)]
  length(idx)
  while (!is.na(idx[1])) {
    z = c(z, idx[1])
    y = c(y, idx[2])
    x = c(x, idx[3])
    idx = idx[4:length(idx)]
  }
  triple = c(BinToDec(x), BinToDec(y), BinToDec(z))
  return(triple)
}

DecToBin = function(x) {
  paste(sapply(strsplit(paste(rev(intToBits(x))),""),`[[`,2),collapse="")
}

BinToDec = function(x) {
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
}

#############################
# obtain sub-brain network
# "KKI2009_113_1_bg.graphml" from
# http://openconnecto.me/data/public/MR/m2g_v1_1_1/KKI2009/derivatives/bg/
Brain = read.graph("KKI2009_113_1_bg.graphml", format = "graphml")

# get only one largest connected component
cl <- clusters(Brain)
SubBrain <- induced_subgraph(Brain, which(cl[[1]] == which.max(cl[[2]] == 95)))

# load sub-brain network 
load("../Data/SubBrain.RData")

# return a triple of (X,Y,Z) location 
position = matrix(0, nrow = length(V(SubBrain)), ncol = 3)
for(i in 1:length(V(SubBrain))){
  position[i,] = mortonXYZ(V(SubBrain)$spatial_id[i])
}

A = as.matrix(get.adjacency(SubBrain)) # symmetric! yay!
P = A / pmax(rowSums(A),1)
X = position
Dy = as.matrix(dist(X), diag = TRUE, upper = TRUE)

## dimensional choice of network metrics
SVD.A = svd(A, nu = 95, nv = 95, LINPACK = FALSE)
fh.k0 =  min( max(getElbows(SVD.A$d, n = 3, plot = FALSE)), 95-1 )
diffusion.q0  = max(getElbows(print.lambda(P, times = 3)[[1]], plot = FALSE, n = 3, threshold = 0))


## network test without contamination
mgc.pval0 = NetworkTest.q(SubBrain, X, option = 1, diffusion = TRUE, dstep = 5, 
                            n.perm = 500, q = diffusion.q0)[[1]][[1]]

mcorr.pval0 = NetworkTest.q(SubBrain, X, option = 2, diffusion = TRUE, dstep = 5, 
                               n.perm = 500, q = diffusion.q0)

hhg.pval0 = NetworkTest.q(SubBrain, X, option = 3, diffusion = TRUE, dstep = 5, 
                              n.perm = 500, q = diffusion.q0)

fh.pval0 = FH_test(A, X, k.range = fh.k0, n.iter = 500)


############### contamination process ######################
mgc.results = list(); mcorr.results = list(); hhg.results = list()
fh.k = c(); fh.pval = c()
contam.G.adj = list(); multi_brainanalysis = list()

for(ii in 1:100){
  ## ii-th iteration
  for(nr in 1:10){
    # the number of deleted edges
    n.delete = as.integer(nr*0.1*length(E(brain1))) 
  
    # randomly sample deleted edges and deleted them from the graph
    set.seed(ii*nr)
  
    delete = sample(E(brain1), n.delete)
    delete.G = brain1 - delete
    delete.adj = as.matrix(get.adjacency(delete.G))
  
    # create the dummy adjacency matrix which just flips every edge
    dummy.adj = matrix(0, 95, 95)
    for(i in 1:95){
      for(j in i:95){
        dummy.adj[i,j] = ifelse(A[i,j] == 1, 0, 1)
        dummy.adj[j,i] = dummy.adj[i,j]
        dummy.adj[i,i] = 0
      }
    }
    dummy.G = graph.adjacency(dummy.adj, "undirected") 
  
    # the number of not-added edges from the flipped graph ((10*(10-nr))%)
    n.add = as.integer(0.1*(10-nr)*length(E(dummy.G))) 
  
    # randomly sampled not-added edges and deleted them from the flipped graph.
    set.seed(ii*nr)
    add = sample(E(dummy.G), n.add)
    dummy.G = dummy.G - add
    dummy.adj = as.matrix(get.adjacency(dummy.G))
  
    # sum of the (1) (10*nr) % deleted edges and
    # (2) remaining (10*nr) % edges from flipped graph. 
    contam.G.adj[[nr]] = dummy.adj + delete.adj

    G = graph.adjacency(contam.G.adj[[nr]], "undirected")
  
    # if G is not connected implement the above procedure again
    while(sum(colSums(contam.G.adj[[nr]]) ==0 ) > 0){
      for(i in 1:95){
        for(j in i:95){
          dummy.adj[i,j] = ifelse(A[i,j] == 1, 0, 1)
          dummy.adj[j,i] = dummy.adj[i,j]
          dummy.adj[i,i] = 0
        }
      }
    
      dummy.G = graph.adjacency(dummy.adj, "undirected") 
      n.add = as.integer(0.1*(10-nr)*length(E(dummy.G))) 
      add = sample(E(dummy.G), n.add)
    
      dummy.G = dummy.G - add
      dummy.adj = as.matrix(get.adjacency(dummy.G))
    
      contam.G.adj[[nr]] = dummy.adj + delete.adj
    
      G = graph.adjacency(contam.G.adj[[nr]], "undirected")
    }
  
    # contaminated transition matrix
    contam.P = contam.G.adj[[nr]] / rowSums(contam.G.adj[[nr]]) 
  
    # choose the dimension of diffusion map
    mgc.results[[nr]] =  NetworkTest.diffusion.stat.multi(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
    mcorr.results[[nr]] =  NetworkTest.diffusion.stat.multi(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
    hhg.results[[nr]] =  NetworkTest.diffusion.stat.multi(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
        
    # estimate the matrix M based on SVD
    SVD.A = svd(contam.G.adj[[nr]], nu = 95, nv = 95, LINPACK = FALSE)
  
    # FH
    fh.k[nr] =  min( max(getElbows(SVD.A$d, n = 2, plot = FALSE, threshold = 0)), 95-1)
    fh.pval[nr] =  FH_test(contam.G.adj[[nr]], X, k.range = fh.k[nr], n.iter = 500)
  } 

  multi_brainanalysis[[ii]] = list(mgc.results, mcorr.results, hhg.results, fh.pval)

}

############### save the result ############
save(multi_brainanalysis, file = "../Data/multi_brainanalysis.RData")

## arrange by contamination level 
contam.mgc = list(); contam.mcorr = list(); contam.hhg = list(); 
contam.fh = matrix(0, nrow = length(multi_brainanalysis), ncol = 10)
for(r in 1:10){
  contam.mgc[[r]] = list(); contam.mcorr[[r]] = list();
  contam.hhg[[r]] = list();
  for(i in 1:length(multi_brainanalysis)){
    contam.mgc[[r]][[i]] = multi_brainanalysis[[i]][[1]][[r]]
    contam.mcorr[[r]][[i]] = multi_brainanalysis[[i]][[2]][[r]]
    contam.hhg[[r]][[i]] = multi_brainanalysis[[i]][[3]][[r]]
  }
}
for(i in 1:length(multi_brainanalysis)){
  contam.fh[i,] = multi_brainanalysis[[i]][[4]]
}

pval.mgc = matrix(0, nrow = length(multi_brainanalysis), ncol = 10)
pval.mcorr = matrix(0, nrow = length(multi_brainanalysis), ncol = 10)
pval.hhg = matrix(0, nrow = length(multi_brainanalysis), ncol = 10)
pval.fh = matrix(0, nrow = length(multi_brainanalysis), ncol = 10)
for(r in 1:10){
  for(i in 1:100){
    pval.mgc[i,r] = print.stat.optimal(contam.mgc[[r]][[i]], 4)$pvalue
    pval.mcorr[i,r] = print.stat.optimal(contam.mcorr[[r]][[i]], 4)$pvalue
    pval.hhg[i,r] = print.stat.optimal(contam.hhg[[r]][[i]], 4)$pvalue
    pval.fh[i,r] = contam.fh[i,r]
  }
}

############### making a figure ######################
pdf("../Figure/multi_brain.pdf", width = 15, height = 6)
par(mfrow = c(1,1), cex.lab = 5, cex.axis = 3,
    mar = c(8,10,3,20), tcl = 0.5)
plot(seq(0,100,10), c(1, colMeans(pval.mgc <= 0.05)), col = "red", 
     lty = 1, lwd = 5, 
     ylim = c(0, 1.0), type = "l", mgp = c(6,2,0),
     xlab = "contamination (%)", yaxt = "n", ylab = "Power")
axis(side = 2, at = c(0.0, 0.25, 0.5, 0.75, 1.0), 
     labels = c(0.0, 0.25, 0.5, 0.75,  1.0), 
     tck = 0.05)
lines(seq(0,100,10),c(1, colMeans(pval.mcorr <= 0.05)), col = "dodgerblue", 
      lty = 2, lwd = 5, yaxis = NULL, type = "l")
lines(seq(0,100,10), c(1, colMeans(pval.hhg <= 0.05)), col = "lightsalmon4", 
      lty = 5, lwd = 5, yaxis = NULL, type = "l")
lines(seq(0,100,10), c(1, colMeans(pval.fh <= 0.05)), col = "darkgreen", 
      lty = 4, lwd = 5, yaxis = NULL,  type = "l")
legend("topright", inset=c(-0.43, 0.4), 
       c(expression(MGC %.% DM), expression(dCorr %.% DM), 
         expression(HHG %.% DM), "FH Test"), 
       seg.len = 2,
       col = c("red", "dodgerblue", "lightsalmon4", "darkgreen"),
       lty = c(1,2,5,4), lwd = 4, bty = 'n', xpd = NA, cex = 2.5)
dev.off() 