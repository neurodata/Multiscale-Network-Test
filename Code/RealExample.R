source("diffmaps.R")
source("NetworkTest.R")
source("elbowmap.R")

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
mgc.pval = matrix(0, nrow = 100, ncol = 10)
mcorr.pval = matrix(0, nrow = 100, ncol = 10)
hhg.pval = matrix(0, nrow = 100, ncol = 10)
fh.pval = matrix(0, nrow = 100, ncol = 10)

q.choice = matrix(0, nrow = 100, ncol = 10)
k.choice = matrix(0, nrow = 100, ncol = 10)
contam.G.adj = list()

for(ii in 1:100){
  ## ii-th iteration
  for(nr in 1:10){

  # the number of deleted edges
  n.delete = as.integer(nr*0.1*length(E(SubBrain))) 
  
  # randomly sample deleted edges and deleted them from the graph
  set.seed(ii*nr)
  delete = sample(E(SubBrain), n.delete)
  delete.G = SubBrain - delete
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
  while(sum(colSums(contam.G.adj[[nr]]) == 0 ) > 0){
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
  diffusion.q[ii,nr]  = max(getElbows(print.lambda(contam.P, times = 3)[[1]], plot = FALSE, n = 3, threshold = 0))
  
  mgc.pval[ii,nr] =  NetworkTest.q(G, X, option = 1, diffusion = TRUE, dstep = 3, 
                              n.perm = 500, q = diffusion.q[nr])[[1]][[1]]
  
  mcorr.pval[ii,nr] =  NetworkTest.q(G, X, option = 2, diffusion = TRUE, dstep = 3, 
                                n.perm = 500, q = diffusion.q[nr])
  
  hhg.pval[ii,nr] =  NetworkTest.q(G, X, option = 3, diffusion = TRUE, dstep = 3, 
                              n.perm = 500, q = diffusion.q[nr])
  
  SVD.A = svd(contam.G.adj[[nr]], nu = 95, nv = 95, LINPACK = FALSE)
  
  # FH
  fh.k[nr] =  min( max(getElbows(SVD.A$d, n = 3, plot = FALSE)), 95-1 )
  fh.pval[nr] =  FH_test(contam.G.adj[[nr]], X, k.range = fh.k[nr], n.iter = 500)

  }  

}

############### save the result ############
RealAnalysis[[1]] = mgc.pval
RealAnalysis[[2]] = mcorr.pval
RealAnalysis[[3]] = hhg.pval
RealAnalysis[[4]] = fh.pval

save(RealAnalysis, file = "../Data/RealAnalysis.RData")


mgc.pval.mean = colMeans(mgc.pval)
mcorr.pval.mean = colMeans(mcorr.pval)
hhg.pval.mean = colMeans(hhg.pval)
fh.pval.mean = colSums(fh.pval) 
############### making a figure ######################
pdf("../Figure/Elbow3_t3.pdf", width = 15, height = 6)
par(mfrow = c(1,1), cex.lab = 5, cex.axis = 3,
    mar = c(8,10,3,20), tcl = 0.5)
plot(seq(0,100,10), c(mgc.pval0, mgc.pval.mean), col = "red", 
     lty = 1, lwd = 5, 
     ylim = c(0, 1.0), type = "l", mgp = c(6,2,0),
     xlab = "contamination (%)", yaxt = "n", ylab = "p-value")
axis(side = 2, at = c(0.0, 0.5, 1.0), 
     labels = c(0.0, 0.5, 1.0), 
     tck = 0.05)
lines(seq(0,100,10), c(dcorr.pval0, mcorr.pval.mean), col = "dodgerblue", 
      lty = 2, lwd = 5, yaxis = NULL, type = "l")
lines(seq(0,100,10), c(hhg.pval0, hhg.pval.mean), col = "hotpink", 
      lty = 3, lwd = 5, yaxis = NULL, type = "l")
lines(seq(0,100,10), c(fh.pval0, fh.pval.mean), col = "darkgreen", 
      lty = 4, lwd = 5, yaxis = NULL,  type = "l")
legend("topright", inset=c(-0.43, 0.4), 
       c(expression(MGC %.% DM), expression(mCorr %.% DM), 
         expression(HHG %.% DM), "FH Test"), 
       seg.len = 2,
       col = c("red", "dodgerblue", "hotpink", "darkgreen"),
       lty = c(1,2,3,4), lwd = 4, bty = 'n', xpd = NA, cex = 2.5)
dev.off() 
