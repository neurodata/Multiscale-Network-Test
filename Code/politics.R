########## LOAD DATA ##########
library(igraph)
library(gplots)
library(ColorPalette)
library(RColorBrewer)
library(DTMCPack)
library(Matrix)
library(energy)
library(HHG)
library(ecodist)
library(SDMTools)
library(amen)
library(stats)
library(car)

#############################
## FH test with multivariate-covariates ####
FH_test2 = function(A, X, k.range, n.iter){
  # this function test independence between network topology and nodal attributes
  # via Fosdick & Hoff (2015) 's method of testing where A is binary
  ################################################
  # # input
  # A : binary adjacency matrix of network 
  # X : any dimensional matrix of nodal attributes
  # k.range : range of dimension of multiplicative factor model
  # n.iter : the number of iteration beyond burn-in
  ################################################
  # require : amen
  #################################################
  # # output
  # pvalues : pvalues of each model of different dimension of factors
  #           pavlues are based on Wilks' statistic 
  ################################################
  pvalues <- c()
  
  for(k in 1:length(k.range)){
    
    fit <- ame(A, burn = 1, nscan = n.iter, odens = 1,
               model = "bin", print = FALSE, rvar = TRUE, cvar = TRUE,
               R = k.range[k] , plot = FALSE, symmetric = TRUE) 
    factors <- cbind(fit$APM, fit$U)
    
    # regress using Multivariate ANOVE
    anov <-lm(factors ~ as.factor(X))
    
    
    if(ncol(factors) == 1){
      pvalues[k] <- Anova(anov)[1,4]
    }else{    
      outtests <- car:::print.Anova.mlm
      
      # allow the function to return the results and disable print
      body(outtests)[[16]] <- quote(invisible(tests))
      body(outtests)[[15]] <- NULL
      
      # Now run the regression
      mod <- lm(factors ~ as.factor(X))
      
      # Run the Anova over all tests  
      tab <- lapply(c("Pillai", "Wilks", "Roy"), 
                    function(i)  outtests(Anova(mod, test.statistic=i)))
      
      tab <- do.call(rbind, tab)
      row.names(tab) <- c("Pillai", "Wilks", "Roy")
      pvalues[k] =  tab[2,6]  
    }
    
    
  }
  
  return(pvalues)}

collab0 <- as.matrix(read.table(file = "../Data/politics/climate9500-collab.txt", 
                               header = TRUE, row.names = 1, sep = ";")) ## 34 x 34 matrix (not symmetric)

# type of organization; vector with five character types
types0 <- as.character(read.table(file="../Data/politics/climate9500-type.txt", 
                                 header = TRUE, row.names = 1, sep = ";")[,2]) ## categorical 5 types

# collaboration; directed network
collab <- as.matrix(read.table(file = "../Data/politics/climate0205-collab.txt", 
                               header = TRUE, row.names = 1, sep = ";")) ## 34 x 34 matrix (not symmetric)

# type of organization; vector with five character types
types <- as.character(read.table(file="../Data/politics/climate0205-type.txt", 
                                 header = TRUE, row.names = 1, sep = ";")[,2]) ## categorical 5 types


######## transform directed network into undirected network #########
undirected.G = graph.adjacency(collab, "undirected")
Dy.types = matrix(0, 34, 34)
for(i in 1:34){
  for(j in 1:34){
    Dy.types[i,j] = ifelse(types[i] == types[j], 0, 1)
  }
}


dstep = c(1:10) # diffusion times
n.perm = 500 # the number of permutations
n.iter = n.perm # the number of iterations (FH test)
k.range = c(0:10) # the dimension of multivariate factors in FH model

G.adj = as.matrix(get.adjacency(undirected.G))
G.P = G.adj / rowSums(G.adj)
G.dmaps = dmap(G.P ,dstep)

Dx.G = list();
mcorr.types <- c(); hhg.types <- c(); fh.types <- c()
mgc.types <- list(); 

for(s in 1:length(dstep)){
  print(s)
  Dx.G[[s]] = as.matrix(dist(G.dmaps[[s]], diag = TRUE, upper = TRUE)) 
  ## English graph vs. French text
  mcorr.types[s] <- dcor.ttest(Dx.G[[s]], Dy.types, distance = TRUE)$p.value
  hhg.types[s] <- hhg.test(Dx.G[[s]], Dy.types, nr.perm = n.perm)$perm.pval.hhg.sl
  mgc.types[[s]] <- MGCPermutationTest(Dx.G[[s]], Dy.types, rep = n.perm, option = 'mcor')
}

#### Latent space model
fh.types <- FH_test2(G.adj, types, k.range = c(0:10), n.iter = 300)

s = 1
mgc.results<- MGCPermutationTest(Dx.G[[s]], Dy.types, n.perm, option = 'mcor')
ind = mgc.results[[5]][1]
if(is.na(ind)){
  node.weight <- MGCLocalCorr(Dx.G[[s]], Dy.types, option = 'mcor', ind = length(mgc.results[[4]]))}else{
  node.weight <- MGCLocalCorr(Dx.G[[s]], Dy.types, option = 'mcor', ind = ind)}

if(length(node.weight)  > 1){
  col.weight <- colSums(node.weight)
  rankings <- order(abs(col.weight), decreasing = TRUE)
}

####### with node ccontribution #######################
pdf("../Figure/politic.pdf", width = 8, height = 5)
par(mfrow = c(1,1), mar = c(0,0,2,2), cex.main = 2)
set.seed(123)
igraph.options(vertex.size = 20 - rankings/2  , edge.arrow.size = 0.5)
V(undirected.G)$label = V(undirected.G)$name
V(undirected.G)$label.cex = rep(2,34)
V(undirected.G)$color[types == "private"] = "red"
V(undirected.G)$color[types == "ngo"] = "darkgreen"
V(undirected.G)$color[types == "party"] = "hotpink"
V(undirected.G)$color[types == "gov"] = "dodgerblue"
V(undirected.G)$color[types == "science"] = "gold"
#V(undirected.G)$color[rankings] =  RedPalette
plot(undirected.G, 
     layout= layout.fruchterman.reingold(undirected.G),
     main = "Political Network and organization type", vertex.label = "")
legend("topright", inset=c(0,0.5), cex = 1.2, bty = 'n', xpd = NA, 
       col = c("red", "darkgreen", "hotpink", "dodgerblue", "gold"), 
       c("Business", "NGO", "Party", "State", "Science"), pch = 19, pt.cex = 2)
dev.off()

###### without node contribution ####################
pdf("../Figure/introplot.pdf", width = 8, height = 5)
par(mfrow = c(1,1), mar = c(0,0,2,2), cex.main = 2)
set.seed(123)
igraph.options(vertex.size = 12, edge.arrow.size = 0.5)
V(undirected.G)$label = V(undirected.G)$name
V(undirected.G)$label.cex = rep(2,34)
V(undirected.G)$color[types == "private"] = "red"
V(undirected.G)$color[types == "ngo"] = "darkgreen"
V(undirected.G)$color[types == "party"] = "hotpink"
V(undirected.G)$color[types == "gov"] = "dodgerblue"
V(undirected.G)$color[types == "science"] = "gold"
#V(undirected.G)$color[rankings] =  RedPalette
plot(undirected.G, 
     layout= layout.fruchterman.reingold(undirected.G),
     main = "Collaborative networks and organization types", vertex.label = "")
legend("topright", inset=c(0,0.5), cex = 1.5, bty = 'n', xpd = NA, 
       col = c("red", "darkgreen", "hotpink", "dodgerblue", "gold"), 
       c("Business", "NGO", "Party", "State", "Science"), pch = 19, pt.cex = 2)
dev.off()


undirected.H = graph.adjacency(collab0, "undirected")
Dy.types0 = matrix(0, 34, 34)
for(i in 1:34){
  for(j in 1:34){
    Dy.types0[i,j] = ifelse(types0[i] == types0[j], 0, 1)
  }
}

ind = which(components(undirected.H)$membership == 1)
H.sub = induced_subgraph(undirected.H, ind)

H.sub.adj = as.matrix(get.adjacency(H.sub))


H.sub.P = H.sub.adj / rowSums(H.sub.adj)
H.sub.dmaps = dmap(H.sub.P ,dstep)

Dx.H = list();
mcorr.types0 <- c(); hhg.types0 <- c(); fh.types0 <- c()
mgc.types0 <- list(); 

Dy.sub.types0 = Dy.types0[ind,ind]

for(s in 1:length(dstep)){
  print(s)
  Dx.H[[s]] = as.matrix(dist(H.sub.dmaps[[s]], diag = TRUE, upper = TRUE)) 
  ## English graph vs. French text
  mcorr.types0[s] <- dcor.ttest(Dx.H[[s]], Dy.sub.types0, distance = TRUE)$p.value
  hhg.types0[s] <- hhg.test(Dx.H[[s]], Dy.sub.types0, nr.perm = n.perm)$perm.pval.hhg.sl
  mgc.types0[[s]] <- MGCPermutationTest(Dx.H[[s]], Dy.sub.types0, rep = n.perm, option = 'mcor')
}

#### Latent space model
sub.types0 = types0[ind]
fh.types0 <- FH_test2(H.sub.adj, sub.types0, k.range = c(0:10), n.iter = 300)

s = 1
mgc.results <- MGCPermutationTest(Dx.H[[s]], Dy.sub.types0, n.perm, option = 'mcor')
ind = mgc.results[[5]][1]
if(is.na(ind)){
  node.weight0 <- MGCLocalCorr(Dx.H[[s]], Dy.sub.types0, option = 'mcor', ind = length(mgc.results[[4]]))}else{
    node.weight0 <- MGCLocalCorr(Dx.H[[s]], Dy.sub.types0, option = 'mcor', ind = ind)}

if(length(node.weight0)  > 1){
  col.weight0 <- colSums(node.weight0) 
  rankings0 <- order(abs(col.weight0), decreasing = TRUE)
}
##############################################################
re.rankings = rep(34, 34)
ind = which(components(undirected.H)$membership == 1)
re.rankings[ind] = rankings0

####### collaboration network H  ##################
RedPalette <- colorRampPalette(brewer.pal(9,"Reds"))(34)
pdf("../Figure/politic2.pdf", width = 8, height = 5)
par(mfrow = c(1,1), mar = c(0,0,2,2), cex.main = 2)
set.seed(123)
igraph.options(vertex.size = 20 - re.rankings/2  , edge.arrow.size = 0.5)
V(undirected.H)$label = V(undirected.H)$name
V(undirected.H)$label.cex = rep(2,34)
V(undirected.H)$color[types0 == "private"] = "red"
V(undirected.H)$color[types0 == "ngo"] = "darkgreen"
V(undirected.H)$color[types0 == "party"] = "hotpink"
V(undirected.H)$color[types0 == "gov"] = "dodgerblue"
V(undirected.H)$color[types0 == "science"] = "gold"
#V(undirected.G)$color[rankings] =  RedPalette
plot(undirected.H, 
     layout= layout.fruchterman.reingold(undirected.H),
     main = "Political Network and organization type", vertex.label = "")
legend("topright", inset=c(0,0.5), cex = 1.2, bty = 'n', xpd = NA, 
       col = c("red", "darkgreen", "hotpink", "dodgerblue", "gold"), 
       c("Business", "NGO", "Party", "State", "Science"), pch = 19, pt.cex = 2)
dev.off()


####### two political networks together ##############
pdf("../Figure/two_politics.pdf", width = 10, height = 5)
par(mfrow = c(1,2), mar = c(2,0,0,0), cex.main = 2,
    oma = c(3, 3, 3, 10), mai = c(0, 0, 0, 0), xpd = NA)
set.seed(123)

igraph.options(vertex.size = 20 - re.rankings/2  , edge.arrow.size = 0.5, xpd = NA)
V(undirected.H)$label = V(undirected.H)$name
V(undirected.H)$label.cex = rep(2,34)
V(undirected.H)$color[types0 == "private"] = "red"
V(undirected.H)$color[types0 == "ngo"] = "darkgreen"
V(undirected.H)$color[types0 == "party"] = "hotpink"
V(undirected.H)$color[types0 == "gov"] = "dodgerblue"
V(undirected.H)$color[types0 == "science"] = "gold"
#V(undirected.G)$color[rankings] =  RedPalette
plot(undirected.H, 
     layout= layout.fruchterman.reingold(undirected.H),
     xlab = "first period", vertex.label = "")

igraph.options(vertex.size = 20 - rankings/2  , edge.arrow.size = 0.5, xpd = NA)
V(undirected.G)$label = V(undirected.G)$name
V(undirected.G)$label.cex = rep(2,34)
V(undirected.G)$color[types == "private"] = "red"
V(undirected.G)$color[types == "ngo"] = "darkgreen"
V(undirected.G)$color[types == "party"] = "hotpink"
V(undirected.G)$color[types == "gov"] = "dodgerblue"
V(undirected.G)$color[types == "science"] = "gold"
#V(undirected.G)$color[rankings] =  RedPalette
plot(undirected.G, 
     layout= layout.fruchterman.reingold(undirected.G),
     xlab = "second period", vertex.label = "")

legend("topright", inset=c(-0.5,0.5), cex = 1.2, bty = 'n', xpd = NA, 
       col = c("red", "darkgreen", "hotpink", "dodgerblue", "gold"), 
       c("Business", "NGO", "Party", "State", "Science"), pch = 19, pt.cex = 2)
mtext("Collaborative networks and organization types", side = 3, outer = TRUE, cex = 2, xpd = NA)
mtext("First Period", adj = 0.2, side = 1, outer = TRUE, cex = 2, xpd = NA)
mtext("Second Period", adj = 1, side = 1, outer = TRUE, cex = 2, xpd = NA)
dev.off()