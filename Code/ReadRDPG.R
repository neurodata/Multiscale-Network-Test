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

source("notsimple3block.R")
source("diffmaps.R")
source("FH_test.R")
source("getElbows.R")
source("NetworkTest.R")
source("optimalT.R")

M = 100; popn = 50
n.perm = 500; n.iter = n.perm

multi_linear50 = list(); multi_exp50 = list()
multi_cubic50 = list(); multi_jointN50 = list()
multi_step50 = list(); multi_quad50 = list()
multi_W50 = list(); multi_spiral50 = list()
multi_bern50 = list(); multi_log50 = list()
multi_root50 = list(); multi_sine4pi50 = list()
multi_sine16pi50 = list(); multi_square50 = list()
multi_twopara50 = list(); multi_circle50 = list()
multi_ellipse50 = list(); multi_diamond50 = list()
multi_multi50 = list(); multi_ind50 = list()

## 1. linear
for(i in 1:M){
  set.seed(i)
  W = runif(popn, 0, 1)
  X = W + rnorm(popn, 0, 0.5)
  W = (W- min(W)) / (max(W) - min(W))
  X = (X - min(X)) / (max(X) - min(X))
  A = matrix(0, popn, popn)
  for(i in 1: (popn-1)){
    for(j in (i+1):popn){
      p =  sum(W[i]*W[j])
      A[i,j] = rbinom(1,1,p); A[j,i] = A[i,j]
  }}
  G = graph.adjacency(A, "undirected")
  mgc.results =  Networktest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  mcorr.results =  Networktest.diffusion.stat(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  hhg.results =  Networktest.diffusion.stat(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  # estimate the matrix M based on SVD
  SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
  # FH
  fh.k = min(max(abs(getElbows(SVD.A$d, n = 2, plot = FALSE))), popn-10)
  fh.results = try(FH_test(A, X, k.range = fh.k, n.iter = n.iter), silent = TRUE)
  multi_linear50[[i]] =  list(mgc.results, mcorr.results, hhg.results, fh.results) 
}

## 2. exponential
for(i in 1:M){
  set.seed(i)
  W = runif(popn, 0, 3)
  X = exp(W) + rnorm(popn, 0, 5)
  W = (W- min(W)) / (max(W) - min(W))
  X = (X - min(X)) / (max(X) - min(X))
  A = matrix(0, popn, popn)
  for(i in 1: (popn-1)){
    for(j in (i+1):popn){
      p =  sum(W[i]*W[j])
      A[i,j] = rbinom(1,1,p); A[j,i] = A[i,j]
  }}
  G = graph.adjacency(A, "undirected")
  mgc.results =  Networktest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  mcorr.results =  Networktest.diffusion.stat(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  hhg.results =  Networktest.diffusion.stat(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  # estimate the matrix M based on SVD
  SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
  # FH
  fh.k = min(max(abs(getElbows(SVD.A$d, n = 2, plot = FALSE))), popn-10)
  fh.results = try(FH_test(A, X, k.range = fh.k, n.iter = n.iter), silent = TRUE)
  multi_exp50[[i]] =  list(mgc.results, mcorr.results, hhg.results, fh.results) 
}

## 3. cubic
for(i in 1:M){
  set.seed(i)
  W = runif(popn, 0, 1)
  X = 20*(W - 1/2)^3 + 2*(W - 1/2)^2 - 1*(W - 1/2) + rnorm(popn, 0, 0.5)
  W = (W- min(W)) / (max(W) - min(W))
  X = (X - min(X)) / (max(X) - min(X))
  A = matrix(0, popn, popn)
  for(i in 1: (popn-1)){
    for(j in (i+1):popn){
      p =  sum(W[i]*W[j])
      A[i,j] = rbinom(1,1,p); A[j,i] = A[i,j]
  }}
  G = graph.adjacency(A, "undirected")
  mgc.results =  Networktest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  mcorr.results =  Networktest.diffusion.stat(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  hhg.results =  Networktest.diffusion.stat(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  # estimate the matrix M based on SVD
  SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
  # FH
  fh.k = min(max(abs(getElbows(SVD.A$d, n = 2, plot = FALSE))), popn-10)
  fh.results = try(FH_test(A, X, k.range = fh.k, n.iter = n.iter), silent = TRUE)
  multi_cubic50[[i]] =  list(mgc.results, mcorr.results, hhg.results, fh.results) 
}


## 4. joint normal
for(i in 1:M){
  set.seed(i)
  tmp = mvrnorm(popn, mu = c(0,0), Sigma = matrix(c(0.7, 0.5, 0.5, 0.7), 2,2)) 
  W = tmp[,1]; X = tmp[,2]
  W = (W- min(W)) / (max(W) - min(W))
  X = (X - min(X)) / (max(X) - min(X))
  A = matrix(0, popn, popn)
  for(i in 1: (popn-1)){
    for(j in (i+1):popn){
      p =  sum(W[i]*W[j])
      A[i,j] = rbinom(1,1,p); A[j,i] = A[i,j]
  }}
  G = graph.adjacency(A, "undirected")
  mgc.results =  Networktest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  mcorr.results =  Networktest.diffusion.stat(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  hhg.results =  Networktest.diffusion.stat(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  # estimate the matrix M based on SVD
  SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
  # FH
  fh.k = min(max(abs(getElbows(SVD.A$d, n = 2, plot = FALSE))), popn-10)
  fh.results = try(FH_test(A, X, k.range = fh.k, n.iter = n.iter), silent = TRUE)
  multi_jointN50[[i]] =  list(mgc.results, mcorr.results, hhg.results, fh.results) 
}

## 5. step normal
for(i in 1:M){
  set.seed(i)
  W = runif(popn, -1, 1)
  X = (W > 0) + rnorm(popn, 0, 0.5)
  W = (W- min(W)) / (max(W) - min(W))
  X = (X - min(X)) / (max(X) - min(X))
  A = matrix(0, popn, popn)
  for(i in 1: (popn-1)){
    for(j in (i+1):popn){
      p =  sum(W[i]*W[j])
      A[i,j] = rbinom(1,1,p); A[j,i] = A[i,j]
  }}
  G = graph.adjacency(A, "undirected")
  mgc.results =  Networktest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  mcorr.results =  Networktest.diffusion.stat(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  hhg.results =  Networktest.diffusion.stat(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  # estimate the matrix M based on SVD
  SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
  # FH
  fh.k = min(max(abs(getElbows(SVD.A$d, n = 2, plot = FALSE))), popn-10)
  fh.results = try(FH_test(A, X, k.range = fh.k, n.iter = n.iter), silent = TRUE)
  multi_step50[[i]] =  list(mgc.results, mcorr.results, hhg.results, fh.results) 
}


## 6. quadratic function
for(i in 1:M){
  set.seed(i)
  W = runif(popn, -1, 1)
  X = W^2 + rnorm(popn, 0, 0.3)
  W = (W- min(W)) / (max(W) - min(W))
  X = (X - min(X)) / (max(X) - min(X))
  A = matrix(0, popn, popn)
  for(i in 1: (popn-1)){
    for(j in (i+1):popn){
      p =  sum(W[i]*W[j])
      A[i,j] = rbinom(1,1,p); A[j,i] = A[i,j]
  }}
  G = graph.adjacency(A, "undirected")
  mgc.results =  Networktest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  mcorr.results =  Networktest.diffusion.stat(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  hhg.results =  Networktest.diffusion.stat(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  # estimate the matrix M based on SVD
  SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
  # FH
  fh.k = min(max(abs(getElbows(SVD.A$d, n = 2, plot = FALSE))), popn-10)
  fh.results = try(FH_test(A, X, k.range = fh.k, n.iter = n.iter), silent = TRUE)
  multi_quad50[[i]] =  list(mgc.results, mcorr.results, hhg.results, fh.results) 
}

## 7. W function
for(i in 1:M){
  set.seed(i)
  W = runif(popn, -1, 1)
  X = 4*( (W^2 - 0.5)^2 )
  W = (W- min(W)) / (max(W) - min(W))
  X = (X - min(X)) / (max(X) - min(X))
  A = matrix(0, popn, popn)
  for(i in 1: (popn-1)){
    for(j in (i+1):popn){
      p =  sum(W[i]*W[j])
      A[i,j] = rbinom(1,1,p); A[j,i] = A[i,j]
  }}
  G = graph.adjacency(A, "undirected")
  mgc.results =  Networktest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  mcorr.results =  Networktest.diffusion.stat(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  hhg.results =  Networktest.diffusion.stat(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  # estimate the matrix M based on SVD
  SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
  # FH
  fh.k = min(max(abs(getElbows(SVD.A$d, n = 2, plot = FALSE))), popn-10)
  fh.results = try(FH_test(A, X, k.range = fh.k, n.iter = n.iter), silent = TRUE)
  multi_W50[[i]] =  list(mgc.results, mcorr.results, hhg.results, fh.results) 
}


## 8. spiral
for(i in 1:M){
  set.seed(i)
  U1 = runif(popn, 0, 5)
  W = U1*cos(pi*U1)
  X = U1*sin(pi*U1) + rnorm(popn, 0, 0.1)
  W = (W- min(W)) / (max(W) - min(W))
  X = (X - min(X)) / (max(X) - min(X))
  A = matrix(0, popn, popn)
  for(i in 1: (popn-1)){
    for(j in (i+1):popn){
      p =  sum(W[i]*W[j])
      A[i,j] = rbinom(1,1,p); A[j,i] = A[i,j]
  }}
  G = graph.adjacency(A, "undirected")
  mgc.results =  Networktest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  mcorr.results =  Networktest.diffusion.stat(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  hhg.results =  Networktest.diffusion.stat(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  # estimate the matrix M based on SVD
  SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
  # FH
  fh.k = min(max(abs(getElbows(SVD.A$d, n = 2, plot = FALSE))), popn-10)
  fh.results = try(FH_test(A, X, k.range = fh.k, n.iter = n.iter), silent = TRUE)
  multi_spiral50[[i]] =  list(mgc.results, mcorr.results, hhg.results, fh.results) 
}

## 9. Bernoulli
for(i in 1:M){
  set.seed(i)
  W = rbinom(popn, 1, 0.5)
  X = (2*rbinom(popn, 1, 0.5) - 1)*W + rnorm(popn, 0, 1)
  W = (W- min(W)) / (max(W) - min(W))
  X = (X - min(X)) / (max(X) - min(X))
  A = matrix(0, popn, popn)
  for(i in 1: (popn-1)){
    for(j in (i+1):popn){
      p =  sum(W[i]*W[j])
      A[i,j] = rbinom(1,1,p); A[j,i] = A[i,j]
  }}
  G = graph.adjacency(A, "undirected")
  mgc.results =  Networktest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  mcorr.results =  Networktest.diffusion.stat(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  hhg.results =  Networktest.diffusion.stat(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  # estimate the matrix M based on SVD
  SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
  # FH
  fh.k = min(max(abs(getElbows(SVD.A$d, n = 2, plot = FALSE))), popn-10)
  fh.results = try(FH_test(A, X, k.range = fh.k, n.iter = n.iter), silent = TRUE)
  multi_bern50[[i]] =  list(mgc.results, mcorr.results, hhg.results, fh.results) 
}


## 10. log
for(i in 1:M){
  set.seed(i)
  W = runif(popn, -1, 1)
  X = 5*log2(abs(W)) + rnorm(popn, 0, 5) 
  W = (W- min(W)) / (max(W) - min(W))
  X = (X - min(X)) / (max(X) - min(X))
  A = matrix(0, popn, popn)
  for(i in 1: (popn-1)){
    for(j in (i+1):popn){
      p =  sum(W[i]*W[j])
      A[i,j] = rbinom(1,1,p); A[j,i] = A[i,j]
  }}
  G = graph.adjacency(A, "undirected")
  mgc.results =  Networktest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  mcorr.results =  Networktest.diffusion.stat(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  hhg.results =  Networktest.diffusion.stat(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  # estimate the matrix M based on SVD
  SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
  # FH
  fh.k = min(max(abs(getElbows(SVD.A$d, n = 2, plot = FALSE))), popn-10)
  fh.results = try(FH_test(A, X, k.range = fh.k, n.iter = n.iter), silent = TRUE)
  multi_log50[[i]] =  list(mgc.results, mcorr.results, hhg.results, fh.results) 
}

## 11. root
for(i in 1:M){
  set.seed(i)
  W = runif(popn, 0, 1)
  X = (abs(W + rnorm(popn, 0, 0.5)))^(1/4) 
  W = (W- min(W)) / (max(W) - min(W))
  X = (X - min(X)) / (max(X) - min(X))
  A = matrix(0, popn, popn)
  for(i in 1: (popn-1)){
    for(j in (i+1):popn){
      p =  sum(W[i]*W[j])
      A[i,j] = rbinom(1,1,p); A[j,i] = A[i,j]
  }}
  G = graph.adjacency(A, "undirected")
  mgc.results =  Networktest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  mcorr.results =  Networktest.diffusion.stat(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  hhg.results =  Networktest.diffusion.stat(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  # estimate the matrix M based on SVD
  SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
  # FH
  fh.k = min(max(abs(getElbows(SVD.A$d, n = 2, plot = FALSE))), popn-10)
  fh.results = try(FH_test(A, X, k.range = fh.k, n.iter = n.iter), silent = TRUE)
  multi_root50[[i]] =  list(mgc.results, mcorr.results, hhg.results, fh.results) 
}

## 12. sine 4pi
for(i in 1:M){
  set.seed(i)
  W = runif(popn, -1, 1) 
  X = sin(4*pi*W) + 0.01*rnorm(popn, 0, 1)
  W = (W- min(W)) / (max(W) - min(W))
  X = (X - min(X)) / (max(X) - min(X))
  A = matrix(0, popn, popn)
  for(i in 1: (popn-1)){
    for(j in (i+1):popn){
      p =  sum(W[i]*W[j])
      A[i,j] = rbinom(1,1,p); A[j,i] = A[i,j]
  }}
  G = graph.adjacency(A, "undirected")
  mgc.results =  Networktest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  mcorr.results =  Networktest.diffusion.stat(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  hhg.results =  Networktest.diffusion.stat(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  # estimate the matrix M based on SVD
  SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
  # FH
  fh.k = min(max(abs(getElbows(SVD.A$d, n = 2, plot = FALSE))), popn-10)
  fh.results = try(FH_test(A, X, k.range = fh.k, n.iter = n.iter), silent = TRUE)
  multi_sine4pi50[[i]] =  list(mgc.results, mcorr.results, hhg.results, fh.results) 
}

## 13. sine 16pi
for(i in 1:M){
  set.seed(i)
  W = runif(popn, -1, 1) 
  X = sin(16*pi*W) + 0.01*rnorm(popn, 0, 1)
  W = (W- min(W)) / (max(W) - min(W))
  X = (X - min(X)) / (max(X) - min(X))
  A = matrix(0, popn, popn)
  for(i in 1: (popn-1)){
    for(j in (i+1):popn){
      p =  sum(W[i]*W[j])
      A[i,j] = rbinom(1,1,p); A[j,i] = A[i,j]
  }}
  G = graph.adjacency(A, "undirected")
  mgc.results =  Networktest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  mcorr.results =  Networktest.diffusion.stat(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  hhg.results =  Networktest.diffusion.stat(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  # estimate the matrix M based on SVD
  SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
  # FH
  fh.k = min(max(abs(getElbows(SVD.A$d, n = 2, plot = FALSE))), popn-10)
  fh.results = try(FH_test(A, X, k.range = fh.k, n.iter = n.iter), silent = TRUE)
  multi_sine16pi50[[i]] =  list(mgc.results, mcorr.results, hhg.results, fh.results) 
}

## 14. square
for(i in 1:M){
  set.seed(i)
  U1 = runif(popn, -1, 1); U2 = runif(popn, -1, 1)
  W = U1*cos(-pi/8) + U2*sin(-pi/8) 
  X = -U1*sin(-pi/8) + U2*cos(-pi/8)
  W = (W- min(W)) / (max(W) - min(W))
  X = (X - min(X)) / (max(X) - min(X))
  A = matrix(0, popn, popn)
  for(i in 1: (popn-1)){
    for(j in (i+1):popn){
      p =  sum(W[i]*W[j])
      A[i,j] = rbinom(1,1,p); A[j,i] = A[i,j]
  }}
  G = graph.adjacency(A, "undirected")
  mgc.results =  Networktest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  mcorr.results =  Networktest.diffusion.stat(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  hhg.results =  Networktest.diffusion.stat(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  # estimate the matrix M based on SVD
  SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
  # FH
  fh.k = min(max(abs(getElbows(SVD.A$d, n = 2, plot = FALSE))), popn-10)
  fh.results = try(FH_test(A, X, k.range = fh.k, n.iter = n.iter), silent = TRUE)
  multi_square50[[i]] =  list(mgc.results, mcorr.results, hhg.results, fh.results) 
}


## 15. two parabola
for(i in 1:M){
  set.seed(i)
  W = runif(popn, 0, 1) 
  X = (W^2 + rnorm(popn, 0.5, 0.3))*(rbinom(popn, 1, 0.3) - 0.5)
  W = (W- min(W)) / (max(W) - min(W))
  X = (X - min(X)) / (max(X) - min(X))
  A = matrix(0, popn, popn)
  for(i in 1: (popn-1)){
    for(j in (i+1):popn){
      p =  sum(W[i]*W[j])
      A[i,j] = rbinom(1,1,p); A[j,i] = A[i,j]
  }}
  G = graph.adjacency(A, "undirected")
  mgc.results =  Networktest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  mcorr.results =  Networktest.diffusion.stat(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  hhg.results =  Networktest.diffusion.stat(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  # estimate the matrix M based on SVD
  SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
  # FH
  fh.k = min(max(abs(getElbows(SVD.A$d, n = 2, plot = FALSE))), popn-10)
  fh.results = try(FH_test(A, X, k.range = fh.k, n.iter = n.iter), silent = TRUE)
  multi_twopara50[[i]] =  list(mgc.results, mcorr.results, hhg.results, fh.results) 
}

## 16. circle
for(i in 1:M){
  set.seed(i)
  U1 = runif(popn, -1, 1)
  W = cos(pi*U1) 
  X = sin(pi*U1) + rnorm(popn, 0, 0.05)
  W = (W- min(W)) / (max(W) - min(W))
  X = (X - min(X)) / (max(X) - min(X))
  A = matrix(0, popn, popn)
  for(i in 1: (popn-1)){
    for(j in (i+1):popn){
      p =  sum(W[i]*W[j])
      A[i,j] = rbinom(1,1,p); A[j,i] = A[i,j]
  }}
  G = graph.adjacency(A, "undirected")
  mgc.results =  Networktest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  mcorr.results =  Networktest.diffusion.stat(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  hhg.results =  Networktest.diffusion.stat(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  # estimate the matrix M based on SVD
  SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
  # FH
  fh.k = min(max(abs(getElbows(SVD.A$d, n = 2, plot = FALSE))), popn-10)
  fh.results = try(FH_test(A, X, k.range = fh.k, n.iter = n.iter), silent = TRUE)
  multi_circle50[[i]] =  list(mgc.results, mcorr.results, hhg.results, fh.results) 
}

## 17. ellipse
for(i in 1:M){
  set.seed(i)
  U1 = runif(popn, -1, 1)
  Y = 5*cos(pi*U1); X = sin(pi*U1) 
  W = (Y - min(Y)) / (max(Y) - min(Y))
  X = X / (max(Y) - min(Y)) + 0.5
  A = matrix(0, popn, popn)
  for(i in 1: (popn-1)){
    for(j in (i+1):popn){
      p =  sum(W[i]*W[j])
      A[i,j] = rbinom(1,1,p); A[j,i] = A[i,j]
  }}
  G = graph.adjacency(A, "undirected")
  mgc.results =  Networktest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  mcorr.results =  Networktest.diffusion.stat(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  hhg.results =  Networktest.diffusion.stat(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  # estimate the matrix M based on SVD
  SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
  # FH
  fh.k = min(max(abs(getElbows(SVD.A$d, n = 2, plot = FALSE))), popn-10)
  fh.results = try(FH_test(A, X, k.range = fh.k, n.iter = n.iter), silent = TRUE)
  multi_ellipse50[[i]] =  list(mgc.results, mcorr.results, hhg.results, fh.results) 
}

## 18. diamond
for(i in 1:M){
  set.seed(i)
  U1 = runif(popn, -1, 1); U2 = runif(popn, -1, 1)
  W = U1*cos(-pi/4) + U2*sin(-pi/4) 
  X = -U1*sin(-pi/4) + U2*cos(-pi/4)
  W = (W- min(W)) / (max(W) - min(W))
  X = (X - min(X)) / (max(X) - min(X))
  A = matrix(0, popn, popn)
  for(i in 1: (popn-1)){
    for(j in (i+1):popn){
      p =  sum(W[i]*W[j])
      A[i,j] = rbinom(1,1,p); A[j,i] = A[i,j]
  }}
  G = graph.adjacency(A, "undirected")
  mgc.results =  Networktest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  mcorr.results =  Networktest.diffusion.stat(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  hhg.results =  Networktest.diffusion.stat(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  # estimate the matrix M based on SVD
  SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
  # FH
  fh.k = min(max(abs(getElbows(SVD.A$d, n = 2, plot = FALSE))), popn-10)
  fh.results = try(FH_test(A, X, k.range = fh.k, n.iter = n.iter), silent = TRUE)
  multi_diamond50[[i]] =  list(mgc.results, mcorr.results, hhg.results, fh.results) 
}

## 19. multiplicative
for(i in 1:M){
  set.seed(i)
  W = rnorm(popn, 0.5, 1)
  X = rnorm(popn, 0.5, 1)*W
  W = (W- min(W)) / (max(W) - min(W))
  X = (X - min(X)) / (max(X) - min(X))
  A = matrix(0, popn, popn)
  for(i in 1: (popn-1)){
    for(j in (i+1):popn){
      p =  sum(W[i]*W[j])
      A[i,j] = rbinom(1,1,p); A[j,i] = A[i,j]
  }}
  G = graph.adjacency(A, "undirected")
  mgc.results =  Networktest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  mcorr.results =  Networktest.diffusion.stat(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  hhg.results =  Networktest.diffusion.stat(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  # estimate the matrix M based on SVD
  SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
  # FH
  fh.k = min(max(abs(getElbows(SVD.A$d, n = 2, plot = FALSE))), popn-10)
  fh.results = try(FH_test(A, X, k.range = fh.k, n.iter = n.iter), silent = TRUE)
  multi_multi50[[i]] =  list(mgc.results, mcorr.results, hhg.results, fh.results) 
}

## 20. independent
for(i in 1:M){
  set.seed(i)
  W = rnorm(popn, 0, 1)
  X = runif(popn, 0, 1)
  W = (W- min(W)) / (max(W) - min(W))
  X = (X - min(X)) / (max(X) - min(X))
  A = matrix(0, popn, popn)
  for(i in 1: (popn-1)){
    for(j in (i+1):popn){
      p =  sum(W[i]*W[j])
      A[i,j] = rbinom(1,1,p); A[j,i] = A[i,j]
  }}
  G = graph.adjacency(A, "undirected")
  mgc.results =  Networktest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  mcorr.results =  Networktest.diffusion.stat(G, X, option = 2, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  hhg.results =  Networktest.diffusion.stat(G, X, option = 3, diffusion = TRUE, t.range = c(0:10), n.perm = n.perm)
  # estimate the matrix M based on SVD
  SVD.A = svd(A, nu = popn, nv = popn, LINPACK = FALSE)
  # FH
  fh.k = min(max(abs(getElbows(SVD.A$d, n = 2, plot = FALSE))), popn-10)
  fh.results = try(FH_test(A, X, k.range = fh.k, n.iter = n.iter), silent = TRUE)
  multi_ind50[[i]] =  list(mgc.results, mcorr.results, hhg.results, fh.results) 
}

save(RDPGdata, file = "../Data/RDPGdata.RData")

multi_linear50 = RDPGdata[[1]]; multi_exp50 = RDPGdata[[2]]
multi_cubic50 = RDPGdata[[3]]; multi_jointN50 = RDPGdata[[4]]
multi_step50 = RDPGdata[[5]]; multi_quad50 = RDPGdata[[6]]
multi_W50 = RDPGdata[[7]]; multi_spiral50 = RDPGdata[[8]]
multi_bern50 = RDPGdata[[9]]; multi_log50 = RDPGdata[[10]]
multi_root50 = RDPGdata[[11]]; multi_sine4pi50 = RDPGdata[[12]]
multi_sine16pi50 = RDPGdata[[13]]; multi_square50 = RDPGdata[[14]]
multi_twopara50 = RDPGdata[[15]]; multi_circle50 = RDPGdata[[16]]
multi_ellipse50 = RDPGdata[[17]]; multi_diamond50 = RDPGdata[[18]]
multi_multi50 = RDPGdata[[19]]; multi_ind50 = RDPGdata[[20]]

mgc.pval = c(); mcorr.pval = c(); hhg.pval = c(); fh.pval = c()
for(i in 1:length(multi_linear50)){
  mgc.pval[i] = print.stat.optimal(multi_linear50[[i]][[1]], 4)$pvalue
  mcorr.pval[i] = print.stat.optimal(multi_linear50[[i]][[2]], 4)$pvalue
  hhg.pval[i] = print.stat.optimal(multi_linear50[[i]][[3]], 4)$pvalue
  fh.pval[i] = multi_linear50[[i]][[4]]
}
linear.mgc = mean(mgc.pval <= 0.05); linear.mcorr = mean(mcorr.pval <= 0.05)
linear.hhg = mean(hhg.pval <= 0.05); linear.fh = mean(fh.pval <= 0.05)

mgc.pval = c(); mcorr.pval = c(); hhg.pval = c(); fh.pval = c()
for(i in 1:length(multi_exp50)){
  mgc.pval[i] = print.stat.optimal(multi_exp50[[i]][[1]], 4)$pvalue
  mcorr.pval[i] = print.stat.optimal(multi_exp50[[i]][[2]], 4)$pvalue
  hhg.pval[i] = print.stat.optimal(multi_exp50[[i]][[3]], 4)$pvalue
  fh.pval[i] = multi_exp50[[i]][[4]]
}
exp.mgc = mean(mgc.pval <= 0.05); exp.mcorr = mean(mcorr.pval <= 0.05)
exp.hhg = mean(hhg.pval <= 0.05); exp.fh = mean(fh.pval <= 0.05)

mgc.pval = c(); mcorr.pval = c(); hhg.pval = c(); fh.pval = c()
for(i in 1:length(multi_cubic50)){
  mgc.pval[i] = print.stat.optimal(multi_cubic50[[i]][[1]], 4)$pvalue
  mcorr.pval[i] = print.stat.optimal(multi_cubic50[[i]][[2]], 4)$pvalue
  hhg.pval[i] = print.stat.optimal(multi_cubic50[[i]][[3]], 4)$pvalue
  fh.pval[i] = multi_cubic50[[i]][[4]]
}
cubic.mgc = mean(mgc.pval <= 0.05); cubic.mcorr = mean(mcorr.pval <= 0.05)
cubic.hhg = mean(hhg.pval <= 0.05); cubic.fh = mean(fh.pval <= 0.05)

mgc.pval = c(); mcorr.pval = c(); hhg.pval = c(); fh.pval = c()
for(i in 1:length(multi_jointN50)){
  mgc.pval[i] = print.stat.optimal(multi_jointN50[[i]][[1]], 4)$pvalue
  mcorr.pval[i] = print.stat.optimal(multi_jointN50[[i]][[2]], 4)$pvalue
  hhg.pval[i] = print.stat.optimal(multi_jointN50[[i]][[3]], 4)$pvalue
  fh.pval[i] = multi_jointN50[[i]][[4]]
}
jointN.mgc = mean(mgc.pval <= 0.05); jointN.mcorr = mean(mcorr.pval <= 0.05)
jointN.hhg = mean(hhg.pval <= 0.05); jointN.fh = mean(fh.pval <= 0.05)

mgc.pval = c(); mcorr.pval = c(); hhg.pval = c(); fh.pval = c()
for(i in 1:length(multi_step50)){
  mgc.pval[i] = print.stat.optimal(multi_step50[[i]][[1]], 4)$pvalue
  mcorr.pval[i] = print.stat.optimal(multi_step50[[i]][[2]], 4)$pvalue
  hhg.pval[i] = print.stat.optimal(multi_step50[[i]][[3]], 4)$pvalue
  fh.pval[i] = multi_step50[[i]][[4]]
}
step.mgc = mean(mgc.pval <= 0.05); step.mcorr = mean(mcorr.pval <= 0.05)
step.hhg = mean(hhg.pval <= 0.05); step.fh = mean(fh.pval <= 0.05)

mgc.pval = c(); mcorr.pval = c(); hhg.pval = c(); fh.pval = c()
for(i in 1:length(multi_quad50)){
  mgc.pval[i] = print.stat.optimal(multi_quad50[[i]][[1]], 4)$pvalue
  mcorr.pval[i] = print.stat.optimal(multi_quad50[[i]][[2]], 4)$pvalue
  hhg.pval[i] = print.stat.optimal(multi_quad50[[i]][[3]], 4)$pvalue
  fh.pval[i] = multi_quad50[[i]][[4]]
}
quad.mgc = mean(mgc.pval <= 0.05); quad.mcorr = mean(mcorr.pval <= 0.05)
quad.hhg = mean(hhg.pval <= 0.05); quad.fh = mean(fh.pval <= 0.05)


mgc.pval = c(); mcorr.pval = c(); hhg.pval = c(); fh.pval = c()
for(i in 1:length(multi_W50)){
  mgc.pval[i] = print.stat.optimal(multi_W50[[i]][[1]], 4)$pvalue
  mcorr.pval[i] = print.stat.optimal(multi_W50[[i]][[2]], 4)$pvalue
  hhg.pval[i] = print.stat.optimal(multi_W50[[i]][[3]], 4)$pvalue
  fh.pval[i] = multi_W50[[i]][[4]]
}
W.mgc = mean(mgc.pval <= 0.05); W.mcorr = mean(mcorr.pval <= 0.05)
W.hhg = mean(hhg.pval <= 0.05); W.fh = mean(fh.pval <= 0.05)

mgc.pval = c(); mcorr.pval = c(); hhg.pval = c(); fh.pval = c()
for(i in 1:length(multi_spiral50)){
  mgc.pval[i] = print.stat.optimal(multi_spiral50[[i]][[1]], 4)$pvalue
  mcorr.pval[i] = print.stat.optimal(multi_spiral50[[i]][[2]], 4)$pvalue
  hhg.pval[i] = print.stat.optimal(multi_spiral50[[i]][[3]], 4)$pvalue
  fh.pval[i] = multi_spiral50[[i]][[4]]
}
spiral.mgc = mean(mgc.pval <= 0.05); spiral.mcorr = mean(mcorr.pval <= 0.05)
spiral.hhg = mean(hhg.pval <= 0.05); spiral.fh = mean(fh.pval <= 0.05)

mgc.pval = c(); mcorr.pval = c(); hhg.pval = c(); fh.pval = c()
for(i in 1:length(multi_bern50)){
  mgc.pval[i] = print.stat.optimal(multi_bern50[[i]][[1]], 4)$pvalue
  mcorr.pval[i] = print.stat.optimal(multi_bern50[[i]][[2]], 4)$pvalue
  hhg.pval[i] = print.stat.optimal(multi_bern50[[i]][[3]], 4)$pvalue
  fh.pval[i] = multi_bern50[[i]][[4]]
}
bern.mgc = mean(mgc.pval <= 0.05); bern.mcorr = mean(mcorr.pval <= 0.05)
bern.hhg = mean(hhg.pval <= 0.05); bern.fh = mean(fh.pval <= 0.05)


mgc.pval = c(); mcorr.pval = c(); hhg.pval = c(); fh.pval = c()
for(i in 1:length(multi_log50)){
  mgc.pval[i] = print.stat.optimal(multi_log50[[i]][[1]], 4)$pvalue
  mcorr.pval[i] = print.stat.optimal(multi_log50[[i]][[2]], 4)$pvalue
  hhg.pval[i] = print.stat.optimal(multi_log50[[i]][[3]], 4)$pvalue
  fh.pval[i] = multi_log50[[i]][[4]]
}
log.mgc = mean(mgc.pval <= 0.05); log.mcorr = mean(mcorr.pval <= 0.05)
log.hhg = mean(hhg.pval <= 0.05); log.fh = mean(fh.pval <= 0.05)

mgc.pval = c(); mcorr.pval = c(); hhg.pval = c(); fh.pval = c()
for(i in 1:length(multi_root50)){
  mgc.pval[i] = print.stat.optimal(multi_root50[[i]][[1]], 4)$pvalue
  mcorr.pval[i] = print.stat.optimal(multi_root50[[i]][[2]], 4)$pvalue
  hhg.pval[i] = print.stat.optimal(multi_root50[[i]][[3]], 4)$pvalue
  fh.pval[i] = multi_root50[[i]][[4]]
}
root.mgc = mean(mgc.pval <= 0.05); root.mcorr = mean(mcorr.pval <= 0.05)
root.hhg = mean(hhg.pval <= 0.05); root.fh = mean(fh.pval <= 0.05)

mgc.pval = c(); mcorr.pval = c(); hhg.pval = c(); fh.pval = c()
for(i in 1:length(multi_sine4pi50)){
  mgc.pval[i] = print.stat.optimal(multi_sine4pi50[[i]][[1]], 4)$pvalue
  mcorr.pval[i] = print.stat.optimal(multi_sine4pi50[[i]][[2]], 4)$pvalue
  hhg.pval[i] = print.stat.optimal(multi_sine4pi50[[i]][[3]], 4)$pvalue
  fh.pval[i] = multi_sine4pi50[[i]][[4]]
}
sine4pi.mgc = mean(mgc.pval <= 0.05); sine4pi.mcorr = mean(mcorr.pval <= 0.05)
sine4pi.hhg = mean(hhg.pval <= 0.05); sine4pi.fh = mean(fh.pval <= 0.05)


mgc.pval = c(); mcorr.pval = c(); hhg.pval = c(); fh.pval = c()
for(i in 1:length(multi_sine16pi50)){
  mgc.pval[i] = print.stat.optimal(multi_sine16pi50[[i]][[1]], 4)$pvalue
  mcorr.pval[i] = print.stat.optimal(multi_sine16pi50[[i]][[2]], 4)$pvalue
  hhg.pval[i] = print.stat.optimal(multi_sine16pi50[[i]][[3]], 4)$pvalue
  fh.pval[i] = multi_sine16pi50[[i]][[4]]
}
sine16pi.mgc = mean(mgc.pval <= 0.05); sine16pi.mcorr = mean(mcorr.pval <= 0.05)
sine16pi.hhg = mean(hhg.pval <= 0.05); sine16pi.fh = mean(fh.pval <= 0.05)

mgc.pval = c(); mcorr.pval = c(); hhg.pval = c(); fh.pval = c()
for(i in 1:length(multi_square50)){
  mgc.pval[i] = print.stat.optimal(multi_square50[[i]][[1]], 4)$pvalue
  mcorr.pval[i] = print.stat.optimal(multi_square50[[i]][[2]], 4)$pvalue
  hhg.pval[i] = print.stat.optimal(multi_square50[[i]][[3]], 4)$pvalue
  fh.pval[i] = multi_square50[[i]][[4]]
}
square.mgc = mean(mgc.pval <= 0.05); square.mcorr = mean(mcorr.pval <= 0.05)
square.hhg = mean(hhg.pval <= 0.05); square.fh = mean(fh.pval <= 0.05)


mgc.pval = c(); mcorr.pval = c(); hhg.pval = c(); fh.pval = c()
for(i in 1:length(multi_twopara50)){
  mgc.pval[i] = print.stat.optimal(multi_twopara50[[i]][[1]], 4)$pvalue
  mcorr.pval[i] = print.stat.optimal(multi_twopara50[[i]][[2]], 4)$pvalue
  hhg.pval[i] = print.stat.optimal(multi_twopara50[[i]][[3]], 4)$pvalue
  fh.pval[i] = multi_twopara50[[i]][[4]]
}
twopara.mgc = mean(mgc.pval <= 0.05); twopara.mcorr = mean(mcorr.pval <= 0.05)
twopara.hhg = mean(hhg.pval <= 0.05); twopara.fh = mean(fh.pval <= 0.05)


mgc.pval = c(); mcorr.pval = c(); hhg.pval = c(); fh.pval = c()
for(i in 1:length(multi_circle50)){
  mgc.pval[i] = print.stat.optimal(multi_circle50[[i]][[1]], 4)$pvalue
  mcorr.pval[i] = print.stat.optimal(multi_circle50[[i]][[2]], 4)$pvalue
  hhg.pval[i] = print.stat.optimal(multi_circle50[[i]][[3]], 4)$pvalue
  fh.pval[i] = multi_circle50[[i]][[4]]
}
circle.mgc = mean(mgc.pval <= 0.05); circle.mcorr = mean(mcorr.pval <= 0.05)
circle.hhg = mean(hhg.pval <= 0.05); circle.fh = mean(fh.pval <= 0.05)


mgc.pval = c(); mcorr.pval = c(); hhg.pval = c(); fh.pval = c()
for(i in 1:length(multi_ellipse50)){
  mgc.pval[i] = print.stat.optimal(multi_ellipse50[[i]][[1]], 4)$pvalue
  mcorr.pval[i] = print.stat.optimal(multi_ellipse50[[i]][[2]], 4)$pvalue
  hhg.pval[i] = print.stat.optimal(multi_ellipse50[[i]][[3]], 4)$pvalue
  fh.pval[i] = multi_ellipse50[[i]][[4]]
}
ellipse.mgc = mean(mgc.pval <= 0.05); ellipse.mcorr = mean(mcorr.pval <= 0.05)
ellipse.hhg = mean(hhg.pval <= 0.05); ellipse.fh = mean(fh.pval <= 0.05)

mgc.pval = c(); mcorr.pval = c(); hhg.pval = c(); fh.pval = c()
for(i in 1:length(multi_diamond50)){
  mgc.pval[i] = print.stat.optimal(multi_diamond50[[i]][[1]], 4)$pvalue
  mcorr.pval[i] = print.stat.optimal(multi_diamond50[[i]][[2]], 4)$pvalue
  hhg.pval[i] = print.stat.optimal(multi_diamond50[[i]][[3]], 4)$pvalue
  fh.pval[i] = multi_diamond50[[i]][[4]]
}
diamond.mgc = mean(mgc.pval <= 0.05); diamond.mcorr = mean(mcorr.pval <= 0.05)
diamond.hhg = mean(hhg.pval <= 0.05); diamond.fh = mean(fh.pval <= 0.05)

mgc.pval = c(); mcorr.pval = c(); hhg.pval = c(); fh.pval = c()
for(i in 1:length(multi_multi50)){
  mgc.pval[i] = print.stat.optimal(multi_multi50[[i]][[1]], 4)$pvalue
  mcorr.pval[i] = print.stat.optimal(multi_multi50[[i]][[2]], 4)$pvalue
  hhg.pval[i] = print.stat.optimal(multi_multi50[[i]][[3]], 4)$pvalue
  fh.pval[i] = multi_multi50[[i]][[4]]
}
multi.mgc = mean(mgc.pval <= 0.05); multi.mcorr = mean(mcorr.pval <= 0.05)
multi.hhg = mean(hhg.pval <= 0.05); multi.fh = mean(fh.pval <= 0.05)


mgc.pval = c(); mcorr.pval = c(); hhg.pval = c(); fh.pval = c()
for(i in 1:length(multi_ind50)){
  mgc.pval[i] = print.stat.optimal(multi_ind50[[i]][[1]], 4)$pvalue
  mcorr.pval[i] = print.stat.optimal(multi_ind50[[i]][[2]], 4)$pvalue
  hhg.pval[i] = print.stat.optimal(multi_ind50[[i]][[3]], 4)$pvalue
  fh.pval[i] = multi_ind50[[i]][[4]]
}
ind.mgc = mean(mgc.pval <= 0.05); ind.mcorr = mean(mcorr.pval <= 0.05)
ind.hhg = mean(hhg.pval <= 0.05); ind.fh = mean(fh.pval <= 0.05)

###
linear50.power = c(linear.mgc, linear.mcorr, linear.hhg, linear.fh)
exp50.power = c(exp.mgc, exp.mcorr, exp.hhg, exp.fh)
cubic50.power = c(cubic.mgc, cubic.mcorr, cubic.hhg, cubic.fh)
jointN50.power = c(jointN.mgc, jointN.mcorr, jointN.hhg, jointN.fh)
step50.power = c(step.mgc, step.mcorr, step.hhg, step.fh)
quad50.power = c(quad.mgc, quad.mcorr, quad.hhg, quad.fh)
W50.power = c(W.mgc, W.mcorr, W.hhg, W.fh)
spiral50.power = c(spiral.mgc, spiral.mcorr, spiral.hhg, spiral.fh)
bern50.power = c(bern.mgc, bern.mcorr, bern.hhg, bern.fh)
log50.power = c(log.mgc, log.mcorr, log.hhg, log.fh)
root50.power = c(root.mgc, root.mcorr, root.hhg, root.fh)
sine4pi50.power = c(sine4pi.mgc, sine4pi.mcorr, sine4pi.hhg, sine4pi.fh)
sine16pi50.power = c(sine16pi.mgc, sine16pi.mcorr, sine16pi.hhg, sine16pi.fh)
square50.power = c(square.mgc, square.mcorr, square.hhg, square.fh)
twopara50.power = c(twopara.mgc, twopara.mcorr, twopara.hhg, twopara.fh)
circle50.power = c(circle.mgc, circle.mcorr, circle.hhg, circle.fh)
ellipse50.power = c(ellipse.mgc, ellipse.mcorr, ellipse.hhg, ellipse.fh)
diamond50.power = c(diamond.mgc, diamond.mcorr, diamond.hhg, diamond.fh)
multi50.power = c(multi.mgc, multi.mcorr, multi.hhg, multi.fh)
ind50.power = c(ind.mgc, ind.mcorr, ind.hhg, ind.fh)

RDPG.mat = matrix(0, nrow = 20, ncol = 4)
RDPG.mat = rbind(linear50.power, exp50.power, cubic50.power, jointN50.power,
                 step50.power, quad50.power, W50.power, spiral50.power,
                 bern50.power, log50.power, root50.power, sine4pi50.power,
                 sine16pi50.power, square50.power, twopara50.power, circle50.power,
                 ellipse50.power, diamond50.power, multi50.power, ind50.power)
colnames(RDPG.mat) = c("MGC", "dCorr", "HHG", "FH")
rownames(RDPG.mat) = c("Linear", "Exponential", "Cubic", "Joint Normal",
                       "Step", "Quadratic", "W Shape", "Spiral",
                       "Bernoulli", "Logarithm", "Fourth Root", "Sine (4pi)",
                       "Sine (16pi)", "Square", "Two Parabolas", "Circle",
                       "Ellipse", "Diamond", "Multiplicative", "Independent")

pdf("../Figure/multi_RDPG_optimal.pdf", width = 35, height = 20)
par(mfrow = c(2,1), cex.lab = 4, cex.axis = 3,
    mar = c(10,10,3,26), tcl = 0.5, xpd=TRUE)
plot(c(1:10), rep(-1,10), ylim = c(0,1), pch = NULL,
     main = "", ylab = "Power", 
     xaxt = 'n', yaxt = 'n', mgp = c(4,2,0), xlab = "")
points(rep(1,4), RDPG.mat[1,], type = "p",
       pch = c(1,4,5,6), col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"),
       lwd = 10, cex = 4)
points(rep(2,4), RDPG.mat[2,], type = "p",
       pch = c(1,4,5,6), col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"),
       lwd = 10, cex = 4)
points(rep(3,4), RDPG.mat[3,], type = "p",
       pch = c(1,4,5,6), col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"),
       lwd = 10, cex = 4)
points(rep(4,4), RDPG.mat[4,], type = "p",
       pch = c(1,4,5,6), col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"),
       lwd = 10, cex = 4)
points(rep(5,4), RDPG.mat[5,], type = "p",
       pch = c(1,4,5,6), col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"),
       lwd = 10, cex = 4)
points(rep(6,4), RDPG.mat[6,], type = "p",
       pch = c(1,4,5,6), col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"),
       lwd = 10, cex = 4)
points(rep(7,4), RDPG.mat[7,], type = "p",
       pch = c(1,4,5,6), col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"),
       lwd = 10, cex = 4)
points(rep(8,4), RDPG.mat[8,], type = "p",
       pch = c(1,4,5,6), col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"),
       lwd = 10, cex = 4)
points(rep(9,4), RDPG.mat[9,], type = "p",
       pch = c(1,4,5,6), col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"),
       lwd = 10, cex = 4)
points(rep(10,4), RDPG.mat[10,], type = "p",
       pch = c(1,4,5,6), col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"),
       lwd = 10, cex = 4)
axis(1, at=1:10,  labels = FALSE, tck = 0.05)
text(c(1-0.2, 2-0.4, 3-0.2, 4-0.5, 5-0.1, 6-0.3, 7-0.3, 8-0.2, 9-0.2, 10-0.2), rep(-0.2,10), labels = c("Linear", "Exponential", "Cubic", "Joint Normal",
                                                                                                        "Step", "Quadratic", "W Shape", "Spiral",
                                                                                                        "Bernoulli", "Logarithm"), cex= 3, pos=4, srt = 0, xpd=TRUE)
axis(2, at=c(0, 0.25, 0.5, 0.75, 1), labels=c(0, 0.25, 0.5, 0.75, 1), cex = 2,  tck = 0.05)

plot(c(1:10), rep(-1,10), ylim = c(0,1), pch = NULL,
     main = "", ylab = "Power", xlab = "",
     xaxt = 'n', yaxt = 'n', mgp = c(4,2,0))
points(rep(1,4), RDPG.mat[11,], type = "p",
       pch = c(1,4,5,6), col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"),
       lwd = 10, cex = 4)
points(rep(2,4), RDPG.mat[12,], type = "p",
       pch = c(1,4,5,6), col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"),
       lwd = 10, cex = 4)
points(rep(3,4), RDPG.mat[13,], type = "p",
       pch = c(1,4,5,6), col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"),
       lwd = 10, cex = 4)
points(rep(4,4), RDPG.mat[14,], type = "p",
       pch = c(1,4,5,6), col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"),
       lwd = 10, cex = 4)
points(rep(5,4), RDPG.mat[15,], type = "p",
       pch = c(1,4,5,6), col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"),
       lwd = 10, cex = 4)
points(rep(6,4), RDPG.mat[16,], type = "p",
       pch = c(1,4,5,6), col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"),
       lwd = 10, cex = 4)
points(rep(7,4), RDPG.mat[17,], type = "p",
       pch = c(1,4,5,6), col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"),
       lwd = 10, cex = 4)
points(rep(8,4), RDPG.mat[18,], type = "p",
       pch = c(1,4,5,6), col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"),
       lwd = 10, cex = 4)
points(rep(9,4), RDPG.mat[19,], type = "p",
       pch = c(1,4,5,6), col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"),
       lwd = 10, cex = 4)
points(rep(10,4), RDPG.mat[20,], type = "p",
       pch = c(1,4,5,6), col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"),
       lwd = 10, cex = 4)
text(c(1-0.5, c(2:4)-0.3, 5-0.5, 6-0.2, 7-0.2, 8-0.2, 9-0.2, 10), rep(-0.2,10), c("Fourth Root", "Sine (4pi)",
                                                                                  "Sine (16pi)", "Square", "Two Parabolas", "Circle",
                                                                                  "Ellipse", "Diamond", "Multiplicative", "Indep"), cex= 3, pos=4, srt = 0, xpd=TRUE)
axis(1, at=c(1:10), labels=FALSE, cex = 2, tck = 0.05)
axis(2, at=c(0, 0.25, 0.5, 0.75, 1), labels=c(0, 0.25, 0.5, 0.75, 1), cex = 2,  tck = 0.05)

legend("topright", inset=c(-0.18, 0), 
       c( expression(MGC %.% DM), expression(dCorr %.% DM), 
          expression(HHG %.% DM), "FH Test"), 
       col = c("red",  "dodgerblue", "lightsalmon4", "darkgreen"),
       pch = c(1,4,5,6), bty = 'n', lwd = 10, xpd = NA, cex = 3)
dev.off()

