NetworkTest.diffusion.stat = function(G, X, option, diffusion, t.range, n.perm){
  
  A = as.matrix(get.adjacency(G))
  D = diag(pmin( (rowSums(A))^(-1/2) , 1))
  P = D %*% A %*% D # a nomalized graph Laplacian
  
  if(diffusion == FALSE){ 
    Dx = as.matrix(dist((A)), diag = TRUE, upper = TRUE) 
    Dy = as.matrix(dist((X)), diag = TRUE, upper = TRUE)
    
    if(option == 1){
      pvalues = mgc.test(Dx, Dy, rep = n.perm, option = 'mgc')$pMGC
      return(pvalues)
    }else if(option == 2){
      pvalues = dcov.test(Dx, Dy, index = 1.0, R = n.perm)$p.value
      return(pvalues)
    }else if(option == 3){
      pvalues = hhg.test(Dx, Dy, nr.perm = n.perm)$perm.pval.hhg.sl
      return(pvalues)
    }
  }
  
  mgc.stat = c(); dcor.stat = c(); hhg.stat = c()
  
  for(s in 1:length(t.range)){

    # as a dimensional choice, use the second elbow of absolute eivengalues from a diffusion map at t = 1.
    diffusion.q  =  min(max(getElbows(abs(print.lambda(P, times = 1)[[1]]), plot = FALSE, n = 2)), nrow(A)-1)
    U  =  dmap.q(P, t.range[s], diffusion.q)[[1]] # dmap.q(Laplacian, Markov time t, dimension q)
    Dx = as.matrix(dist((U)), diag = TRUE, upper = TRUE) 
    Dy = as.matrix(dist((X)), diag = TRUE, upper = TRUE)
    
    if(option == 1){
      mgc = mgc.sample(Dx, Dy, option = 'mgc')[[1]]
      mgc.stat[s] = mgc
    }
    
    if(option == 2){
      dcor = dcor.ttest(Dx, Dy, distance = TRUE)
      dcor.stat[s] = dcor$statistic
    }
      
    if(option == 3){
      hhg = hhg.test(Dx, Dy, nr.perm = n.perm)
      hhg.stat[s] = hhg$sum.lr
    }
  }
  
  tmp = matrix(0, nrow = n.perm, ncol = length(t.range))
  
  for(r in 1:n.perm){
    
    for(s in 1:length(t.range)){
      
      diffusion.q  =  min( max(getElbows(abs(print.lambda(P, times = 1)[[1]]), plot = FALSE, n = 2)), nrow(A)-1)
      U  =  dmap.q(P, t.range[s], diffusion.q)[[1]]

      if(class(X) == "numeric"){
        per = sample(length(X));
        newX = X[per]
      }else if(class(X) == "matrix"){
        per=sample(nrow(X));
        newX = X[per,]
      }


      Dx = as.matrix(dist(U), diag = TRUE, upper = TRUE) 
      Dy = as.matrix(dist(newX), diag = TRUE, upper = TRUE)
    
      if(option == 1){
        mgc = mgc.sample(Dx, Dy, option = 'mgc')[[1]]
        tmp[r,s] = mgc
      }
      
      if(option == 2){
        dcor = dcor.ttest(Dx, Dy, distance = TRUE)
        tmp[r,s] = dcor$statistic
      }
      
      if(option == 3){
        hhg = hhg.test(Dx, Dy, nr.perm = n.perm)
        tmp[r,s] = hhg$sum.lr
      }
    } 
  }
  
  if(option == 1) return(list(mgc.stat, tmp))
  if(option == 2) return(list(dcor.stat, tmp))
  if(option == 3) return(list(hhg.stat, tmp))
}






TwoGraphs.diffusion.stat = function(G1, G2, option, diffusion, t.range, n.perm){

    A1 = as.matrix(get.adjacency(G1))
    A2 = as.matrix(get.adjacency(G2))

    D1 = diag(pmin( (rowSums(A1))^(-1/2) , 1))
    P1 = D1 %*% A1 %*% D1

    D2 = diag(pmin( (rowSums(A2))^(-1/2) , 1))
    P2 = D2 %*% A2 %*% D2
  
  if(diffusion == FALSE){ 
      Dx = as.matrix(dist((A1)), diag = TRUE, upper = TRUE) 
      Dy = as.matrix(dist((A2)), diag = TRUE, upper = TRUE)
    
    if(option == 1){
      result = MGCPermutationTest(Dx, Dy, n.perm, option = 'mcor')
      return(result)
    }else if(option == 2){
      pvalues = dcov.test(Dx, Dy, index = 1.0, R = n.perm)$p.value
      return(pvalues)
    }else if(option == 3){
      pvalues = hhg.test(Dx, Dy, nr.perm = n.perm)$perm.pval.hhg.sl
      return(pvalues)
    }
  }

  mgc.stat = c(); dcor.stat = c(); hhg.stat = c()
  
  for(s in 1:length(t.range)){

    diffusion.q1  =  min(max(getElbows(abs(print.lambda(P1, times = 1)[[1]]), plot = FALSE, n = 2)), nrow(A1) - 1)
    diffusion.q2  =  min(max(getElbows(abs(print.lambda(P2, times = 1)[[1]]), plot = FALSE, n = 2)), nrow(A2) - 1)

    dmap1 = dmap.q(P1, times = t.range[s], q = diffusion.q1)[[1]]
    dmap2 = dmap.q(P2, times = t.range[s], q = diffusion.q2)[[1]] 
    
    Dx = as.matrix(dist(dmap1), diag = TRUE, upper = TRUE)
    Dy = as.matrix(dist(dmap2), diag = TRUE, upper = TRUE)
    
    if(option == 1){
      mgc = MGCSampleStat(Dx, Dy, option = 'mcor')[[1]]
      mgc.stat[s] = mgc
    }
    
    if(option == 2){
      dcor = dcor.ttest(Dx, Dy, distance = TRUE)
      dcor.stat[s] = dcor$statistic
    }
    
    
    if(option == 3){
      hhg = hhg.test(Dx, Dy, nr.perm = n.perm)
      hhg.stat[s] = hhg$sum.lr
    }
  }

  tmp = matrix(0, nrow = n.perm, ncol = length(t.range))

  for(r in 1:n.perm){

    for(s in 1:length(t.range)){

     perm.idx = sample(nrow(A2))
     new.A2 = A2[perm.idx, perm.idx]
     new.P2 = new.A2 / pmax(rowSums(new.A2), 1)

     diffusion.q1  =  min(max(getElbows(abs(print.lambda(P1, times = 1)[[1]]), plot = FALSE, n = 2)), nrow(A1) - 1)
     diffusion.q2  =  min(max(getElbows(abs(print.lambda(new.P2, times = 1)[[1]]), plot = FALSE, n = 2)), nrow(A2) - 1)

     
     dmap1 = dmap.q(P1, times = t.range[s], q = diffusion.q1)[[1]]
     dmap2 = dmap.q(new.P2, times = t.range[s], q = diffusion.q2)[[1]] 
    
     Dx = as.matrix(dist(dmap1), diag = TRUE, upper = TRUE)
     Dy = as.matrix(dist(dmap2), diag = TRUE, upper = TRUE)
    
      if(option == 1){
        mgc = MGCSampleStat(Dx, Dy, option = 'mcor')[[1]]
        tmp[r,s] = mgc
      }
    
      if(option == 2){
        dcor = dcor.ttest(Dx, Dy, distance = TRUE)
        tmp[r,s] = dcor$statistic
      }
    
      if(option == 3){
        hhg = hhg.test(Dx, Dy, nr.perm = n.perm)
        tmp[r,s] = hhg$sum.lr
      }
    } 
  }

  if(option == 1) return(list(mgc.stat, tmp))
  if(option == 2) return(list(dcor.stat, tmp))
  if(option == 3) return(list(hhg.stat, tmp))
}