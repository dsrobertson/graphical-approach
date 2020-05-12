################################################################################
# Function: kgMCPaug                                                             
# Implements the augmented graphical approach for k-FWER control             
################################################################################
# Inputs                                                                       
# ######                                                                      
#                                                                              
#       k - number of false rejections for k-FWER control                     
#       graph - list of transition matrix (G) and weights
#       pvalues - vector of p-values
#       alpha - desired k-FWER, defaults to 0.05
#       delta - determines how many of the 'free' rejections are used
#       
################################################################################
# Outputs                                                                      
# #######                                                                      
#                                                                              
# List of p-values, alpha and the rejected hypotheses                                           
################################################################################

kgMCPaug = function(k, graph, pvalues, alpha = 0.05, delta = 1) {
  
  G = graph$G
  n = ncol(G)
  h = rep(0,n)
  
  names(h) = rownames(G)
  
  a = alpha*graph$weights
  
  crit = 0
  
  while(crit == 0){
    test = (pvalues <= a)
    
    if(any(test)){
      
      rej = which.max(test)
      h[rej] = 1
      Gtemp = matrix(0, ncol = n, nrow = n)
      
      for(i in 1:n){
        
        a[i] = a[i] + a[rej]*G[rej,i]
        GG = G[i,rej]*G[rej,i]
        
        if(G[i,rej]*G[rej,i]<1){
          
          Gi = G[i,rej]
          Gtemp[i,] = (G[i,] + Gi*G[rej,])/(1 - GG)
        }
      }
      G = Gtemp
      G[rej,] = G[,rej] = 0
      a[rej] = 0
      
    } else {
      crit = 1
    }
  }
  
  Iset = which(h==0)
  R.aug = integer(0)
  
  if(k > 1 && length(Iset)>=1){
    
    w = graphWeight(Iset, graph)
    
    while(length(R.aug) < k-1 && length(Iset) >= 1 && max(w) > 0) {
      
      j = Iset[which.min(pvalues[Iset]/w[Iset])]
      
      rej = (pvalues[j] <= k*w[j]*delta)
      
      if(rej){
        
        h[j] = 1
        R.aug = union(R.aug, j)
        Iset = setdiff(Iset, j)
        
        w = graphWeight(Iset, graph)
        
      } else {
        w = 0
      }
    }
  }
  
  return(list(pvalues = pvalues, alpha = alpha, rejected = h))
  
}


################################################################################
# Function: kgMCP                                                              
# Implements the generalised graphical approach for k-FWER control             
################################################################################
# Inputs                                                                       
# ######                                                                      
#                                                                              
#       k - number of false rejections for k-FWER control                     
#       graph - list of transition matrix (G) and weights
#       pvalues - vector of p-values
#       alpha - desired k-FWER, defaults to 0.05
#       delta - determines how many of the 'free' rejections are used
#       
################################################################################
# Outputs                                                                      
# #######                                                                      
#                                                                              
# List of p-values, alpha and the rejected hypotheses                                           
################################################################################

kgMCP = function(k, graph, pvalues, alpha = 0.05, delta = 1) {
  
  G = graph$G
  w = graph$weights
  n = ncol(G)
  
  Iset = 1:n
  
  rej = (pvalues <= w*k*alpha)
  R = which(rej)
  Iset = setdiff(Iset, R)
  
  if(sum(rej)>=k){
    
    while(length(Iset) >= 1){
      
      J = combn(R, k-1)
      
      K = rbind(J, matrix(Iset, nrow = length(Iset), ncol = dim(J)[2]))
      
      wset = apply(K, 2, graphWeight, graph=graph)
      
      rej[Iset] = (pvalues[Iset] <= apply(wset[Iset,,drop=F], 1, min)*k*alpha)
      
      if(sum(rej[Iset])==0){
        Iset = integer(0)
      } else {
        R = union(R, which(rej))
        Iset = setdiff(Iset, which(rej))
      }
    }
  } else if(sum(rej) < k-1) {
    
    w = graphWeight(Iset, graph)
    
    while(sum(rej) < k-1 && length(Iset) >= 1 && max(w) > 0) {
      
      j = Iset[which.min(pvalues[Iset]/w[Iset])]
      
      rej[j] = (pvalues[j] <= k*w[j]*delta)
      
      if(rej[j]){
        
        R = union(R, j)
        Iset = setdiff(Iset, j)
        
        w = graphWeight(Iset, graph)
        
      } else {
        w = 0
      }
    }
  }
  
  h = rep(0,n)
  # names(h) = paste0('H',1:n)
  names(h) = rownames(G)
  h[R] = 1
  
  return(list(pvalues = pvalues, alpha = alpha, rejected = h))
}