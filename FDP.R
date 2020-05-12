################################################################################
# Function: FDP_gMCPaug                                                           
# Implements the augments graphical approach for FDP control             
################################################################################
# Inputs                                                                       
# ######                                                                      
#                                                                              
#       graph - list of transition matrix (G) and weights
#       pvalues - vector of p-values
#       gamma.fdp - tail probability bound for the FDP, defaults to 0.1
#       alpha - desired FDP level, defaults to 0.05
#       delta - determines how many of the 'free' rejections are used
#       
################################################################################
# Outputs                                                                      
# #######                                                                      
#                                                                              
# List of p-values, alpha and the rejected hypotheses                                           
################################################################################

FDP_gMCPaug = function(graph, pvalues, gamma.fdp = 0.1, alpha = 0.05,
                       delta=1) {
  
  G = graph$G
  n = ncol(G)
  h = rep(0,n)

  names(h) = rownames(G)
  
  a = alpha*graph$weights
  
  crit = 0
  
  while(crit == 0){
    test = pvalues <= a
    
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
  
  D = floor(sum(h)*gamma.fdp/(1-gamma.fdp))
  
  if(D >= 1 && length(Iset)>=1){
    
    w = graphWeight(Iset, graph)
    
    while(length(R.aug) < D && length(Iset) >= 1 && max(w) > 0) {
      
      j = Iset[which.min(pvalues[Iset]/w[Iset])]
      
      rej = (pvalues[j] <= w[j]*delta)
      
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
# Function: FDP_gMCP                                                              
# Implements the generalised graphical approach for FDP control             
################################################################################
# Inputs                                                                       
# ######                                                                      
#                                                                              
#       graph - list of transition matrix (G) and weights
#       pvalues - vector of p-values
#       gamma.fdp - tail probability bound for the FDP, defaults to 0.1
#       alpha - desired FDP level, defaults to 0.05
#       delta - determines how many of the 'free' rejections are used
#       
################################################################################
# Outputs                                                                      
# #######                                                                      
#                                                                              
# List of p-values, alpha and the rejected hypotheses                                           
################################################################################

FDP_gMCP = function(graph, pvalues, gamma.fdp = 0.1, alpha = 0.05,
                    delta=1) {
  
  n = length(pvalues)
  
  for(k in 1:n){
    
    tmpk = kgMCP(k, graph, pvalues, alpha, delta)
    
    if(sum(tmpk$rejected) < (k/gamma.fdp-1)){
      return(tmpk)
    }
  }
  return(tmpk)
}

