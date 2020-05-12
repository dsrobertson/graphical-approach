###############################################################################
#### Gatekeeping: reject at least 6 out of 9 hypotheses
#### With a single secondary hypothesis
###############################################################################


################################################################################
# Function: gate_weight                                                           
# Calculates the weights wehn using the gatekeeping strategy               
################################################################################
# Inputs                                                                       
# ######                                                                      
#       
#       J - subset of hypotheses to calculate the weights for
#       
################################################################################
# Outputs                                                                      
# #######                                                                      
#                                                                              
# Vector of weights                                           
################################################################################

gate_weight = function(J) {
  
  w = rep(0, 10)
  
  J.primary = setdiff(J, 10)
  n.primary = length(J.primary)
  
  secondary = (10 %in% J)
  
  if(n.primary > 3){
    
    w[J.primary] = rep(1/n.primary, n.primary)
    
  } else if(n.primary == 3){
    
    w[J.primary] = rep(83/252, 3)
    
    if(secondary){
      w[10] = 1/84 
    }
  } else if(n.primary == 2){
    
    w[J.primary] = rep(11/24, 2)
    
    if(secondary){
      w[10] = 1/12
    }
  } else if(n.primary == 1){
    
    w[J.primary] = 2/3
    
    if(secondary){
      w[10] = 1/3
    }
  } else {
    w[10] = 1
  }
  return(w)
}

################################################################################
# Function: gate_kHolmAug                                                           
# Implements the augmented k-FWER Holm procedure for the gatekeeping strategy              
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

gate_kHolmAug = function(k, pvalues, alpha = 0.05, delta = 1) {
  
  if(length(pvalues) != 10){
    stop('There should be 10 p-values.')
  }
  
  h = rep(0,10)
  names(h) = paste0('H',1:10)
  
  Iset = 1:10
  
  w = c(rep(1/9, 9), 0)
  
  rej = (pvalues <= w*k*alpha)
  R = which(rej)
  Iset = setdiff(Iset, R)
  
  while(length(Iset) >= 1){
    
    wset = gate_weight(Iset)
    
    rej[Iset] = (pvalues[Iset] <= wset[Iset]*k*alpha)
    
    if(sum(rej[Iset])==0){
      Iset = integer(0)
    } else {
      R = union(R, which(rej))
      Iset = setdiff(Iset, which(rej))
    }
  }
  
  h[R] = 1
  
  Iset = which(!rej)
  R.aug = integer(0)
  
  if(k > 1 && length(Iset)>=1){
    
    w = gate_weight(Iset)
    
    while(length(R.aug) < k-1 && length(Iset) >= 1 && max(w) > 0) {
      
      j = Iset[which.min(pvalues[Iset]/w[Iset])]
      
      rej = (pvalues[j] <= k*w[j]*delta)
      
      if(rej){
        
        h[j] = 1
        R.aug = union(R.aug, j)
        Iset = setdiff(Iset, j)
        
        w = gate_weight(Iset)
        
      } else {
        w = 0
      }
    }
  }
  
  return(list(pvalues = pvalues, alpha = alpha, rejected = h))
}


################################################################################
# Function: gate_kHolm                                                           
# Implements the generalised k-FWER Holm procedure for the gatekeeping strategy              
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

gate_kHolm = function(k, pvalues, alpha = 0.05, delta = 1) {
  
  if(length(pvalues) != 10){
    stop('There should be 10 p-values.')
  }
  
  Iset = 1:10
  
  w = c(rep(1/9, 9), 0)
  
  rej = (pvalues <= w*k*alpha)
  R = which(rej)
  Iset = setdiff(Iset, R)
  
  if(sum(rej)>=k){
    
    while(length(Iset) >= 1){
      
      J = combn(R, k-1)
      
      K = rbind(J, matrix(Iset, nrow = length(Iset), ncol = dim(J)[2]))
      
      wset = apply(K, 2, gate_weight)
      
      rej[Iset] = (pvalues[Iset] <= apply(wset[Iset,,drop=F], 1, min)*k*alpha)
      
      if(sum(rej[Iset])==0){
        Iset = numeric(0)
      } else {
        R = union(R, which(rej))
        Iset = setdiff(Iset, which(rej))
      }
    }
  } else if(sum(rej) < k-1) {
    
    w = gate_weight(Iset)
    
    while(sum(rej) < k-1 && length(Iset) >= 1 && max(w) > 0) {
      
      j = Iset[which.min(pvalues[Iset]/w[Iset])]
      
      rej[j] = (pvalues[j] <= k*w[j]*delta)
      
      if(rej[j]){
        
        R = union(R, j)
        Iset = setdiff(Iset, j)
        
        w = gate_weight(Iset)
        
      } else {
        w = 0
      }
    }
  }
  
  h = rep(0,10)
  names(h) = paste0('H',1:10)
  h[R] = 1
  
  return(list(pvalues = pvalues, alpha = alpha, rejected = h))
}

################################################################################
# Function: FDP_Holm                                                           
# Implements the augmented FDP Holm procedure for the gatekeeping strategy             
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

gate_FDP_HolmAug = function(pvalues, gamma.fdp = 0.1, alpha = 0.05,
                            delta=1) {
  
  if(length(pvalues) != 10){
    stop('There should be 10 p-values.')
  }
  
  h = w = rep(0,10)
  names(h) = paste0('H',1:10)
  
  h = rep(0,10)
  names(h) = paste0('H',1:10)
  
  Iset = 1:10
  
  w = c(rep(1/9, 9), 0)
  
  rej = (pvalues <= w*alpha)
  R = which(rej)
  Iset = setdiff(Iset, R)
  
  while(length(Iset) >= 1){
    
    wset = gate_weight(Iset)
    
    rej[Iset] = (pvalues[Iset] <= wset[Iset]*alpha)
    
    if(sum(rej[Iset])==0){
      Iset = integer(0)
    } else {
      R = union(R, which(rej))
      Iset = setdiff(Iset, which(rej))
    }
  }
  
  h[R] = 1
  
  Iset = which(!rej)
  R.aug = integer(0)
  
  D = floor(sum(h)*gamma.fdp/(1-gamma.fdp))
  
  if(D >= 1 && length(Iset)>=1){
    
    w = gate_weight(Iset)
    
    while(length(R.aug) < D && length(Iset) >= 1 && max(w) > 0) {
      
      j = Iset[which.min(pvalues[Iset]/w[Iset])]
      
      rej = (pvalues[j] <= w[j]*delta)
      
      if(rej){
        
        h[j] = 1
        R.aug = union(R.aug, j)
        Iset = setdiff(Iset, j)
        
        w = gate_weight(Iset)
        
      } else {
        w = 0
      }
    }
  }
  return(list(pvalues = pvalues, alpha = alpha, rejected = h))
}


################################################################################
# Function: FDP_Holm                                                           
# Implements the generalised FDP Holm procedure for the gatekeeping strategy             
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

gate_FDP_Holm = function(pvalues, gamma.fdp = 0.1, alpha = 0.05,
                         delta=1) {
  
  if(length(pvalues) != 10){
    stop('There should be 10 p-values.')
  }
  
  n = length(pvalues)
  
  for(k in 1:n){
    
    tmpk = gate_kHolm(k, pvalues, alpha, delta)
    
    if(sum(tmpk$rejected) < (k/gamma.fdp-1)){
      return(tmpk)
    }
  }
  return(tmpk)
}


################################################################################
# Function: denomT                                                           
# Estimate the variance of difference between means. Based on code from the
# rPowerSampleSize package
################################################################################
# Inputs                                                                       
# ######                                                                      
#                                                                              
#       XE - endpoints for experimental treatment
#       XC - endpoints for control treatment
#       rho - correlation coefficient
#       
################################################################################
# Outputs                                                                      
# #######                                                                      
#                                                                              
# Estimate of variance                                           
################################################################################

denomT = function (XE, XC, rho = 0) 
{
  m = ncol(XE)
  nE = nrow(XE)
  nC = nrow(XC)
  XEbar = colMeans(XE)
  XCbar = colMeans(XC)
  
  VartildeE = matrix(0, nrow = m, ncol = m)
  for (i in 1:nE) {
    VartildeE = VartildeE + (XE[i, ] - XEbar) %*% 
      t(XE[i, ] - XEbar)
  }
  VartildeC = matrix(0, nrow = m, ncol = m)
  for (i in 1:nC) {
    VartildeC = VartildeC + (XC[i, ] - XCbar) %*% 
      t(XC[i, ] - XCbar)
  }
  
  VartildeE = VartildeE/nE
  VarhatE = nE * VartildeE/(nE - 1)
  VartildeC = VartildeC/nC
  VarhatC = nC * VartildeC/(nC - 1)
  res = diag(VarhatE)/nE + diag(VarhatC)/nC
  
  return(res)
}
