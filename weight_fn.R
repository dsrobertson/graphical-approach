################################################################################
# Function: graphWeight                                                           
# Implements the graphical weighting strategy of Bretz et al. (2011)               
################################################################################
# Inputs                                                                       
# ######                                                                      
#       
#       J - subset of hypotheses to calculated the weights for
#       graph - list of transition matrix (G) and weights
#       Iset - initial set of hypotheses
#       
################################################################################
# Outputs                                                                      
# #######                                                                      
#                                                                              
# List of weights                                           
################################################################################

graphWeight = function(J, graph, Iset = 1:ncol(graph$G)) {
  
  G = graph$G
  w = graph$weights
  n = ncol(G)
  
  Gtemp = matrix(0, ncol = n, nrow = n)
  
  Jc = setdiff(Iset, J)
  
  while(length(Jc) >= 1){
    
    j = Jc[1]
    
    for(i in Iset){
      
      w[i] = w[i] + w[j]*G[j,i]
      
      GijGji = G[i,j]*G[j,i]
      
      if(GijGji<1){
        Gij = G[i,j]
        Gtemp[i,] = (G[i,] + Gij*G[j,])/(1 - GijGji)
      }
      
      Gtemp[i,i] = 0
    }
    
    G = Gtemp
    Jc = Jc[-1]
    Iset = Iset[-which(Iset==j)]
  }
  
  return(w)
}
