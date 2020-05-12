### Code for Section 6.1: an analysis of the Pharmacodynamic study reported by
### Ferber et al. (2011)
################################################################################

# Source required functions
source('weight_fn.R')
source('kFWER.R')
source('FDP.R')

# p-values as reported by Ferber et al.
pvalues = c(0.7808, 0.9433, 0.9993,
            0.06, 0.0053, 1e-05,
            0.0137, 6.5e-06, 1.7e-11,
            0.0724, 2.8e-6, 9.1e-8,
            0.0162, 9.1e-8, 8.1e-13)

# Transition matrix
G = rbind( T1D1=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
           T1D2=c(0.5, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
           T1D3=c(0, 0.5, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
           T2D1=c(0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
           T2D2=c(0, 0, 0.5, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0),
           T2D3=c(0, 0, 0.5, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0),
           T3D1=c(0, 0, 0, 0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
           T3D2=c(0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0.5, 0, 0, 0, 0, 0),
           T3D3=c(0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.5, 0, 0, 0, 0),
           T4D1=c(0, 0, 0, 0, 0, 0, 0.5, 0.5, 0, 0, 0, 0, 0, 0, 0),
           T4D2=c(0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0.5, 0, 0),
           T4D3=c(0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.5, 0),
           T5D1=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.5, 0, 0, 0, 0),
           T5D2=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.5, 0, 0),
           T5D3=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0.5, 0))


################################################################################
### Table 2: Rejected hypotheses for the pharmacodynamic study of Ferber et al.,
### with initial weights of 1/3 for T4D3; T5D2 and T5D3.


# Intial weights
w = c(rep(0,11), 1/3, 0, 1/3, 1/3)

# Set up graph
graph = list(G=G, weights=w)


### k-FWER control

# Generalised
names(which(kgMCP(k=1, graph, pvalues)$rejected == 1))
names(which(kgMCP(k=2, graph, pvalues)$rejected == 1))
names(which(kgMCP(k=3, graph, pvalues)$rejected == 1))

# Augmented
names(which(kgMCPaug(k=1, graph, pvalues)$rejected == 1))
names(which(kgMCPaug(k=2, graph, pvalues)$rejected == 1))
names(which(kgMCPaug(k=3, graph, pvalues)$rejected == 1))

### FDP control

# Generalised
names(which(FDP_gMCP(graph, pvalues, gamma.fdp = 0.1)$rejected == 1))
names(which(FDP_gMCP(graph, pvalues, gamma.fdp = 0.2)$rejected == 1))
names(which(FDP_gMCP(graph, pvalues, gamma.fdp = 0.3)$rejected == 1))

# Augmented
names(which(FDP_gMCPaug(graph, pvalues, gamma.fdp = 0.1)$rejected == 1))
names(which(FDP_gMCPaug(graph, pvalues, gamma.fdp = 0.2)$rejected == 1))
names(which(FDP_gMCPaug(graph, pvalues, gamma.fdp = 0.3)$rejected == 1))


################################################################################
### Table 3: Rejected hypotheses for the pharmacodynamic study of Ferber et al.
### with initial weights of 1/15 for each hypothesis.

# Intial weights
w = rep(1/15, 15)

# Set up graph
graph = list(G=G, weights=w)


### k-FWER control

# Generalised
names(which(kgMCP(k=1, graph, pvalues)$rejected == 1))
names(which(kgMCP(k=2, graph, pvalues)$rejected == 1))
names(which(kgMCP(k=3, graph, pvalues)$rejected == 1))

# Augmented
names(which(kgMCPaug(k=1, graph, pvalues)$rejected == 1))
names(which(kgMCPaug(k=2, graph, pvalues)$rejected == 1))
names(which(kgMCPaug(k=3, graph, pvalues)$rejected == 1))

### FDP control

# Generalised
names(which(FDP_gMCP(graph, pvalues, gamma.fdp = 0.1)$rejected == 1))
names(which(FDP_gMCP(graph, pvalues, gamma.fdp = 0.2)$rejected == 1))
names(which(FDP_gMCP(graph, pvalues, gamma.fdp = 0.3)$rejected == 1))

# Augmented
names(which(FDP_gMCPaug(graph, pvalues, gamma.fdp = 0.1)$rejected == 1))
names(which(FDP_gMCPaug(graph, pvalues, gamma.fdp = 0.2)$rejected == 1))
names(which(FDP_gMCPaug(graph, pvalues, gamma.fdp = 0.3)$rejected == 1))
