### Code for Section 6.2: simulation study of the Pre-RELAX-AHF trial reported
### by Teerlink et al. (2009). Based on code from the rPowerSampleSize package.
################################################################################

# Source required functions
source('gatekeeping.R')

# Empirical means for control treatment
muC = c(23 / 100, 1679, 1 - 21 / 100, -12.0, 44.2, 1 - 17.2 / 100, 1 -
           14.3 / 100, 13 / 100, 7 / 100)

# Empirical means for experimental treatment
muE = c(40 / 100, 2567, 1 - 12 / 100, -10.2, 47.9, 1 - 2.60 / 100, 1,
         21 / 100, 7 / 100)

# Empirical standard deviations for control treatment
sdC = c(sqrt(0.23 * (1 - 0.23)), 2556, sqrt(0.79 * (1 - 0.79)), 7.3,
         14.2, sqrt(0.828 * (1 - 0.828)), sqrt(0.857 * (1 - 0.857)),
         sqrt(0.13 * (1 - 0.13)), sqrt(0.07 * (1 - 0.07)))

# Empirical standard deviations for experimental treatment
sdE = c(sqrt(0.4 * (1 - 0.4))  , 2898, sqrt(0.88 * (1 - 0.88)), 6.1,
         10.1, sqrt(0.974 * (1 - 0.974)), 1e-12, sqrt(0.21 * (1 - 0.21)), 
         sqrt(0.07 * (1 - 0.07)))

# Correlation coefficient
rho = 0.5

# Correlation matrix
cor = matrix(rho, nrow = 9, ncol = 9)
diag(cor) = 1

# Covariance matrix for experimental treatment
SigmaE = diag(sdE) %*% cor %*% diag(sdE)

# Covariance matrix for control treatment
SigmaC = diag(sdC) %*% cor %*% diag(sdC)

# Mean of secondary outcome
mu.secondary = 3

# Function to calculate trace of a matrix M
rtr = function(M){
  sum(diag(M)) 
}

################################################################################
### Run simulation

alpha = 0.1

# Number of subjects
nE = nC = 200   

# Number of simulation relicates
nsim = 10^4

# Vector of p-values
pvals = rep(NA, 10)

# Output matrices
out.kgMCP1 = out.kgMCP.aug1 = out.kgMCP2 = out.kgMCP.aug2 = 
  out.kgMCP3 = out.kgMCP.aug3 = 
  out.FDP.gMCP1 = out.FDP.gMCP.aug1 = out.FDP.gMCP2 = out.FDP.gMCP.aug2 =
  out.FDP.gMCP3 = out.FDP.gMCP.aug3 =
  matrix(0, nrow = nsim, ncol = 10)

set.seed(3)

for (i in 1:nsim) {
  # Endpoints for experimental treatment
  XE = mvtnorm::rmvnorm(nE, mean = muE, sigma = SigmaE)
  
  # Endpoints for control treatment
  XC = mvtnorm::rmvnorm(nC, mean = muC, sigma = SigmaC)
  
  # Estimate of variance of difference between means
  varhatvec = denomT(XE, XC, rho = rho)
  
  XEbar = colMeans(XE)
  XCbar = colMeans(XC)
  statvec = (XEbar - XCbar)/sqrt(varhatvec)
  
  SigmaEC = SigmaE/nE + SigmaC/nC
  
  # Degrees of freedom for t-distribution
  df = (rtr(SigmaEC %*% SigmaEC) +
           (rtr(SigmaEC))^2)/((rtr(SigmaE %*% SigmaE/(nE^2)) +
                                 (rtr(SigmaE/nE))^2)/(nE - 1) + 
                                (rtr(SigmaC %*% SigmaC/(nC^2)) +
                                   (rtr(SigmaC/nC))^2)/(nC - 1))
  
  # Simulated p-values
  pvals[1:9] = 1 - pt(statvec, df = df)
  pvals[10] = 1 - pnorm(rnorm(1, mean=mu.secondary))
  
  ### FWER
  
  # k = 1
  out.kgMCP1[i,] = gate_kHolm(k=1, pvalues = pvals, alpha=alpha)$rejected
  out.kgMCP.aug1[i,] = gate_kHolmAug(k=1, pvalues = pvals, alpha=alpha)$rejected
  
  # k = 2
  out.kgMCP2[i,] = gate_kHolm(k=2, pvalues = pvals, alpha=alpha)$rejected
  out.kgMCP.aug2[i,] = gate_kHolmAug(k=2, pvalues = pvals, alpha=alpha)$rejected
  
  # k = 3
  out.kgMCP3[i,] = gate_kHolm(k=3, pvalues = pvals, alpha=alpha)$rejected
  out.kgMCP.aug3[i,] = gate_kHolmAug(k=3, pvalues = pvals, alpha=alpha)$rejected
  
  
  ### FDP
  
  # gamma = 0.1
  out.FDP.gMCP1[i,] = gate_FDP_Holm(pvalues = pvals, gamma.fdp=0.1,
                                   alpha=alpha)$rejected
  
  out.FDP.gMCP.aug1[i,] = gate_FDP_HolmAug(pvalues = pvals, gamma.fdp=0.1,
                                   alpha=alpha)$rejected
  
  # gamma = 0.2
  out.FDP.gMCP2[i,] = gate_FDP_Holm(pvalues = pvals, gamma.fdp=0.2,
                                    alpha=alpha)$rejected
  
  out.FDP.gMCP.aug2[i,] = gate_FDP_HolmAug(pvalues = pvals, gamma.fdp=0.2,
                                           alpha=alpha)$rejected
  
  # gamma = 0.3
  out.FDP.gMCP3[i,] = gate_FDP_Holm(pvalues = pvals, gamma.fdp=0.3,
                                    alpha=alpha)$rejected
  
  out.FDP.gMCP.aug3[i,] = gate_FDP_HolmAug(pvalues = pvals, gamma.fdp=0.3,
                                           alpha=alpha)$rejected
}

### Calculate average power

# FWER
power.kgMCP1 = colSums(out.kgMCP1)/nsim
power.kgMCP.aug1 = colSums(out.kgMCP.aug1)/nsim
power.kgMCP2 = colSums(out.kgMCP2)/nsim
power.kgMCP.aug2 = colSums(out.kgMCP.aug2)/nsim
power.kgMCP3 = colSums(out.kgMCP3)/nsim
power.kgMCP.aug3 = colSums(out.kgMCP.aug3)/nsim

# FDP
power.FDP.gMCP1 = colSums(out.FDP.gMCP1)/nsim
power.FDP.gMCP.aug1 = colSums(out.FDP.gMCP.aug1)/nsim
power.FDP.gMCP2 = colSums(out.FDP.gMCP2)/nsim
power.FDP.gMCP.aug2 = colSums(out.FDP.gMCP.aug2)/nsim
power.FDP.gMCP3 = colSums(out.FDP.gMCP3)/nsim
power.FDP.gMCP.aug3 = colSums(out.FDP.gMCP.aug3)/nsim


### Display Table 5

# FWER
round(100*power.kgMCP1)
round(100*power.kgMCP.aug1)
round(100*power.kgMCP2)
round(100*power.kgMCP.aug2)
round(100*power.kgMCP3)
round(100*power.kgMCP.aug3)

# FDP
round(100*power.FDP.gMCP1) 
round(100*power.FDP.gMCP.aug1)
round(100*power.FDP.gMCP2)
round(100*power.FDP.gMCP.aug2)
round(100*power.FDP.gMCP3)
round(100*power.FDP.gMCP.aug3) 
