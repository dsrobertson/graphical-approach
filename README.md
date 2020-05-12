# graphical-approach

Code for the paper **Graphical approaches for the control of generalised error rates** (2020),
available at https://arxiv.org/abs/2004.01759


## Description of R files
- **PD_study.R** - code for Section 6.1: an analysis of the Pharmacodynamic study
reported by Ferber *et al.* (2011)

- **Pre_RELAX_study.R** - code for Section 6.2: simulation study of the Pre-RELAX-AHF
trial reported by Teerlink *et al.* (2009)

- **kFWER.R** - implements the augmented and generalised graphical approaches for
k-FWER control

- **FDP.R** - implements the augmented and generalised graphical approaches for FDP
control

- **weight_fn.R** - implements the graphical weighting strategy of Bretz *et al.* (2011)

- **gatekeeping.R** - implements the gatekeeping strategy requiring rejection of at
least 6 out of 9 hypotheses, with a single secondary hypothesis
