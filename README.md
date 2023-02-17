# Propensity score weighting plus regression adjustment does not equal doubly robust
## Reproducibility materials

### What is this?

This contains a series of R scripts to reproduce the simulation results in the above-titled paper. 

### How do I run this code? 

#### Prerequisites
You will need a recent version of R, along with the following packages `c("targets", "tarchetypes", "stdReg", "flexsurv", "survival")`. 

#### Steps

1. Download the repository
2. Source the `_targets.R` script
3. Run `tar_make()` or `tar_make_clustermq()` to run the simulation
4. Inspect the results

