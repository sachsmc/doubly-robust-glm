library(targets)
library(tarchetypes)
tar_option_set(packages = c("stdReg", "flexsurv", "survival"))

source("R/data-generation.R")
source("R/analysis.R")

## add log binomial, linear with binary outcome, poisson with binary outcome

settings <-  data.frame(
    generation = c("linear", "linear", "log_poisson", "linear", 
                   "log_binomial", "log_binomial", "identity_binomial", "logit_binomial", 
                   "survival", "survival", "survival", "survival", "survival", "survival", "survival", "survival"), 
    analysis = c("ols", "ols_weighted", "poisson_weighted", "ols_weighted_standardized", 
                "log_binomial_weighted_standardized", "poisson_weighted_standardized",
                "ols_weighted_standardized",
                "logit_binomial_weighted_standardized", 
                 "survspline",  "weibullPH", "exponential",  "survspline",  "weibullPH", "exponential", 
                 "coxph", "coxph"), 
    coefZ = c(2, 2, .6, 2, .2, .2, .2, .4, 0, 0, 0, .6, .6, .6, 0, .6)
  )
  

target_runs <- tar_map(settings, 
        tar_target(simulate, run_simulation(generation, analysis, coefZ))
        )

combined_runs <- tar_combine(simulation_results, target_runs[["simulate"]])

c(target_runs, combined_runs)

#tar_make()
