library(targets)
library(tarchetypes)
tar_option_set(packages = c("stdReg", "ggplot2", "data.table", "parallel", "sandwich"))

source("R/data-generation.R")
source("R/analysis.R")
source("R/results-summary.R")
source("R/variance-estimation.R")

## add log binomial, linear with binary outcome, poisson with binary outcome

settings <-  data.frame(
    generation = c("linear", ## ols settings
                   "log_poisson", "logit_binomial", "inverse_gaussian", ## glm canon settings
                   "logit_binomial",  ## non canon setting
                   "linear", "logit_binomial", ## funk comparison
                   "linear", "logit_binomial", ## funk comparison
                   "linear", "logit_binomial", ## funk comparison
                   "linearodd"
                   ), 
    analysis = c("ols_weighted",  # ols settings
                "poisson_weighted_standardized", "logit_binomial_weighted_standardized", "inverse_gaussian_weighted_standardized",
                "log_binomial_weighted_standardized", 
                "linear_compare", "logit_compare", 
                "linear_compare", "logit_compare", 
                "linear_compare", "logit_compare", 
                "ols_weighted_standardized_odd"),
    n = c(rep(2000,5), 1000, 1000, 500,500, 100,100, 2000),
    coefZ = c(2, 2,2, 200, 20, 2, 5,2,5,2,5, 2)
  )
  
#settings <- subset(settings, grepl("compare", analysis))
target_runs <- tar_map(settings, 
        tar_target(simulate, run_simulation(generation, analysis, coefZ, n))
        )

combined_runs <- tar_combine(simulation_results, target_runs[["simulate"]])

tables <- list(
  tar_target(linear_summary, results_table(subset(simulation_results, analysis %in% c("ols", "ols_weighted", "ols_weighted_standardized", "ols_weighted_standardized")))), 
  tar_target(glm_summary, results_table(subset(simulation_results, grepl("binomial|poisson|inverse", analysis))))
)

c(target_runs, combined_runs, tables)
#target_runs


## linear regression with no interaction is dr, under any data generating mech including continuous, binary, and poisson, maybe survival with pseudo observations 
## canonical link glm is dr, poisson, binomial, and inverse gaussian (add)
## noncanonical links do not work: log binomial with logit generation

## survival stuff, estimand should be risk difference or rmean difference
## add generation for non PH, add tsiatis and chen method for rmean
## cox and parametric are not dr except under the null,
## pseudo obs or direct binomial works under independent censoring
