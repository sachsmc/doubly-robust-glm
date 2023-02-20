library(targets)
library(tarchetypes)
tar_option_set(packages = c("stdReg", "ggplot2", "data.table"), 
               cue = tar_cue("always"))

source("R/data-generation.R")
source("R/analysis.R")
source("R/results-summary.R")

## add log binomial, linear with binary outcome, poisson with binary outcome

settings <-  data.frame(
    generation = c("linear", "linear", "log_poisson", "logit_binomial", ## ols settings
                   "log_poisson", "logit_binomial", "inverse_gaussian", ## glm canon settings
                   "logit_binomial" ## non canon setting
                   ), 
    analysis = c("ols", "ols_weighted", "ols_weighted_standardized", "ols_weighted_standardized", # ols settings
                "poisson_weighted_standardized", "logit_binomial_weighted_standardized", "inverse_gaussian_weighted_standardized",
                "log_binomial_weighted_standardized"), 
    coefZ = c(2, 2, 2, 5, 2,2, 200, 20)
  )
  
target_runs <- tar_map(settings, 
        tar_target(simulate, run_simulation(generation, analysis, coefZ))
        )

combined_runs <- tar_combine(simulation_results, target_runs[["simulate"]])

tables <- list(
  tar_target(linear_summary, results_table(subset(simulation_results, grepl("ols", analysis)))), 
  tar_target(glm_summary, results_table(subset(simulation_results, grepl("binomial|poisson|inverse\\.gaussian", analysis))))
)

c(target_runs, combined_runs, tables)

## linear regression with no interaction is dr, under any data generating mech including continuous, binary, and poisson, maybe survival with pseudo observations 
## canonical link glm is dr, poisson, binomial, and inverse gaussian (add)
## noncanonical links do not work: log binomial with logit generation

## survival stuff, estimand should be risk difference or rmean difference
## add generation for non PH, add tsiatis and chen method for rmean
## cox and parametric are not dr except under the null,
## pseudo obs or direct binomial works under independent censoring
