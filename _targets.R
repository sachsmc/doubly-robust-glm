library(targets)
library(tarchetypes)
tar_option_set(packages = c("stdReg", "ggplot2", "data.table", "parallel", "sandwich"))

source("R/data-generation.R")
source("R/analysis.R")
source("R/results-summary.R")
source("R/variance-estimation.R")


B <- 20  ## number of simulation replicates
cores <- 10 ## number of cpus to use
boots <- 20 ## number of bootstrap replicates

## add log binomial, linear with binary outcome, poisson with binary outcome

settings <-  data.frame(
    generation = c("linear", ## ols settings
                   "linear",
                   "log_poisson", "logit_binomial", "inverse_gaussian", ## glm canon settings
                   "logit_binomial",  ## non canon setting
                   "linear", "logit_binomial", ## funk comparison
                   "linear", "logit_binomial", ## funk comparison
                   "linear", "logit_binomial", ## funk comparison
                   "linear", "logit_binomial", ## funk comparison
                   "linearodd"
                   ), 
    analysis = c("ols_weighted", "ols_weighted_standardized",  # ols settings
                "poisson_weighted_standardized", "logit_binomial_weighted_standardized", "inverse_gaussian_weighted_standardized",
                "log_binomial_weighted_standardized",
                "linear_compare", "logit_compare", 
                "linear_compare", "logit_compare", 
                "linear_compare", "logit_compare", 
                "linear_compare", "logit_compare", 
                "ols_weighted_standardized_odd"),
    n = c(rep(2000,6),2000, 2000, 1000, 1000, 500,500, 100,100, 2000),
    coefZ = c(2, 2,2,2, 200, 20, 2,5,2, 5,2,5,2,5, 2)
  )
  

target_runs <- tar_map(settings, 
        tar_target(simulate, run_simulation(generation, analysis, coefZ, n))
        )

combined_runs <- tar_combine(simulation_results, target_runs[["simulate"]])

tables <- list(
  tar_target(sim_summary, results_table(subset(simulation_results, analysis %in% 
                                                    c("ols_weighted", "ols_weighted_standardized", "poisson_weighted_standardized", 
                                                      "logit_binomial_weighted_standardized", 
                                                      "inverse_gaussian_weighted_standardized", "log_binomial_weighted_standardized")))), 
  tar_target(funk_compare, stddev_table(subset(simulation_results, grepl("compare", analysis)))), 
  tar_target(odd_text, results_table(subset(simulation_results, analysis == "ols_weighted_standardized_odd"))),
  tar_target(sdev_table, eif_table(subset(simulation_results, analysis %in% 
                                           c("ols_weighted_standardized", "poisson_weighted_standardized", 
                                             "logit_binomial_weighted_standardized", 
                                             "inverse_gaussian_weighted_standardized", "log_binomial_weighted_standardized")))),
  tar_quarto(example, "example-analysis.qmd")
  
)

c(target_runs, combined_runs, tables)

