library(targets)
tar_option_set(packages = c("stdReg", "flexsurv", "survival"))

source("R/data-generation.R")

settings <- function() {
  
  data.frame(
    generation = c("linear", "linear", "log_poisson", "linear", "logit_binomial", 
                   "survival", "survival", "survival", "survival", "survival", "survival", "survival", "survival"), 
    analysis = c("ols", "ols_weighted", "glm_weighted", "ols_weighted_standardized", "glm_weighted_standardized",
                 "survspline",  "weibullPH", "exponential",  "survspline",  "weibullPH", "exponential", 
                 "coxph", "coxph"), 
    coefZ = c(2, 2, .6, 2, .4, 0, 0, 0, .6, .6, .6, 0, .6)
  )
  
}

list(
  tar_target(settings, settings()), 
  tar_target(numeric_results, 
             run_simulation(settings$generation, settings$analyis, settings$coefZ), 
             pattern = map(settings$generation, settings$analysis, settings$coefZ), 
             iteration = "vector"
             )
)

