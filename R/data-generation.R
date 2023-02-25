#' Generate data under a linear regression model
#' @param n Sample size
#' @param coefZ Numeric value for the Z coefficient


generate_linear <- function(n = 2000,coefZ = 2) {
  
  
  C <- rnorm(n)
  D <- rnorm(n, mean = 1)
  Z <- rbinom(n, 1, prob = plogis(.2*(-2 + 2 * C + 1.4 * C^2 + 2 * D)))
  
  linpred <- -2 + coefZ * Z + 1 * C + .4 * C^2 + 1.5 * D
  
  Y <- rnorm(n,  linpred, sd = 1)
  
  data.frame(Y, Z, C, D)
  
}




generate_linearodd <- function(n = 2000,coefZ = 2) {
  
  
  C <- rnorm(n)
  D <- rbinom(n, 1, .35)
  Z <- rbinom(n, 1, prob = plogis(.2*(-2 + 20 * C + 1.4 * C^2 + 2 * D)))
  
  linpred <- -2 + coefZ * Z + 1 * C + 4.5 * D
  
  Y <- rnorm(n,  linpred, sd = 1)
  
  data.frame(Y, Z, C, D)
  
}

truev_linear <- function(coefZ = 2) {
  
  coefZ
  
}


truev_linearodd <- function(coefZ = 2) {
  
  coefZ
  
}

generate_inverse_gaussian <- function(n = 2000,coefZ = 200) {
  
  C <- rnorm(n)
  D <- rnorm(n, mean = 1)
  Z <- rbinom(n, 1, prob = plogis(.2*(-2 + 2 * C + 1.4 * C^2 + 2 * D)))
  
  linpred <- 50 + coefZ * Z + 4 * C + 10 * C^2 + 5 * D
  meanparm <- 1/sqrt(linpred)
  
  shape <- 2
  nu2 <- rnorm(n)^2
  xx <- meanparm + meanparm ^ 2 * nu2 / (2 * shape) - 
    meanparm / (2 * shape) * sqrt(4 * meanparm * shape * nu2 + meanparm^2 * nu2^2)
  
  zz <- runif(n)
  Y <- ifelse(zz <= meanparm / (meanparm + xx), xx, meanparm^2 / xx)
  
  data.frame(Y, Z, C, D)
  
}

truev_inverse_gaussian <- function(coefZ = 200) {
  
  C <- rnorm(1e6)
  D <- rnorm(1e6, mean = 1)
  Z <- rbinom(1e6, 1, prob = plogis(.2*(-2 + 2 * C + 1.4 * C^2 + 2 * D)))
  
  linpred1 <- 50 + coefZ * 1 + 4 * C + 10 * C^2 + 5 * D
  linpred0 <- 50 + coefZ * 0 + 4 * C + 10 * C^2 + 5 * D
  
  mean(1/sqrt(linpred1) - 1/sqrt(linpred0))
  
  
}


#' Generate data for a binomial glm with log link
#' 
#' @param n Sample size
#' @param coefZ Numeric value for the Z coefficient

generate_log_binomial <- function(n = 2000,coefZ = .4) {
  
  
  C <- rnorm(n)
  D <- rnorm(n, mean = 1)
  Z <- rbinom(n, 1, prob = plogis(.2*(-2 + 2 * C + 1.4 * C^2 + 2 * D)))
  
  linpred <- -4 + coefZ * Z + .1 * C + .04 * C^2 + .8 * D
  meanparm <- pmin(1, exp(linpred)) # just in case
  
  Y <- rbinom(n, 1, meanparm)
  
  data.frame(Y, Z, C, D)
  
}

truev_log_binomial <- function(coefZ = .4) {
  
  n <- 1e6
  C <- rnorm(n)
  D <- rnorm(n, mean = 1)
  
  linpred <- exp(-4 + coefZ  + .1 * C + .04 * C^2 + .8 * D) - 
    exp(-4 + .1 * C + .04 * C^2 + .8 * D)
  
  mean(linpred)
  
  
}

#' Generate data for a binomial glm with identity link
#' 
#' @param n Sample size
#' @param coefZ Numeric value for the Z coefficient

generate_identity_binomial <- function(n = 2000,coefZ = .2) {
  
  
  C <- rnorm(n)
  D <- rnorm(n, mean = 1)
  Z <- rbinom(n, 1, prob = plogis(.2*(-2 + 2 * C + 1.4 * C^2 + 2 * D)))
  
  linpred <- .3 + coefZ * Z + .02 * C + .005 * C^2 + .08 * D
  meanparm <- pmax(0, pmin(1, linpred)) # just in case
  
  Y <- rbinom(n, 1, meanparm)
  
  data.frame(Y, Z, C, D)
  
}

truev_identity_binomial <- function(coefZ = .2) {
  
 coefZ
  
}


#' Generate data for a binomial glm with logit link
#' 
#' @param n Sample size
#' @param coefZ Numeric value for the Z coefficient

generate_logit_binomial <- function(n = 2000,coefZ = 2) {
  
  
  C <- rnorm(n)
  D <- rnorm(n, mean = 1)
  Z <- rbinom(n, 1, prob = plogis(.2*(-2 + 2 * C + 1.4 * C^2 + 2 * D)))
  
  linpred <- -2 + coefZ * Z + 1 * C + 4 * D + 1 * D * Z
  meanparm <- plogis(linpred)
  
  Y <- rbinom(n, 1, meanparm)
  
  data.frame(Y, Z, C, D)
  
}


truev_logit_binomial <- function(coefZ = 2) {
  
  n <- 1e6
  C <- rnorm(n)
  D <- rnorm(n, mean = 1)
  
  linpred <- plogis(-2 + coefZ + 1 * C + 4 * D + D) - plogis(-2 + 1 * C + 4 * D)
  
  mean(linpred)
  
  
}


#' Generate data for a poisson glm with log link
#' 
#' @param n Sample size
#' @param coefZ Numeric value for the Z coefficient

generate_log_poisson <- function(n = 2000,coefZ = 1) {
  
  
  C <- rnorm(n)
  D <- rnorm(n, mean = 1)
  Z <- rbinom(n, 1, prob = plogis(.2*(-2 + 2 * C + 1.4 * C^2 + 2 * D)))
  
  linpred <- 0 + coefZ * Z + .1 * C + .05 * C^2 + .4 * D
  meanparm <- exp(linpred)
  
  Y <- rpois(n, meanparm)
  
  data.frame(Y, Z, C, D)
  
}

truev_log_poisson <- function(coefZ = 1) {
  
  n <- 1e6
  C <- rnorm(n)
  D <- rnorm(n, mean = 1)
  
  mean(exp(0 + coefZ * 1 + .1 * C + .05 * C^2 + .4 * D) - exp(0 + coefZ * 0 + .1 * C + .05 * C^2 + .4 * D))
  
}
