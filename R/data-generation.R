#' Generate data under a linear regression model
#' @param n Sample size
#' @param coefZ Numeric value for the Z coefficient

generate_linear <- function(n = 500, coefZ = 2) {
  
  
  C <- rnorm(n)
  D <- rnorm(n)
  Z <- rbinom(n, 1, prob = plogis(-1 + 1 * C + .4 * C^2 - 1 * D))
  
  linpred <- -2 + coefZ * Z + 1 * C + 1 * D
  
  Y <- rnorm(n,  linpred, sd = 1)
  
  data <- data.frame(Y, Z, C, D)
  
}

#' Generate data for a binomial glm with log link
#' 
#' @param n Sample size
#' @param coefZ Numeric value for the Z coefficient

generate_log_binomial <- function(n = 500, coefZ = .4) {
  
  
  C <- rnorm(n)
  D <- rnorm(n)
  Z <- rbinom(n, 1, prob = plogis(-1 + 1 * C + .4 * C^2 - 1 * D))
  
  linpred <- -2 + coefZ * Z + .1 * C + .2 * D
  meanparm <- exp(linpred)
  
  Y <- rbinom(n, 1, meanparm)
  
  data <- data.frame(Y, Z, C, D)
  
}

#' Generate data for a binomial glm with logit link
#' 
#' @param n Sample size
#' @param coefZ Numeric value for the Z coefficient

generate_logit_binomial <- function(n = 500, coefZ = 2) {
  
  
  C <- rnorm(n)
  D <- rnorm(n)
  Z <- rbinom(n, 1, prob = plogis(-1 + 1 * C + .4 * C^2 - 1 * D))
  
  linpred <- -2 + coefZ * Z + 1 * C + 1 * D
  meanparm <- plogis(linpred)
  
  Y <- rbinom(n, 1, meanparm)
  
  data <- data.frame(Y, Z, C, D)
  
}


#' Generate data for a poisson glm with log link
#' 
#' @param n Sample size
#' @param coefZ Numeric value for the Z coefficient

generate_log_poisson <- function(n = 500, coefZ = .4) {
  
  
  C <- rnorm(n)
  D <- rnorm(n)
  Z <- rbinom(n, 1, prob = plogis(-1 + 1 * C + .4 * C^2 - 1 * D))
  
  linpred <- -2 + coefZ * Z + 1 * C + 2 * D
  meanparm <- exp(linpred)
  
  Y <- rpois(n, meanparm)
  
  data <- data.frame(Y, Z, C, D)
  
}



#' Generate data for a survival outcomes
#' 
#' @param n Sample size
#' @param coefZ Numeric value for the Z coefficient
generate_data <- function(n = 500, coefZ = log(2)) {
  
  C <- rnorm(n)
  D <- rnorm(n)
  Z <- rbinom(n, 1, prob = plogis(-1 + 1 * C + .4 * C^2 - 1 * D))
  scale <- exp(loghr.Z * Z + C + D)
  
  Y.true <- rweibullPH(n, shape = 2, scale = scale)
  Cens <- rweibull(n, 1.4 * (1 / 1.7))
  
  
  E<-rexp(n, exp(1+coefZ * Z + C + D))
  
  Y.obs <- pmin(Y.true, Cens)
  D.ind <- 1.0 * (Y.true <= Cens)
  Ye.obs<-pmin(E, Cens)
  E.ind<-1.0 * (E <= Cens)
  data.frame(Y.obs, D.ind, Z, C, D,Ye.obs,E.ind)
  
}

