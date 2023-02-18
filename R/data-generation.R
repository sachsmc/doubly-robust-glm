#' Generate data under a linear regression model
#' @param n Sample size
#' @param coefZ Numeric value for the Z coefficient


generate_linear <- function(n = 500, coefZ = 2) {
  
  
  C <- rnorm(n)
  D <- rnorm(n)
  Z <- rbinom(n, 1, prob = plogis(-1 + 1 * C + .4 * C^2 - 1 * D))
  
  linpred <- -2 + coefZ * Z + 1 * C + 1 * D
  
  Y <- rnorm(n,  linpred, sd = 1)
  
  data.frame(Y, Z, C, D)
  
}

truev_linear <- function(coefZ = 2) {
  
  coefZ
  
}


generate_linear_interaction <- function(n = 500, coefZ = 2) {
  
  
  C <- rnorm(n)
  D <- rnorm(n)
  Z <- rbinom(n, 1, prob = plogis(-1 + 1 * C + .4 * C^2 - 1 * D))
  
  linpred <- -2 + coefZ * Z + .5 * Z * C + 1 * C + 1 * D
  
  Y <- rnorm(n,  linpred, sd = 1)
  
  data.frame(Y, Z, C, D)
  
}

truev_linear_interaction <- function(coefZ = 2) {
  
  n <- 1e6
  C <- rnorm(n)
  D <- rnorm(n)
  
  linpred <-  -2 + coefZ * 1 + .5 * 1 * C + 1 * C + 1 * D - 
    (-2  + 1 * C + 1 * D)
  
  mean(linpred)
  
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
  
  data.frame(Y, Z, C, D)
  
}

truev_log_binomial <- function(coefZ = .4) {
  
  n <- 1e6
  C <- rnorm(n)
  D <- rnorm(n)
  
  linpred <- exp(-2 + coefZ + .1 * C + .2 * D) - exp(-2 + .1 * C + .2 * D)
  
  mean(linpred)
  
  
}

#' Generate data for a binomial glm with identity link
#' 
#' @param n Sample size
#' @param coefZ Numeric value for the Z coefficient

generate_identity_binomial <- function(n = 500, coefZ = .2) {
  
  
  C <- rnorm(n)
  D <- rnorm(n)
  Z <- rbinom(n, 1, prob = plogis(-1 + 1 * C + .4 * C^2 - 1 * D))
  
  linpred <- .3 + coefZ * Z + .05 * C + .02 * D
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

generate_logit_binomial <- function(n = 500, coefZ = 2) {
  
  
  C <- rnorm(n)
  D <- rnorm(n)
  Z <- rbinom(n, 1, prob = plogis(-1 + 1 * C + .4 * C^2 - 1 * D))
  
  linpred <- -2 + coefZ * Z + 1 * C + 1 * D
  meanparm <- plogis(linpred)
  
  Y <- rbinom(n, 1, meanparm)
  
  data.frame(Y, Z, C, D)
  
}


truev_logit_binomial <- function(coefZ = .4) {
  
  n <- 1e6
  C <- rnorm(n)
  D <- rnorm(n)
  
  linpred <- plogis(-2 + coefZ + .1 * C + .2 * D) - plogis(-2 + .1 * C + .2 * D)
  
  mean(linpred)
  
  
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
  
  data.frame(Y, Z, C, D)
  
}

truev_log_poisson <- function(coefZ = .4) {
  
  n <- 1e6
  C <- rnorm(n)
  D <- rnorm(n)
  
  mean(exp(-2 + coefZ * 1 + 1 * C + 2 * D) - exp(-2  + 1 * C + 2 * D))
  
}

#' Generate data for a survival outcomes
#' 
#' @param n Sample size
#' @param coefZ Numeric value for the Z coefficient
generate_survival <- function(n = 500, coefZ = log(2)) {
  
  C <- rnorm(n)
  D <- rnorm(n)
  Z <- rbinom(n, 1, prob = plogis(-1 + 1 * C + .4 * C^2 - 1 * D))
  scale <- exp(coefZ * Z + C + D)
  
  Y.true <- rweibullPH(n, shape = 2, scale = scale)
  Cens <- rweibull(n, 1.4 * (1 / 1.7))
  
  
  E<-rexp(n, exp(1+coefZ * Z + C + D))
  
  Y.obs <- pmin(Y.true, Cens)
  D.ind <- 1.0 * (Y.true <= Cens)
  Ye.obs<-pmin(E, Cens)
  E.ind<-1.0 * (E <= Cens)
  data.frame(Y.obs, D.ind, Z, C, D,Ye.obs,E.ind)
  
}


truev_survival <- function(coefZ = log(2)) {

  n <- 1e6
  C <- rnorm(n)
  D <- rnorm(n)
  Z <- rbinom(n, 1, prob = plogis(-1 + 1 * C + .4 * C^2 - 1 * D))
  scale1 <- exp(coefZ * 1 + C + D)
  scale0 <- exp(C + D)
  
  c(weib = mean(pweibullPH(1, shape = 2, scale = scale1, lower.tail = FALSE) - pweibullPH(1, shape = 2, scale = scale0, lower.tail = FALSE)),
    exp = mean(pexp(1, exp(1 + coefZ + C + D), lower.tail = FALSE) - pexp(1, exp(1 + C + D), lower.tail = FALSE)))
  
}



