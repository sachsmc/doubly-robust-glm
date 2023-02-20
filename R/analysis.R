#' Run the simulation according to the settings
#' 
#' @param generation linear, log_poisson, log_binomial, logit_binomial, or survival
#' @param analysis ols, ols_weighted, ols_weighted_standardized, 
#'                  glm_weighted, glm_weighted_standardized, survspline, weibullPH, exponential, or coxph
#' @param coefZ numeric value for the regression coefficient
#' @param B number of simulation replicates
#' @param n sample size in each generation
#'                  

B <- 200
n <- 2000

run_simulation <- function(generation, analysis, coefZ) {
  
  ## number of simulation replicates
  
  generate_data <- get(paste0("generate_", generation))
  analyze_data <- get(paste0("analyze_", analysis))
  true_values <- get(paste0("truev_", generation))
  
  truev <- true_values(coefZ)
  
  if(generation == "survival") {
    truev <- if(analysis == "exponential") truev[2] else truev[1]
  }
  
  results <- NULL
  for(i in 1:B) {
    
    datin <- generate_data(n = n, coefZ = coefZ) 
    res <- tryCatch(analyze_data(datin), error = function(e) data.frame(est = NA, type = "failed"))
    results <- rbind(results, res)
    
  }
  
  results$setting <- generation
  results$analysis <- analysis
  results$true_value <- truev
  results$coefZ <- coefZ
  results
  
}


analyze_ols <- function(data) {
  
  fit1 <- lm(Y ~ Z + C + I(C^2) + D, data = data)
  fit2 <- lm(Y ~ Z + C, data = data)
  data.frame(est = c(coefficients(fit1)[2], coefficients(fit2)[2]), 
             type = c("right outcome no weights", "wrong outcome no weights"))
  
}

analyze_ols_weighted <- function(data) {
  
  
  zmod1 <- glm(Z ~ C + I(C^2) + D, data = data, family = binomial)
  phat <- predict(zmod1, type = "response")
  #phat <- plogis(-1 + 1 * data$C + .4 * data$C^2 - 1 * data$D)
  data$W1 <- data$Z / phat + (1 - data$Z) / (1 - phat)
  
  phatWR <- predict(glm(Z ~ C, data= data, family = binomial), type = "response")
  #phatWR <- plogis(-1 + 0.5 * data$C)
  data$W2 <- data$Z / phatWR + (1 - data$Z) / (1 - phatWR)
  
  fit1 <- lm(Y ~ Z + C, data = data, weights = W1)
  fit2 <- lm(Y ~ Z + C + I(C^2) + D, data = data, weights = W2)
  
  data.frame(est = c(coefficients(fit1)[2], coefficients(fit2)[2]), 
             type = c("wrong outcome right weights", "right outcome wrong weights"))
  
}


analyze_ols_weighted_standardized <- function(data) {
  
  
  zmod1 <- glm(Z ~ C + I(C^2) + D, data = data, family = binomial)
  phat <- predict(zmod1, type = "response")
  #phat <- plogis(-1 + 1 * data$C + .4 * data$C^2 - 1 * data$D)
  data$W1 <- data$Z / phat + (1 - data$Z) / (1 - phat)
  
  phatWR <- predict(glm(Z ~ C, data= data, family = binomial), type = "response")
  #phatWR <- plogis(-1 + 0.5 * data$C)
  data$W2 <- data$Z / phatWR + (1 - data$Z) / (1 - phatWR)
  
  fit1 <- glm(Y ~ Z + C, data = data, weights = W1, family = "gaussian")
  fit2 <- glm(Y ~ Z + C + I(C^2) + D, data = data, weights = W2, family = "gaussian")
  fit3 <- glm(Y ~ Z + C, data = data, weights = W2, family = "gaussian")
  
  fit1$prior.weights <- fit2$prior.weights <- rep(1, nrow(data))
  
  ests <- sapply(list(fit1, fit2, fit3), 
                 \(fit) {
                   summary(stdGlm(fit, data, "Z"), contrast = "difference", reference = 0)$est.table[2, 1]
                 })
  
  
  
  data.frame(est = ests, 
             type = c("wrong outcome right weights", "right outcome wrong weights", "wrong both"))
  
}


analyze_poisson_weighted <- function(data) {
  
  zmod1 <- glm(Z ~ C + I(C^2) + D, data = data, family = binomial)
  phat <- predict(zmod1, type = "response")
  #phat <- plogis(-1 + 1 * data$C + .4 * data$C^2 - 1 * data$D)
  data$W1 <- data$Z / phat + (1 - data$Z) / (1 - phat)
  
  phatWR <- predict(glm(Z ~ C, data= data, family = binomial), type = "response")
  #phatWR <- plogis(-1 + 0.5 * data$C)
  data$W2 <- data$Z / phatWR + (1 - data$Z) / (1 - phatWR)
  
  fit1 <- glm(Y ~ Z + C, data = data, weights = W1, family = "poisson")
  fit2 <- glm(Y ~ Z + C + I(C^2) + D, data = data, weights = W2, family = "poisson")
  
  data.frame(est = sapply(list(fit1, fit2), \(x) coefficients(x)[2]), 
             type = c("wrong outcome right weights", "right outcome wrong weights"))
  
}


analyze_poisson_weighted_standardized <- function(data) {
  
  zmod1 <- glm(Z ~ C + I(C^2) + D, data = data, family = binomial)
  phat <- predict(zmod1, type = "response")
  #phat <- plogis(-1 + 1 * data$C + .4 * data$C^2 - 1 * data$D)
  data$W1 <- data$Z / phat + (1 - data$Z) / (1 - phat)
  
  phatWR <- predict(glm(Z ~ C, data= data, family = binomial), type = "response")
  #phatWR <- plogis(-1 + 0.5 * data$C)
  data$W2 <- data$Z / phatWR + (1 - data$Z) / (1 - phatWR)
  
  fit1 <- glm(Y ~ Z + C, data = data, weights = W1, family = "poisson")
  fit2 <- glm(Y ~ Z + C + I(C^2) + D, data = data, weights = W2, family = "poisson")
  fit3 <- glm(Y ~ Z + C, data = data, weights = W2, family = "poisson")
  
  
  fit1$prior.weights <- fit2$prior.weights <- rep(1, nrow(data))
  
  ests <- sapply(list(fit1, fit2, fit3), 
                 \(fit) {
                   summary(stdGlm(fit, data, "Z"), contrast = "difference", reference = 0)$est.table[2, 1]
                 })
  
  data.frame(est = ests, 
             type = c("wrong outcome right weights", "right outcome wrong weights", "wrong both"))
  
}


analyze_inverse_gaussian_weighted_standardized <- function(data) {
  
  zmod1 <- glm(Z ~ C + I(C^2) + D, data = data, family = binomial)
  phat <- predict(zmod1, type = "response")
  #phat <- plogis(-1 + 1 * data$C + .4 * data$C^2 - 1 * data$D)
  data$W1 <- data$Z / phat + (1 - data$Z) / (1 - phat)
  
  phatWR <- predict(glm(Z ~ C, data= data, family = binomial), type = "response")
  #phatWR <- plogis(-1 + 0.5 * data$C)
  data$W2 <- data$Z / phatWR + (1 - data$Z) / (1 - phatWR)
  
  fit1 <- glm(Y ~ Z + C, data = data, weights = W1, family = "inverse.gaussian")
  fit2 <- glm(Y ~ Z + C + I(C^2) + D, data = data, weights = W2, family = "inverse.gaussian")
  fit3 <- glm(Y ~ Z + C, data = data, weights = W2, family = "inverse.gaussian")
  
  
  fit1$prior.weights <- fit2$prior.weights <- rep(1, nrow(data))
  
  ests <- sapply(list(fit1, fit2, fit3), 
                 \(fit) {
                   summary(stdGlm(fit, data, "Z"), contrast = "difference", reference = 0)$est.table[2, 1]
                 })
  
  data.frame(est = ests, 
             type = c("wrong outcome right weights", "right outcome wrong weights", "wrong both"))
  
}




analyze_log_binomial_weighted_standardized <- function(data) {
  
  #data <- generate_log_binomial()
  zmod1 <- glm(Z ~ C + I(C^2) + D, data = data, family = binomial)
  phat <- predict(zmod1, type = "response")
  #phat <- plogis(-1 + 1 * data$C + .4 * data$C^2 - 1 * data$D)
  data$W1 <- data$Z / phat + (1 - data$Z) / (1 - phat)
  
  phatWR <- predict(glm(Z ~ C, data= data, family = binomial), type = "response")
  #phatWR <- plogis(-1 + 0.5 * data$C)
  data$W2 <- data$Z / phatWR + (1 - data$Z) / (1 - phatWR)
  
  fit1 <- glm(Y ~ Z + C, data = data, weights = W1, family = binomial(link = "log"), 
              start = c(log(mean(data$Y)), 0, 0))
  fit2 <- glm(Y ~ Z + C + I(C^2) + D, data = data, weights = W2, 
              family = binomial(link = "log"), start = c(log(mean(data$Y)), 0, 0, 0, 0))
  
  fit1$prior.weights <- fit2$prior.weights <- rep(1, nrow(data))
  
  ests <- sapply(list(fit1, fit2), 
                 \(fit) {
                   summary(stdGlm(fit, data, "Z"), contrast = "difference", reference = 0)$est.table[2, 1]
                 })
  
  data.frame(est = ests, 
             type = c("wrong outcome right weights", "right outcome wrong weights"))
  
}


analyze_identity_binomial_weighted_standardized <- function(data) {
  
  #data <- generate_identity_binomial(n = 5000)
  zmod1 <- glm(Z ~ C + I(C^2) + D, data = data, family = binomial)
  phat <- predict(zmod1, type = "response")
  #phat <- plogis(-1 + 1 * data$C + .4 * data$C^2 - 1 * data$D)
  data$W1 <- data$Z / phat + (1 - data$Z) / (1 - phat)
  
  phatWR <- predict(glm(Z ~ C, data= data, family = binomial), type = "response")
  #phatWR <- plogis(-1 + 0.5 * data$C)
  data$W2 <- data$Z / phatWR + (1 - data$Z) / (1 - phatWR)
  
  fit1 <- glm(Y ~ Z + C, data = data, weights = W1, family = binomial(link = "identity"), 
              start = c(mean(data$Y), 0, 0))
  fit2 <- glm(Y ~ Z + C + I(C^2) + D, data = data, weights = W2, 
              family = binomial(link = "identity"), start = c(mean(data$Y), 0.1, 0, 0, 0), 
              maxit = 50)
  
  fit1$prior.weights <- fit2$prior.weights <- rep(1, nrow(data))
  
  ests <- sapply(list(fit1, fit2), 
                 \(fit) {
                   summary(stdGlm(fit, data, "Z"), contrast = "difference", reference = 0)$est.table[2, 1]
                 })
  
  data.frame(est = ests, 
             type = c("wrong outcome right weights", "right outcome wrong weights"))
  
}



analyze_logit_binomial_weighted_standardized <- function(data) {
  
  zmod1 <- glm(Z ~ C + I(C^2) + D, data = data, family = binomial)
  phat <- predict(zmod1, type = "response")
  #phat <- plogis(-1 + 1 * data$C + .4 * data$C^2 - 1 * data$D)
  data$W1 <- data$Z / phat + (1 - data$Z) / (1 - phat)
  
  phatWR <- predict(glm(Z ~ C, data= data, family = binomial), type = "response")
  #phatWR <- plogis(-1 + 0.5 * data$C)
  data$W2 <- data$Z / phatWR + (1 - data$Z) / (1 - phatWR)
  
  fit1 <- glm(Y ~ Z + C, data = data, weights = W1, family = binomial())
  fit2 <- glm(Y ~ Z + C + I(C^2) + D, data = data, weights = W2, family = binomial())
  fit3 <- glm(Y ~ Z + C, data = data, weights = W2, family = binomial())
  
  fit1$prior.weights <- fit2$prior.weights <- rep(1, nrow(data))
  
  ests <- sapply(list(fit1, fit2, fit3), 
                 \(fit) {
                   summary(stdGlm(fit, data, "Z"), contrast = "difference", reference = 0)$est.table[2, 1]
                 })
  
  data.frame(est = ests, 
             type = c("wrong outcome right weights", "right outcome wrong weights", "wrong both"))
  
}


analyze_probit_binomial_weighted_standardized <- function(data) {
  
  zmod1 <- glm(Z ~ C + I(C^2) + D, data = data, family = binomial)
  phat <- predict(zmod1, type = "response")
  #phat <- plogis(-1 + 1 * data$C + .4 * data$C^2 - 1 * data$D)
  data$W1 <- data$Z / phat + (1 - data$Z) / (1 - phat)
  
  phatWR <- predict(glm(Z ~ C, data= data, family = binomial), type = "response")
  #phatWR <- plogis(-1 + 0.5 * data$C)
  data$W2 <- data$Z / phatWR + (1 - data$Z) / (1 - phatWR)
  
  fit1 <- glm(Y ~ Z + C, data = data, weights = W1, family = binomial(link = "probit"))
  fit2 <- glm(Y ~ Z + C + I(C^2) + D + D:Z, data = data, weights = W2, family = binomial(link = "probit"))
  fit3 <- glm(Y ~ Z + C, data = data, weights = W2, family = binomial(link = "probit"))
  
  fit1$prior.weights <- fit2$prior.weights <- rep(1, nrow(data))
  
  ests <- sapply(list(fit1, fit2, fit3), 
                 \(fit) {
                   summary(stdGlm(fit, data, "Z"), contrast = "difference", reference = 0)$est.table[2, 1]
                 })
  
  data.frame(est = ests, 
             type = c("wrong outcome right weights", "right outcome wrong weights", "wrong both"))
  
}

