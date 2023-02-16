#' Run the simulation according to the settings
#' 
#' @param generation linear, log_poisson, log_binomial, logit_binomial, or survival
#' @param analysis ols, ols_weighted, ols_weighted_standardized, 
#'                  glm_weighted, glm_weighted_standardized, survspline, weibullPH, exponential, or coxph
#' @param coefZ numeric value for the regression coefficient
#'                  

run_simulation <- function(generation, analyis, coefZ) {
  
  generate_data <- get(paste0("generate_", generation))
  analyze_data <- get(paste0("analyze_", analysis))
  true_values <- get(paste0("truev_", generation))
  
  
  
}


analyze_ols <- function(data) {
  
  fit1 <- lm(Y ~ Z + C + D, data = data)
  fit2 <- lm(Y ~ Z, data = data)
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
  
  fit1 <- lm(Y ~ Z, data = data, weights = W1)
  fit2 <- lm(Y ~ Z + C + D, data = data, weights = W2)
  
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
  
  fit1 <- glm(Y ~ Z * C, data = data, weights = W1, family = "gaussian")
  fit2 <- glm(Y ~ Z + C + D, data = data, weights = W2, family = "gaussian")
  
  fit1$weights <- fit2$weights <- rep(1, nrow(data))
  v
  
  
  data.frame(est = ests, 
             type = c("wrong outcome right weights", "right outcome wrong weights"))
  
}


analyze_glm_weighted <- function(data) {
  
  zmod1 <- glm(Z ~ C + I(C^2) + D, data = data, family = binomial)
  phat <- predict(zmod1, type = "response")
  #phat <- plogis(-1 + 1 * data$C + .4 * data$C^2 - 1 * data$D)
  data$W1 <- data$Z / phat + (1 - data$Z) / (1 - phat)
  
  phatWR <- predict(glm(Z ~ C, data= data, family = binomial), type = "response")
  #phatWR <- plogis(-1 + 0.5 * data$C)
  data$W2 <- data$Z / phatWR + (1 - data$Z) / (1 - phatWR)
  
  fit1 <- glm(Y ~ Z, data = data, weights = W1, family = "poisson")
  fit2 <- glm(Y ~ Z + C + D, data = data, weights = W2, family = "poisson")
  
  data.frame(est = sapply(list(fit1, fit2), \(x) coefficients(x)[2]), 
             type = c("wrong outcome right weights", "right outcome wrong weights"))
  
}

analyze_glm_weighted_standardized <- function(data) {
  
  zmod1 <- glm(Z ~ C + I(C^2) + D, data = data, family = binomial)
  phat <- predict(zmod1, type = "response")
  #phat <- plogis(-1 + 1 * data$C + .4 * data$C^2 - 1 * data$D)
  data$W1 <- data$Z / phat + (1 - data$Z) / (1 - phat)
  
  phatWR <- predict(glm(Z ~ C, data= data, family = binomial), type = "response")
  #phatWR <- plogis(-1 + 0.5 * data$C)
  data$W2 <- data$Z / phatWR + (1 - data$Z) / (1 - phatWR)
  
  fit1 <- glm(Y ~ Z, data = data, weights = W1, family = binomial())
  fit2 <- glm(Y ~ Z + C + D, data = data, weights = W2, family = binomial())
  
  
  ests <- sapply(list(fit1, fit2), 
                 \(fit) {
                   summary(stdGlm(fit, data, "Z"), contrast = "difference", reference = 0)$est.table[2, 1]
                 })
  
  data.frame(est = ests, 
             type = c("wrong outcome right weights", "right outcome wrong weights"))
  
}


analyze_weibullPH <- function(data) {
  
  
  zmod1 <- glm(Z ~ C + I(C^2) + D, data = data, family = binomial)
  phat <- predict(zmod1, type = "response")
  #phat <- plogis(-1 + 1 * data$C + .4 * data$C^2 - 1 * data$D)
  data$W1 <- data$Z / phat + (1 - data$Z) / (1 - phat)
  
  phatWR <- predict(glm(Z ~ C, data= data, family = binomial), type = "response")
  #phatWR <- plogis(-1 + 0.5 * data$C)
  data$W2 <- data$Z / phatWR + (1 - data$Z) / (1 - phatWR)
  
  fit1 <- flexsurvreg(Surv(Y.obs,  D.ind) ~ Z+C+D, data = data,
            dist = "weibullPH", weights = data$W2, 
            inits = c(2, 1, log(2), 1, 1))
  
  fit2 <- flexsurvreg(Surv(Y.obs,  D.ind) ~ Z, data = data,
                      dist = "weibullPH", weights = data$W1, 
                      inits = c(2, 1, log(2)))
  
  
  data1 <- data0 <- data
  data1$Z <- 1
  data0$Z <- 0
  
  ests <- sapply(list(fit1, fit2), \(fit) {
  S110<-predict(fit, newdata = data1, times=1, type="survival")
  S100<-predict(fit, newdata = data0, times=1, type="survival")
   mean(S110$.pred_survival-S100$.pred_survival)
  })
  
  data.frame(est = ests, 
             type = c("wrong outcome right weights", "right outcome wrong weights"))
  
  
}


analyze_survspline <- function(data) {
  
  
  zmod1 <- glm(Z ~ C + I(C^2) + D, data = data, family = binomial)
  phat <- predict(zmod1, type = "response")
  #phat <- plogis(-1 + 1 * data$C + .4 * data$C^2 - 1 * data$D)
  data$W1 <- data$Z / phat + (1 - data$Z) / (1 - phat)
  
  phatWR <- predict(glm(Z ~ C, data= data, family = binomial), type = "response")
  #phatWR <- plogis(-1 + 0.5 * data$C)
  data$W2 <- data$Z / phatWR + (1 - data$Z) / (1 - phatWR)
  
  fit1 <- flexsurvspline(Surv(Y.obs,  D.ind) ~ Z+C+D, data = data,
    scale = "hazard", k = 3, weights = data$W2) 
  
  fit2 <- flexsurvspline(Surv(Y.obs,  D.ind) ~ Z, data = data,
                         scale = "hazard", k = 3, weights = data$W1) 
  
  data1 <- data0 <- data
  data1$Z <- 1
  data0$Z <- 0
  
  ests <- sapply(list(fit1, fit2), \(fit) {
    S110<-predict(fit, newdata = data1, times=1, type="survival")
    S100<-predict(fit, newdata = data0, times=1, type="survival")
    mean(S110$.pred_survival-S100$.pred_survival)
  })
  
  data.frame(est = ests, 
             type = c("wrong outcome right weights", "right outcome wrong weights"))
  
  
}



analyze_exponential <- function(data) {
  
  
  zmod1 <- glm(Z ~ C + I(C^2) + D, data = data, family = binomial)
  phat <- predict(zmod1, type = "response")
  #phat <- plogis(-1 + 1 * data$C + .4 * data$C^2 - 1 * data$D)
  data$W1 <- data$Z / phat + (1 - data$Z) / (1 - phat)
  
  phatWR <- predict(glm(Z ~ C, data= data, family = binomial), type = "response")
  #phatWR <- plogis(-1 + 0.5 * data$C)
  data$W2 <- data$Z / phatWR + (1 - data$Z) / (1 - phatWR)
  
  fit1 <- survreg(Surv(Ye.obs,E.ind) ~ Z+C+D, data = data,
    dist = "exp", weights = data$W2)
  
  fit2 <- survreg(Surv(Ye.obs,E.ind) ~ Z, data = data,
                  dist = "exp", weights = data$W1)
  
  data1 <- data0 <- data
  data1$Z <- 1
  data0$Z <- 0
  
  ests <- sapply(list(fit1, fit2), \(fit) {
    
    SE110<-1-psurvreg(1, mean = predict(fit, newdata=data1, type="lp"), distribution = "exponential")
    SE100<-1-psurvreg(1, mean = predict(fit, newdata=data0, type="lp"), distribution = "exponential")
    
    coef <- exp(-fit$coefficients[2])
    mean(SE110-SE100)
  })
  
  data.frame(est = ests, 
             type = c("wrong outcome right weights", "right outcome wrong weights"))
  
  
}


analyze_coxph <- function(data) {
  
  
  zmod1 <- glm(Z ~ C + I(C^2) + D, data = data, family = binomial)
  phat <- predict(zmod1, type = "response")
  #phat <- plogis(-1 + 1 * data$C + .4 * data$C^2 - 1 * data$D)
  data$W1 <- data$Z / phat + (1 - data$Z) / (1 - phat)
  
  phatWR <- predict(glm(Z ~ C, data= data, family = binomial), type = "response")
  #phatWR <- plogis(-1 + 0.5 * data$C)
  data$W2 <- data$Z / phatWR + (1 - data$Z) / (1 - phatWR)
  
  fit1 <- coxph(Surv(Y.obs, D.ind) ~ Z, data = data, weights = W1, ties = "breslow")
  fit2 <- coxph(Surv(Y.obs, D.ind) ~ Z + C + D, data = data, weights = W2, ties = "breslow")
  
  
  ests <- sapply(list(fit1, fit2), 
                 \(fit) {
                   summary(stdCoxph(fit, data, "Z", 0:1, t = 1), 
                           contrast = "difference", reference = 0)$est.table[[1]][2, 1]
                 })
  
  data.frame(est = ests, 
             type = c("wrong outcome right weights", "right outcome wrong weights"))
  
  
}
