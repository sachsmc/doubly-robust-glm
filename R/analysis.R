#' Run the simulation according to the settings
#' 
#' @param generation linear, log_poisson, log_binomial, logit_binomial, or survival
#' @param analysis ols, ols_weighted, ols_weighted_standardized, 
#'                  glm_weighted, glm_weighted_standardized, survspline, weibullPH, exponential, or coxph
#' @param coefZ numeric value for the regression coefficient
#' @param B number of simulation replicates
#' @param n sample size in each generation
#'                  


run_simulation <- function(generation, analysis, coefZ, n = 2000) {
  
  ## number of simulation replicates
  
  generate_data <- get(paste0("generate_", generation))
  analyze_data <- get(paste0("analyze_", analysis))
  true_values <- get(paste0("truev_", generation))
  
  truev <- true_values(coefZ)
  
  if(generation == "survival") {
    truev <- if(analysis == "exponential") truev[2] else truev[1]
  }
  
  #results <- NULL
  #for(i in 1:B) {
  #  
  #  datin <- generate_data(coefZ = coefZ) 
  #  res <- tryCatch(analyze_data(datin), error = function(e) data.frame(est = NA, lowerCL = NA, upperCL = NA, type = "failed"))
  #  results <- rbind(results, res)
  #  
  #}
  
  results <- do.call(rbind, mclapply(1:B, \(i) {
    
    datin <- generate_data(n = n, coefZ = coefZ) 
    tryCatch(analyze_data(datin), error = function(e) data.frame(est = NA, lowerCL.acm = NA, 
                                                                 upperCL.acm = NA, lowerCL.boot = NA, upperCL.boot = NA, 
                                                                 lowerCL.infl = NA, upperCL.infl = NA, 
                                                                 se.eif = NA,
                                                                 se.infl = NA,
                                                                 type = "failed"))
    
  }, mc.cores = cores))
  
  results$setting <- generation
  results$analysis <- analysis
  results$true_value <- truev
  results$coefZ <- coefZ
  results$sampsize <- n
  results
  
}


analyze_ols <- function(data) {
  
  fit1 <- lm(Y ~ Z + C + I(C^2) + D, data = data)
  fit2 <- lm(Y ~ Z + C, data = data)
  data.frame(est = c(coefficients(fit1)[2], coefficients(fit2)[2]), 
             lowerCL.acm = sapply(list(fit1, fit2),\(x) confint(x)[2, 1]), 
             upperCL.acm = sapply(list(fit1, fit2),\(x) confint(x)[2, 2]),
             lowerCL.boot = NA, upperCL.boot = NA,
             lowerCL.infl = NA, upperCL.infl = NA,
             type = c("right outcome no weights", "wrong outcome no weights"))
  
}

analyze_ols_weighted <- function(data, nboot = boots) {
  
  
  zmod1 <- glm(Z ~ C + I(C^2) + D, data = data, family = binomial)
  phat <- predict(zmod1, type = "response")
  #phat <- plogis(-1 + 1 * data$C + .4 * data$C^2 - 1 * data$D)
  data$W1 <- data$Z / phat + (1 - data$Z) / (1 - phat)
  
  zmod2 <- glm(Z ~ C, data= data, family = binomial)
  phatWR <- predict(zmod2, type = "response")
  #phatWR <- plogis(-1 + 0.5 * data$C)
  data$W2 <- data$Z / phatWR + (1 - data$Z) / (1 - phatWR)
  
  fit10 <- glm(Y ~ Z + C, data = data, weights = W1)
  fit20 <- glm(Y ~ Z + C + I(C^2) + D, data = data, weights = W2)
  fit30 <- glm(Y ~ Z + C, data = data, weights = W2)
  fit40 <- glm(Y ~ Z + C + I(C^2) + D, data = data, weights = W1)
  
  bootb <- matrix(NA, ncol = 4, nrow = nboot)
  for(i in 1:nrow(bootb)) {
    datab <- data[sample(1:nrow(data), nrow(data), replace = TRUE),]
    
    zmod1 <- glm(Z ~ C + I(C^2) + D, data = datab, family = binomial)
    phat <- predict(zmod1, type = "response")
    #phat <- plogis(-1 + 1 * data$C + .4 * data$C^2 - 1 * data$D)
    datab$W1 <- datab$Z / phat + (1 - datab$Z) / (1 - phat)
    
    phatWR <- predict(glm(Z ~ C, data= datab, family = binomial), type = "response")
    #phatWR <- plogis(-1 + 0.5 * data$C)
    datab$W2 <- datab$Z / phatWR + (1 - datab$Z) / (1 - phatWR)
    
    fit1 <- lm(Y ~ Z + C, data = datab, weights = W1)
    fit2 <- lm(Y ~ Z + C + I(C^2) + D, data = datab, weights = W2)
    fit3 <- lm(Y ~ Z + C, data = datab, weights = W2)
    fit4 <- glm(Y ~ Z + C + I(C^2) + D, data = datab, weights = W1)
    
    bootb[i, ] <- sapply(list(fit1, fit2, fit3, fit4), \(x) x$coefficients[2])
  }
  
  bootci <- apply(bootb, 2, quantile, c(.025, .975))
  
  infcis <- sapply(list(list(fit10, zmod1), list(fit20, zmod2), list(fit30, zmod2), list(fit40, zmod1)), 
  \(lp) infunc_confint(lp[[1]], lp[[2]]))
  
  
  data.frame(est = c(coefficients(fit10)[2], coefficients(fit20)[2], coefficients(fit30)[2], coefficients(fit40)[2]), 
             lowerCL.acm = sapply(list(fit10, fit20, fit30, fit40),\(x) confint(x)[2, 1]), 
             upperCL.acm = sapply(list(fit10, fit20, fit30, fit40),\(x) confint(x)[2, 2]),
             lowerCL.boot = bootci[1, ], 
             upperCL.boot = bootci[2, ],
             lowerCL.infl = infcis[1, ],
             upperCL.infl = infcis[2, ],
             se.eif = infcis[3,],
             se.infl = infcis[4,],
             type = c("wrong outcome right weights", "right outcome wrong weights", "wrong both", "right both"))
  
}


stdGlm2 <- function(data, ofit, wfit, nboot = boots, ostart = NULL, noweights = TRUE) {
  
  rewfit0 <- glm(wfit$formula, family = wfit$family, data = data)
  phat <- predict(rewfit0, type = "response")
  data$ww <- data$Z / phat + (1 - data$Z) / (1 - phat)
  reofit0 <- glm(ofit$formula, family = ofit$family, weights = ww, 
                 data = data, start = ostart)
  
  if(noweights) {
    reofit0$prior.weights <- rep(1, length(reofit0$prior.weights))
  }
  mainest <- summary(stdGlm(reofit0, data, "Z"), contrast = "difference", reference = 0)$est.table[2, c(1, 3, 4)]
  
  bootests <- rep(NA, nboot) 
  for(i in 1:nboot) {
    
    thisstart <- if(!is.null(ostart)) ostart else reofit0$coefficients
    
    bdata <- data[sample(1:nrow(data), nrow(data), replace = TRUE),]
    rewfit <- glm(wfit$formula, family = wfit$family, data = bdata, 
                  start = rewfit0$coefficients, x = FALSE, y = FALSE)
    phat <- predict(rewfit, type = "response")
    bdata$ww <- bdata$Z / phat + (1 - bdata$Z) / (1 - phat)
    reofit <- glm(ofit$formula, family = ofit$family, weights = ww, data = bdata, start = thisstart, 
                  x = FALSE, y = FALSE)
    
    reofit$prior.weights <- rep(1, length(reofit$prior.weights))
    bootests[i] <- summary(stdGlm(reofit, bdata, "Z"), contrast = "difference", reference = 0)$est.table[2, c(1)]
    
  }
  
  bootci <- quantile(bootests, c(.025, .975))
  inflci <- infunc_confint(reofit0, rewfit0, start = ostart)
  
  res <- c(mainest, bootci, inflci)
  names(res) <- c("est", "lowerCL.acm", "upperCL.acm", "lowerCL.boot", "upperCL.boot", 
                  "lowerCL.infl", "upperCL.infl", "se.eif", "se.infl")
  
  as.data.frame(as.list(res))
  
}

analyze_ols_weighted_standardized <- function(data) {
  
  res <- rbind.data.frame(
    stdGlm2(data, glm(Y ~ Z + C, data = data, family = "gaussian"), 
                  glm(Z ~ C + I(C^2) + D, data = data, family = binomial)), 
    stdGlm2(data, glm(Y ~ Z + C + I(C^2) + D, data = data, family = "gaussian"), 
                  glm(Z ~ C, data = data, family = binomial)), 
    stdGlm2(data, glm(Y ~ Z + C, data = data, family = "gaussian"), 
                  glm(Z ~ C, data = data, family = binomial)), 
    stdGlm2(data, glm(Y ~ Z + C + I(C^2) + D, data = data, family = "gaussian"), 
            glm(Z ~ C + I(C^2) + D, data = data, family = binomial))
  )
  res$type <- c("wrong outcome right weights", "right outcome wrong weights", "wrong both", "right both")
  res
  
}

analyze_ols_weighted_standardized_odd <- function(data) {
  
  res <- rbind.data.frame(
    stdGlm2(data, glm(Y ~ Z + C, data = data, family = "gaussian"), 
            glm(Z ~ D, data = data, family = binomial), nboot = 1)
  )
  res$type <- c("residual confounding test")
  res
  
}




analyze_poisson_weighted_standardized <- function(data) {
  
  res <- rbind.data.frame(
    stdGlm2(data, glm(Y ~ Z + C, data = data, family = "poisson"), 
            glm(Z ~ C + I(C^2) + D, data = data, family = binomial)), 
    stdGlm2(data, glm(Y ~ Z + C + I(C^2) + D, data = data, family = "poisson"), 
            glm(Z ~ C, data = data, family = binomial)), 
    stdGlm2(data, glm(Y ~ Z + C, data = data, family = "poisson"), 
            glm(Z ~ C, data = data, family = binomial)),
    stdGlm2(data, glm(Y ~ Z + C + I(C^2)+ D, data = data, family = "poisson"), 
            glm(Z ~ C + I(C^2) + D, data = data, family = binomial))
  )
  res$type <- c("wrong outcome right weights", "right outcome wrong weights", "wrong both", "right both")
  res
  
}


analyze_inverse_gaussian_weighted_standardized <- function(data) {
  
  res <- rbind.data.frame(
    stdGlm2(data, glm(Y ~ Z + C, data = data, family = "inverse.gaussian"), 
            glm(Z ~ C + I(C^2) + D, data = data, family = binomial)), 
    stdGlm2(data, glm(Y ~ Z + C + I(C^2) + D, data = data, family = "inverse.gaussian"), 
            glm(Z ~ C, data = data, family = binomial)), 
    stdGlm2(data, glm(Y ~ Z + C, data = data, family = "inverse.gaussian"), 
            glm(Z ~ C, data = data, family = binomial)), 
    stdGlm2(data, glm(Y ~ Z + C + I(C^2) + D, data = data, family = "inverse.gaussian"), 
            glm(Z ~ C + I(C^2) + D, data = data, family = binomial))
  )
  res$type <- c("wrong outcome right weights", "right outcome wrong weights", "wrong both", "right both")
  res
  
}




analyze_log_binomial_weighted_standardized <- function(data) {
  
  lmeanY <- log(mean(data$Y))
  
  res <- rbind.data.frame(
    stdGlm2(data, glm(Y ~ Z + C, data = data, family = binomial(link = "log"), start = c(lmeanY, 0, 0)), 
            glm(Z ~ C + I(C^2) + D, data = data, family = binomial), ostart = c(lmeanY, 0, 0), nboot = boots), 
    stdGlm2(data, glm(Y ~ Z + C + I(C^2) + D, data = data, family = binomial(link = "log"), start = c(lmeanY, 0, 0, 0, 0)), 
            glm(Z ~ C, data = data, family = binomial), ostart = c(lmeanY, 0, 0, 0, 0), nboot = boots), 
    stdGlm2(data, glm(Y ~ Z + C, data = data, family = binomial(link = "log"), start = c(lmeanY, 0, 0)), 
            glm(Z ~ C, data = data, family = binomial), ostart = c(lmeanY, 0, 0), nboot = boots),
    stdGlm2(data, glm(Y ~ Z + C + I(C^2) + D, data = data, family = binomial(link = "log"), start = c(lmeanY, 0, 0, 0, 0)), 
            glm(Z ~ C + I(C^2) + D, data = data, family = binomial), ostart = c(lmeanY, 0, 0, 0, 0), nboot = boots) 
  )
  res$type <- c("wrong outcome right weights", "right outcome wrong weights", "wrong both", "right both")
  res
  
   
}


analyze_identity_binomial_weighted_standardized <- function(data) {
  
  
  lmeanY <- mean(data$Y)
  
  res <- rbind.data.frame(
    stdGlm2(data, glm(Y ~ Z + C, data = data, family = binomial(link = "identity")), 
            glm(Z ~ C + I(C^2) + D, data = data, family = binomial), ostart = c(lmeanY, 0.1, 0)), 
    stdGlm2(data, glm(Y ~ Z + C + I(C^2) + D, data = data, family = binomial(link = "identity")), 
            glm(Z ~ C, data = data, family = binomial), ostart = c(lmeanY, 0.1, 0, 0, 0)), 
    stdGlm2(data, glm(Y ~ Z + C, data = data, family = binomial(link = "identity")), 
            glm(Z ~ C, data = data, family = binomial), ostart = c(lmeanY, 0.1, 0)), 
    stdGlm2(data, glm(Y ~ Z + C + I(C^2) + D, data = data, family = binomial(link = "identity")), 
            glm(Z ~ C + I(C^2) + D, data = data, family = binomial), ostart = c(lmeanY, 0.1, 0)) 
  )
  res$type <- c("wrong outcome right weights", "right outcome wrong weights", "wrong both", "right both")
  res
  
}



analyze_logit_binomial_weighted_standardized <- function(data) {
  
  res <- rbind.data.frame(
    stdGlm2(data, glm(Y ~ Z + C, data = data, family = binomial), 
            glm(Z ~ C + I(C^2) + D, data = data, family = binomial)), 
    stdGlm2(data, glm(Y ~ Z + C + I(C^2) + D, data = data, family = binomial), 
            glm(Z ~ C, data = data, family = binomial)), 
    stdGlm2(data, glm(Y ~ Z + C, data = data, family = binomial), 
            glm(Z ~ C, data = data, family = binomial)), 
    stdGlm2(data, glm(Y ~ Z + C + I(C^2) + D, data = data, family = binomial), 
            glm(Z ~ C + I(C^2) + D, data = data, family = binomial))
  )
  res$type <- c("wrong outcome right weights", "right outcome wrong weights", "wrong both", "right both")
  res
  
}


analyze_probit_binomial_weighted_standardized <- function(data) {
  
  
  res <- rbind.data.frame(
    stdGlm2(data, glm(Y ~ Z + C, data = data, family = binomial(link = "probit")), 
            glm(Z ~ C + I(C^2) + D, data = data, family = binomial)), 
    stdGlm2(data, glm(Y ~ Z + C + I(C^2) + D, data = data, family = binomial(link = "probit")), 
            glm(Z ~ C, data = data, family = binomial)), 
    stdGlm2(data, glm(Y ~ Z + C, data = data, family = binomial(link = "probit")), 
            glm(Z ~ C, data = data, family = binomial)), 
    stdGlm2(data, glm(Y ~ Z + C + I(C^2) + D, data = data, family = binomial(link = "probit")), 
            glm(Z ~ C + I(C^2) + D, data = data, family = binomial)), 
  )
  res$type <- c("wrong outcome right weights", "right outcome wrong weights", "wrong both", "right both")
  res
  
  
}

funk_estimate <- function(data, ofit0, ofit1, wfit, weight = FALSE) {
  
  
  
  data0 <- data1 <- data
  data0$Z <- 0
  data1$Z <- 1
  ps <- predict(wfit, type = "response")
  
  if(weight) {

    data$ww <- data$Z / ps + (1 - data$Z) / (1 - ps)
    ofit1 <- glm(ofit1$formula, family = ofit1$family, weights = ww, 
                   data = data[data$Z == 1,])
    ofit0 <- glm(ofit0$formula, family = ofit0$family, weights = ww, 
                 data = data[data$Z == 0,])
    
  }
  
  dr1 <- data$Y * (data$Z == 1) / ps - predict(ofit1, newdata = data1, type = "response") * ((data$Z == 1) - ps) / ps
  dr0 <- data$Y * (data$Z == 0) / (1 - ps) + predict(ofit0, newdata = data0, type = "response") * ((data$Z == 1) - ps) / (1 - ps)
  
  mean(dr1) - mean(dr0)
  
}

analyze_linear_compare <- function(data) {
  
  
  wofit <- glm(Y ~ Z + C, data = data, family = "gaussian")
  cwfit <- glm(Z ~ C + I(C^2) + D, data = data, family = binomial)
  cofit <- glm(Y ~ Z + C + I(C^2) + D, data = data, family = "gaussian")
  wwfit <- glm(Z ~ C, data = data, family = binomial)
  
  # wofit0 <- glm(Y ~ C, data = data[data$Z == 0,], family = "gaussian")
  # wofit1 <- glm(Y ~ C, data = data[data$Z == 1,], family = "gaussian")
  # cofit0 <- glm(Y ~ C + I(C^2) + D, data = data[data$Z == 0,], family = "gaussian")
  # cofit1 <- glm(Y ~ C + I(C^2) + D, data = data[data$Z == 1,], family = "gaussian")
  
  res <- rbind.data.frame(
    stdGlm2(data, wofit, cwfit, nboot = 1), 
    stdGlm2(data, cofit, wwfit, nboot = 1), 
    stdGlm2(data, wofit, wwfit, nboot = 1), 
    stdGlm2(data, cofit, cwfit, nboot = 1)
  )
  res$type <- c("wrong outcome right weights", "right outcome wrong weights", "wrong both", "right both")
  res$lowerCL.acm <- res$upperCL.acm <- 
    res$lowerCL.boot <- res$upperCL.boot <- 
    res$lowerCL.infl <- res$upperCL.infl <- res$se.infl <-
    NULL
  colnames(res)[1] <- "est.adhoc"
  
  res2 <- rbind.data.frame(
    stdGlm2(data, wofit, cwfit, nboot = 1, noweights = FALSE), 
    stdGlm2(data, cofit, wwfit, nboot = 1, noweights = FALSE), 
    stdGlm2(data, wofit, wwfit, nboot = 1, noweights = FALSE), 
    stdGlm2(data, cofit, cwfit, nboot = 1, noweights = FALSE)
  )
  res2$type <- c("wrong outcome right weights", "right outcome wrong weights", "wrong both", "right both")
  res2$lowerCL.acm <- res2$upperCL.acm <- 
    res2$lowerCL.boot <- res2$upperCL.boot <-
    res$lowerCL.infl <- res$upperCL.infl <-  res$se.infl <-NULL
  colnames(res2)[1] <- "est.wtd"
  
  res$est.wtd <- res2$est.wtd
  
  res$est.funk <- c(
    funk_estimate(data, wofit, wofit, cwfit), 
    funk_estimate(data, cofit, cofit, wwfit), 
    funk_estimate(data, wofit, wofit, wwfit), 
    funk_estimate(data, cofit, cofit, cwfit))
  res
  
}


analyze_logit_compare <- function(data) {
  
  tryCatch({
  wofit <- glm(Y ~ Z + C, data = data, family = "binomial")
  cwfit <- glm(Z ~ C + I(C^2) + D, data = data, family = binomial)
  cofit <- glm(Y ~ Z + C + I(C^2) + D, data = data, family = "binomial")
  wwfit <- glm(Z ~ C, data = data, family = binomial)
  
##  wofit0 <- glm(Y ~ C, data = data[data$Z == 0,], family = "binomial")
#  wofit1 <- glm(Y ~ C, data = data[data$Z == 1,], family = "binomial")
#  cofit0 <- glm(Y ~ C + I(C^2) + D, data = data[data$Z == 0,], family = "binomial")
#  cofit1 <- glm(Y ~ C + I(C^2) + D, data = data[data$Z == 1,], family = "binomial")
  
  
  res <- rbind.data.frame(
    stdGlm2(data, wofit, cwfit, nboot = 1), 
    stdGlm2(data, cofit, wwfit, nboot = 1), 
    stdGlm2(data, wofit, wwfit, nboot = 1), 
    stdGlm2(data, cofit, cwfit, nboot = 1)
  )
  res$type <- c("wrong outcome right weights", "right outcome wrong weights", "wrong both", "right both")
  res$lowerCL.acm <- res$upperCL.acm <- res$lowerCL.boot <- res$upperCL.boot <- 
    res$lowerCL.infl <- res$upperCL.infl <-  res$se.infl <-NULL
  colnames(res)[1] <- "est.adhoc"
  
  
  res2 <- rbind.data.frame(
    stdGlm2(data, wofit, cwfit, nboot = 1, noweights = FALSE), 
    stdGlm2(data, cofit, wwfit, nboot = 1, noweights = FALSE), 
    stdGlm2(data, wofit, wwfit, nboot = 1, noweights = FALSE),
    stdGlm2(dtaa, cofit, cwfit, nboot = 1, noweights = FALSE)
  )
  res2$type <- c("wrong outcome right weights", "right outcome wrong weights", "wrong both", "right both")
  res2$lowerCL.acm <- res2$upperCL.acm <- 
    res2$lowerCL.boot <- res2$upperCL.boot <- 
    res2$lowerCL.infl <- res$upperCL.infl <-  res$se.infl <-NULL
  colnames(res2)[1] <- "est.wtd"
  
  res$est.wtd <- res2$est.wtd
  
  res$est.funk <- c(funk_estimate(data, wofit, wofit, cwfit), 
                    funk_estimate(data, cofit, cofit, wwfit), 
                    funk_estimate(data, wofit, wofit, wwfit), 
                    funk_estimate(data, cofit, cofit, cwfit))
  
  res$est.funk.wtd <- c(
    funk_estimate(data, wofit, wofit, cwfit, weight = TRUE), 
    funk_estimate(data, cofit, cofit, wwfit, weight = TRUE), 
    funk_estimate(data, wofit, wofit, wwfit, weight = TRUE), 
    funk_estimate(data, cofit, cofit, cwfit, weight = TRUE))
  res
  }, error = function(e) {
    data.frame(est.adhoc = NA, type = "failed", est.wtd = NA, est.funk = NA, est.funk.wtd = NA)
  })
  
}