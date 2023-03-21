## variance estimation, 
## start with only the first term E(Y | X = 1)

infunc_confint <- function(ofit, wfit) {
  
  
  data <- ofit$data
  data1 <- data0 <- data
  data1$Z <- 1
  data0$Z <- 0
  
  n <- nrow(data1)
  
  est1f <- mean(predict(ofit, newdata = data1, type = "response"))
  est0f <- mean(predict(ofit, newdata = data0, type = "response"))
  
  # refit using unweighted eeqs
  ofitunwt <- glm(ofit$formula, family = ofit$family, 
                 data = data)

  phat <- predict(wfit, type = "response")
 
  
  est1 <- predict(ofitunwt, newdata = data1, type = "response")
  est0 <- predict(ofitunwt, newdata = data0, type = "response")

  XXw <- model.matrix(wfit)
  XXo <- XXo1 <- XXo0 <- model.matrix(ofitunwt)
  XXo1[,"Z"] <- 1
  XXo0[,"Z"] <- 0
  
  ifweight <- t(vcov(wfit) %*% t(sandwich::estfun(wfit)))
  ifout <- t(vcov(ofit) %*% t(sandwich::estfun(ofitunwt)))
  
  
  eifterms1 <- (data$Z / phat * (data$Y - est1) + (est1 - est1f)) / n
  eifterms0 <- ((1 - data$Z) / (1 - phat) * (data$Y - est0) + (est0 - est0f)) / n
  
  
  hdot <- family(wfit)$mu.eta(predict(wfit, type = "link"))
  
  gdot0 <- family(ofit)$mu.eta(predict(ofitunwt, newdata = data0, type = "link"))
  gdot1 <- family(ofit)$mu.eta(predict(ofitunwt, newdata = data1, type = "link"))
  
  Kterm1 <- (-1/n) * matrix(((data$Z * hdot) / phat^2) * (data$Y - est1), nrow = 1, ncol = n) %*%
    XXw
  Kterm0 <- (1/n) * matrix((((1 - data$Z) * hdot) / (1 - phat)^2) * (data$Y - est0), nrow = 1, ncol = n) %*%
    XXw
  
  Lterm1 <- (1/ n) * matrix(gdot1 * (1 - data$Z/phat), nrow = 1, ncol = n) %*% 
    XXo1
  
  Lterm0 <- (1/ n) * matrix(gdot0 * ((1 - data$Z)/(1 - phat) - 1), nrow = 1, ncol = n) %*% 
    XXo0
  
  fullif1 <- cbind(eifterms1, (ifweight %*% t(Kterm1)), (ifout %*% t(Lterm1)))
  fullif0 <- cbind(eifterms0, (ifweight %*% t(Kterm0)), (ifout %*% t(Lterm0)))
  
  
  estse <- sqrt(sum((rowSums(fullif1) - rowSums(fullif0))^2))
  est <- est1f - est0f
  
  c(lower = est - 1.96 * estse, upper = est + 1.96 * estse)
  
  
}


