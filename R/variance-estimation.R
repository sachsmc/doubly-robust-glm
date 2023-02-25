## variance estimation, 
## start with only the first term E(Y | X = 1)

infunc_confint <- function(ofit, wfit, est) {
  
  ofit <- glm(ofit$formula, family = ofit$family, 
                 data = ofit$data)
  
  data <- ofit$data
  phat <- predict(wfit, type = "response")
  
  data1 <- data0 <- data
  data1$Z <- 1
  data0$Z <- 0
  
  est1 <- predict(ofit, newdata = data1, type = "response")
  est0 <- predict(ofit, newdata = data0, type = "response")
  
  
  eifterms1 <- data$Z / phat * (data$Y - est1) + (est1 - mean(est1))
  ifweight <- t(vcov(wfit) %*% t(sandwich::estfun(wfit)))
  ifout <- t(vcov(ofit) %*% t(sandwich::estfun(ofit)))
  
  eifterms0 <- (1 - data$Z) / (1 - phat) * (data$Y - est0) + (est0 - mean(est0))
  
  XXw <- model.matrix(wfit)
  XXo <- model.matrix(ofit)
  hdot <- XXw * family(wfit)$mu.eta(predict(wfit, type = "link"))
  gdot0 <- model.matrix(terms(ofit), data = data0) * family(ofit)$mu.eta(predict(ofit, newdata = data0, 
                                                                          type = "link"))
  gdot1 <- model.matrix(terms(ofit), data = data1) * family(ofit)$mu.eta(predict(ofit, newdata = data1, 
                                                                          type = "link"))
  
  fullif1 <- cbind(eifterms1, 
                   -diag((data$Z * hdot / phat^2 * (data$Y - est1)) %*% t(ifweight)),
                   diag((gdot1 / phat * (phat - data$Z)) %*% t(ifout)))
  
  fullif0 <- cbind(eifterms0, 
                   diag(((1 - data$Z) * hdot / (1 - phat)^2 * (data$Y - est0)) %*% t(ifweight)),
                   diag((gdot0 / (1 - phat) * (data$Z - phat)) %*% t(ifout)))
  
  
  estse = sd(rowSums(fullif1) - rowSums(fullif0)) / sqrt(nrow(data))
  #est = mean(est1)- mean(est0)
  
  c(lower = est - 1.96 * estse, upper = est + 1.96 * estse)
  
  
}

