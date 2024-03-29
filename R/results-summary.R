
results_barplot <- function(data) {

  data <- data.table(data)
  toplot <- data[type != "failed", .(mean_bias = mean(est - true_value), 
           sd_bias = sd(est - true_value)), 
       by = .(setting, analysis, type, true_value)]
  
    ggplot(toplot, 
           aes(x = factor(round(true_value, 2)), y = mean_bias, ymin = mean_bias - sd_bias, 
               ymax = mean_bias + sd_bias, color = type)) + 
      geom_hline(yintercept = 0, linetype = 2) + 
    geom_crossbar(position = position_dodge(width = 1)) + 
      facet_wrap(~ setting + analysis, scales = "free") + theme_bw() + 
    theme(legend.position = "bottom") + 
      labs(x = "True value", y = "Bias (mean +/- 1 std. dev.)", color = "Model: ")

}


results_table <- function(data) {
  
  data <- data.table(data)
  data[, in_ci.acm := true_value >= lowerCL.acm & true_value <= upperCL.acm]
  data[, in_ci.boot := true_value >= lowerCL.boot & true_value <= upperCL.boot]
  data[, in_ci.infl := true_value >= lowerCL.infl & true_value <= upperCL.infl]
  data[type != "failed", .(mean_bias = round(mean(est - true_value),4), 
                                     sd_bias = sd(est - true_value), 
                           percent_bias = 100 * mean(est - true_value) / mean(true_value), 
                           infl_se = mean(se.infl), 
                           cover_pct.acm = 100 * mean(in_ci.acm), 
                           cover_pct.boot = 100 * mean(in_ci.boot),
                           cover_pct.infl = 100 * mean(in_ci.infl)
                           ), 
                 by = .(setting, analysis, type, true_value)]
  
  
}

stddev_table <- function(data) {
  
  data <- data.table(data)
  res <- data[, .(sd.glm = sd(est.adhoc), sd.funk = sd(est.funk)), keyby = .(sampsize, setting, type)] |> 
    dcast(sampsize + type ~ setting, value.var = c("sd.glm", "sd.funk"))
  res[, .(sampsize, type, sd.glm_linear, sd.funk_linear, sd.glm_logit_binomial, sd.funk_logit_binomial)]
  
  
}



eif_table <- function(data) {
  
  data <- data.table(data)
  data[type != "failed", .(emp.se = sd(est), 
                           eif.se = mean(se.eif), 
                           infl.se = mean(se.infl)), 
       by = .(setting, analysis, type)]
  
  
}