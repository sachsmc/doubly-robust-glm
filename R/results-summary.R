
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
      facet_wrap(~ analysis, scales = "free") + theme_bw() + 
    theme(legend.position = "bottom") + 
      labs(x = "True value", y = "Bias (mean +/- 1 std. dev.)", color = "Model: ")

}