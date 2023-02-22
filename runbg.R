source("_targets.R")

tar_make()

tar_read(linear_summary)
tar_read(glm_summary)
# 
 library(xtable)
 library(gtsummary)
 library(gt)
library(data.table)
# 
library(ggplot2)
library(ggExtra)

test <- tar_read(linear_summary)


teset <- tar_read(linear_summary)[, bias := sprintf("%.3f (%.2f)", mean_bias, sd_bias)]
teset <- teset[, setting := sprintf("Generation: %s; analysis: %s", setting, analysis)][, .(setting, type, true_value, 
                                                                                            bias, cover_pct.acm, cover_pct.boot)]
gt(teset, groupname_col = "setting") |> 
  fmt_number(columns = true_value, decimals = 2) |> 
  fmt_number(columns = starts_with("cover"), decimals = 1) |> 
  as_latex() |> cat()


 
 
 teset <- tar_read(glm_summary)[, bias := sprintf("%.3f (%.2f)", mean_bias, sd_bias)]
 teset <- teset[, setting := sprintf("Generation: %s; analysis: %s", setting, analysis)][, .(setting, type, true_value, 
                                                                                    bias, cover_pct.acm, cover_pct.boot)]
 gt(teset, groupname_col = "setting") |> 
   fmt_number(columns = true_value, decimals = 2) |> 
   fmt_number(columns = starts_with("cover"), decimals = 1) |> 
   as_latex() |> cat()
 
 
testlin <- tar_read(simulate_linear_linear_compare_2)
testlogit <- tar_read(simulate_logit_binomial_logit_compare_5)

ggplot(rbind(testlin, testlogit[, colnames(testlin)]), aes(x = est.adhoc, y = est.funk)) + geom_point(alpha = .1) + 
  geom_hline(aes(yintercept = true_value), linetype = 2) + 
  geom_vline(aes(xintercept = true_value), linetype = 2) + 
  facet_wrap(setting ~ type, scales = "free") + theme_bw() + geom_rug(alpha = .1)+ 
  labs(x = "GLM", y = "Funk")
ggsave("funk-compare.pdf", width = 6.5, height = 13 / 3)

tcomp <-cbind(by(testlin$est.adhoc  , testlin$type, sd),
by(testlin$est.funk   , testlin$type, sd),
by(testlogit$est.adhoc, testlogit$type, sd),
by(testlogit$est.funk , testlogit$type, sd))

colnames(tcomp) <- c("GLM.linear", "Funk.linear", "GLM.logit", "Funk.logit")

xtable(tcomp, digits = 3)

ggplot(testlogit, aes(x = est.adhoc, y = est.wtd)) + geom_point() + 
  geom_hline(aes(yintercept = true_value), linetype = 2) + 
  geom_vline(aes(xintercept = true_value), linetype = 2) + 
  facet_wrap(~ type) + theme_bw() + geom_rug()


testlogit[,c(1, 3, 4,5)] |> plot()


chk <- tar_read(simulate_linearodd_ols_weighted_standardized_odd_2)
by(chk$est, chk$type, mean)                                                           
                                                                                            