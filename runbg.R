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

cmps <- sort(grep("compare", tar_manifest()$name, value = TRUE))[c(1, 4, 3, 6, 2, 5)]
res <- NULL
for(i in c(1, 3, 5)) {
    x <- subset(tar_read_raw(cmps[i]), type != "failed")
    x1 <- subset(tar_read_raw(cmps[i + 1]), type != "failed")
    res <- rbind(res, cbind(
      by(x$est.adhoc, x$type, sd),
      by(x$est.funk, x$type, sd),
      by(x1$est.adhoc, x1$type, sd),
      by(x1$est.funk, x1$type, sd)))
  }

tcomp <-res
colnames(tcomp) <- c("GLM.linear", "Funk.linear", "GLM.logit", "Funk.logit")

xtable(tcomp, digits = 3)

ggplot(testlogit, aes(x = est.adhoc, y = est.wtd)) + geom_point() + 
  geom_hline(aes(yintercept = true_value), linetype = 2) + 
  geom_vline(aes(xintercept = true_value), linetype = 2) + 
  facet_wrap(~ type) + theme_bw() + geom_rug()


testlogit[,c(1, 3, 4,5)] |> plot()


chk <- tar_read(simulate_linearodd_ols_weighted_standardized_odd_2)
by(chk$est, chk$type, mean)                                                           
                                                                                            


test <- tar_read(simulate_linear_ols_weighted_2000_2)
results_table(test)


source("_targets.R")
torun <- tar_manifest()$name[c(4,9,10,11)]

tar_delete(names = any_of(!!torun))
tar_make(names = any_of(!!torun))

ares <-do.call(rbind, lapply(torun, \(ne) {
  res <- tar_read_raw(ne)
}))

results_table(ares)
