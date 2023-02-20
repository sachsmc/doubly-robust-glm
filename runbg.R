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
 