source("_targets.R")

tar_make()



tar_read(linear_summary)
tar_read(glm_summary)

library(xtable)
library(gtsummary)
library(gt)

gt(tar_read(linear_summary)) |> fmt_number(columns = true_value:sd_bias, decimals = 3) |> 
  fmt_number(columns = cover_pct, decimals = 1) |> 
  as_latex() |> cat()
