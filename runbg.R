source("_targets.R")

tar_make_clustermq(workers = 6)



test <- tar_read_raw(tomake[2])
by(test$est, test$type, mean)
mean(test$true_value)
