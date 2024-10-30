set.seed(1)
y <- rexp(100, 1)
delta <- rbinom(100, size = 1, prob = 0.5)
X <- data.frame(rnorm(100), rnorm(100))
landmark_times <- quantile(y, probs = c(0.25, 0.5, 0.75))
approx_times <- sort(c(unique(y), landmark_times))

################################
### no crossfit, no sample split
################################
f_hat <- list(f_hat = matrix(runif(300), nrow = 100, ncol = length(landmark_times)))
fs_hat <- list(fs_hat = matrix(runif(300), nrow = 100, ncol = length(landmark_times)))
S_hat <- list(S_hat = matrix(rep(seq(1, 0.1, length.out = length(approx_times)), 100),
                             nrow = 100,
                             ncol = length(approx_times),
                             byrow = TRUE))
G_hat <- list(G_hat = matrix(rep(seq(1, 0.1, length.out = length(approx_times)), 100),
                             nrow = 100,
                             ncol = length(approx_times),
                             byrow = TRUE))
folds <- rep(1, 100)
ss_folds <- rep(1, 100)

# accuracy
output <- vim_accuracy(time = y,
                       event = delta,
                       approx_times = approx_times,
                       landmark_times = landmark_times,
                       f_hat = f_hat,
                       fs_hat = fs_hat,
                       S_hat = S_hat,
                       G_hat = G_hat,
                       cf_folds = folds,
                       ss_folds = ss_folds,
                       sample_split = FALSE)
test_that("vim_accuracy(). no xfit, no sample split", {
  expect_equal(dim(output)[1], 3)
  expect_equal(dim(output)[2], 9)
  expect_equal(names(output), c("landmark_time", "est", "var_est", "cil", "ciu",
                                "cil_1sided", "p", "large_predictiveness", "small_predictiveness"))
  expect_equal(sum(is.na(output)), 3)
})

# AUC
output <- vim_AUC(time = y,
                  event = delta,
                  approx_times = approx_times,
                  landmark_times = landmark_times,
                  f_hat = f_hat,
                  fs_hat = fs_hat,
                  S_hat = S_hat,
                  G_hat = G_hat,
                  cf_folds = folds,
                  ss_folds = ss_folds,
                  sample_split = FALSE)
test_that("vim_AUC(). no xfit, no sample split", {
  expect_equal(dim(output)[1], 3)
  expect_equal(dim(output)[2], 9)
  expect_equal(names(output), c("landmark_time", "est", "var_est", "cil", "ciu",
                                "cil_1sided", "p", "large_predictiveness", "small_predictiveness"))
  expect_equal(sum(is.na(output)), 3)
})

# Brier score
output <- vim_brier(time = y,
                    event = delta,
                    approx_times = approx_times,
                    landmark_times = landmark_times,
                    f_hat = f_hat,
                    fs_hat = fs_hat,
                    S_hat = S_hat,
                    G_hat = G_hat,
                    cf_folds = folds,
                    ss_folds = ss_folds,
                    sample_split = FALSE)
test_that("vim_brier(). no xfit, no sample split", {
  expect_equal(dim(output)[1], 3)
  expect_equal(dim(output)[2], 9)
  expect_equal(names(output), c("landmark_time", "est", "var_est", "cil", "ciu",
                                "cil_1sided", "p", "large_predictiveness", "small_predictiveness"))
  expect_equal(sum(is.na(output)), 3)
})

# R-squared
output <- vim_rsquared(time = y,
                       event = delta,
                       approx_times = approx_times,
                       landmark_times = landmark_times,
                       f_hat = f_hat,
                       fs_hat = fs_hat,
                       S_hat = S_hat,
                       G_hat = G_hat,
                       cf_folds = folds,
                       ss_folds = ss_folds,
                       sample_split = FALSE)
test_that("vim_rsquared(). no xfit, no sample split", {
  expect_equal(dim(output)[1], 3)
  expect_equal(dim(output)[2], 9)
  expect_equal(names(output), c("landmark_time", "est", "var_est", "cil", "ciu",
                                "cil_1sided", "p", "large_predictiveness", "small_predictiveness"))
  expect_equal(sum(is.na(output)), 3)
})

# RMST MSE
output <- vim_survival_time_mse(time = y,
                                event = delta,
                                approx_times = approx_times,
                                restriction_time = landmark_times[3],
                                f_hat = lapply(f_hat, function(x) x[,1]),
                                fs_hat = lapply(fs_hat, function(x) x[,1]),
                                S_hat = S_hat,
                                G_hat = G_hat,
                                cf_folds = folds,
                                ss_folds = ss_folds,
                                sample_split = FALSE)
test_that("vim_rmst_mse(). no xfit, no sample split", {
  expect_equal(dim(output)[1], 1)
  expect_equal(dim(output)[2], 9)
  expect_equal(names(output), c("restriction_time", "est", "var_est", "cil", "ciu",
                                "cil_1sided", "p", "large_predictiveness", "small_predictiveness"))
  expect_equal(sum(is.na(output)), 1)
})

# C-index
output <- vim_cindex(time = y,
                     event = delta,
                     approx_times = approx_times,
                     restriction_time = landmark_times[3],
                     f_hat = lapply(f_hat, function(x) x[,1]),
                     fs_hat = lapply(fs_hat, function(x) x[,1]),
                     S_hat = S_hat,
                     G_hat = G_hat,
                     cf_folds = folds,
                     ss_folds = ss_folds,
                     sample_split = FALSE)
test_that("vim_cindex(). no xfit, no sample split", {
  expect_equal(dim(output)[1], 1)
  expect_equal(dim(output)[2], 9)
  expect_equal(names(output), c("restriction_time", "est", "var_est", "cil", "ciu",
                                "cil_1sided", "p", "large_predictiveness", "small_predictiveness"))
  expect_equal(sum(is.na(output)), 1)
})


#############################
### no crossfit, sample split
#############################
f_hat <- list(f_hat_1 = matrix(runif(150), nrow = 50, ncol = length(landmark_times)),
              f_hat_2 = matrix(runif(150), nrow = 50, ncol = length(landmark_times)))
fs_hat <- list(fs_hat_1 = matrix(runif(150), nrow = 50, ncol = length(landmark_times)),
               fs_hat_2 = matrix(runif(150), nrow = 50, ncol = length(landmark_times)))
S_hat <- list(S_hat_1 = matrix(rep(seq(1, 0.1, length.out = length(approx_times)), 50),
                               nrow = 50,
                               ncol = length(approx_times),
                               byrow = TRUE),
              S_hat_2 = matrix(rep(seq(1, 0.1, length.out = length(approx_times)), 50),
                               nrow = 50,
                               ncol = length(approx_times),
                               byrow = TRUE))
G_hat <- list(G_hat_1 = matrix(rep(seq(1, 0.1, length.out = length(approx_times)), 50),
                               nrow = 50,
                               ncol = length(approx_times),
                               byrow = TRUE),
              G_hat_2 = matrix(rep(seq(1, 0.1, length.out = length(approx_times)), 50),
                               nrow = 50,
                               ncol = length(approx_times),
                               byrow = TRUE))
folds <- c(rep(1, 50), rep(2, 50))
ss_folds <- c(rep(1, 50), rep(0, 50))

# accuracy
output <- vim_accuracy(time = y,
                       event = delta,
                       approx_times = approx_times,
                       landmark_times = landmark_times,
                       f_hat = f_hat,
                       fs_hat = fs_hat,
                       S_hat = S_hat,
                       G_hat = G_hat,
                       cf_folds = folds,
                       ss_folds = ss_folds,
                       sample_split = TRUE)
test_that("vim_accuracy(). no xfit, sample split", {
  expect_equal(dim(output)[1], 3)
  expect_equal(dim(output)[2], 9)
  expect_equal(names(output), c("landmark_time", "est", "var_est", "cil", "ciu",
                                "cil_1sided", "p", "large_predictiveness", "small_predictiveness"))
  expect_equal(sum(is.na(output)), 0)
})

# AUC
output <- vim_AUC(time = y,
                  event = delta,
                  approx_times = approx_times,
                  landmark_times = landmark_times,
                  f_hat = f_hat,
                  fs_hat = fs_hat,
                  S_hat = S_hat,
                  G_hat = G_hat,
                  cf_folds = folds,
                  ss_folds = ss_folds,
                  sample_split = TRUE)
test_that("vim_AUC(). no xfit, sample split", {
  expect_equal(dim(output)[1], 3)
  expect_equal(dim(output)[2], 9)
  expect_equal(names(output), c("landmark_time", "est", "var_est", "cil", "ciu",
                                "cil_1sided", "p", "large_predictiveness", "small_predictiveness"))
  expect_equal(sum(is.na(output)), 0)
})

# Brier score
output <- vim_brier(time = y,
                    event = delta,
                    approx_times = approx_times,
                    landmark_times = landmark_times,
                    f_hat = f_hat,
                    fs_hat = fs_hat,
                    S_hat = S_hat,
                    G_hat = G_hat,
                    cf_folds = folds,
                    ss_folds = ss_folds,
                    sample_split = TRUE)
test_that("vim_brier(). no xfit, sample split", {
  expect_equal(dim(output)[1], 3)
  expect_equal(dim(output)[2], 9)
  expect_equal(names(output), c("landmark_time", "est", "var_est", "cil", "ciu",
                                "cil_1sided", "p", "large_predictiveness", "small_predictiveness"))
  expect_equal(sum(is.na(output)), 0)
})

# R-squared
output <- vim_rsquared(time = y,
                       event = delta,
                       approx_times = approx_times,
                       landmark_times = landmark_times,
                       f_hat = f_hat,
                       fs_hat = fs_hat,
                       S_hat = S_hat,
                       G_hat = G_hat,
                       cf_folds = folds,
                       ss_folds = ss_folds,
                       sample_split = TRUE)
test_that("vim_rsquared(). no xfit, sample split", {
  expect_equal(dim(output)[1], 3)
  expect_equal(dim(output)[2], 9)
  expect_equal(names(output), c("landmark_time", "est", "var_est", "cil", "ciu",
                                "cil_1sided", "p", "large_predictiveness", "small_predictiveness"))
  expect_equal(sum(is.na(output)), 0)
})

# RMST MSE
output <- vim_survival_time_mse(time = y,
                                event = delta,
                                approx_times = approx_times,
                                restriction_time = landmark_times[3],
                                f_hat = lapply(f_hat, function(x) x[,1]),
                                fs_hat = lapply(fs_hat, function(x) x[,1]),
                                S_hat = S_hat,
                                G_hat = G_hat,
                                cf_folds = folds,
                                ss_folds = ss_folds,
                                sample_split = TRUE)
test_that("vim_rmst_mse(). no xfit, sample split", {
  expect_equal(dim(output)[1], 1)
  expect_equal(dim(output)[2], 9)
  expect_equal(names(output), c("restriction_time", "est", "var_est", "cil", "ciu",
                                "cil_1sided", "p", "large_predictiveness", "small_predictiveness"))
  expect_equal(sum(is.na(output)), 0)
})

# C-index
output <- vim_cindex(time = y,
                     event = delta,
                     approx_times = approx_times,
                     restriction_time = landmark_times[3],
                     f_hat = lapply(f_hat, function(x) x[,1]),
                     fs_hat = lapply(fs_hat, function(x) x[,1]),
                     S_hat = S_hat,
                     G_hat = G_hat,
                     cf_folds = folds,
                     ss_folds = ss_folds,
                     sample_split = TRUE)
test_that("vim_cindex(). no xfit, sample split", {
  expect_equal(dim(output)[1], 1)
  expect_equal(dim(output)[2], 9)
  expect_equal(names(output), c("restriction_time", "est", "var_est", "cil", "ciu",
                                "cil_1sided", "p", "large_predictiveness", "small_predictiveness"))
  expect_equal(sum(is.na(output)), 0)
})

#####################
### main VIM function
#####################
output <- vim(type = "AUC",
              time = y,
              event = delta,
              X = X,
              landmark_times = landmark_times,
              large_feature_vector = 1:2,
              small_feature_vector = 1,
              conditional_surv_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                        bin_size = 0.1,
                                                        V = 2),
              large_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                    V = 2),
              small_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                    V = 2),
              cf_fold_num = 2,
              sample_split = TRUE,
              scale_est = TRUE)
test_that("vim(). AUC, xfit, sample split", {
  expect_equal(dim(output$result)[1], 3)
  expect_equal(dim(output$result)[2], 12)
  expect_equal(names(output$result), c("landmark_time", "est", "var_est", "cil", "ciu",
                                "cil_1sided", "p", "large_predictiveness", "small_predictiveness",
                                "vim", "large_feature_vector", "small_feature_vector"))
  expect_equal(sum(is.na(output$result)), 0)
  expect_equal(names(output$folds), c("cf_folds", "ss_folds"))
  expect_equal(length(output$folds$cf_folds), 100)
  expect_equal(length(output$folds$ss_folds), 100)
  expect_equal(sort(unique(output$folds$cf_folds)), c(1,2,3,4))
  expect_equal(sort(unique(output$folds$ss_folds)), c(0,1))
  expect_equal(output$approx_times, sort(unique(c(quantile(y[delta == 1 & y <= max(landmark_times)],
                                                           probs = seq(0, 1, by = 0.01)),
                                                         landmark_times))))
  expect_equal(names(output$conditional_surv_preds), c("S_hat", "S_hat_train", "G_hat", "G_hat_train"))
  expect_equal(dim(output$conditional_surv_preds$S_hat[[1]]), c(25, length(output$approx_times)))
  expect_equal(dim(output$conditional_surv_preds$S_hat_train[[1]]), c(75, length(output$approx_times)))
  expect_equal(dim(output$conditional_surv_preds$G_hat[[1]]), c(25, length(output$approx_times)))
  expect_equal(dim(output$conditional_surv_preds$G_hat_train[[1]]), c(75, length(output$approx_times)))
  expect_equal(names(output$large_oracle_preds), c("f_hat", "f_hat_train"))
  expect_equal(dim(output$large_oracle_preds$f_hat[[1]]), c(25, length(landmark_times)))
  expect_equal(dim(output$large_oracle_preds$f_hat_train[[1]]), c(75, length(landmark_times)))
  expect_equal(names(output$small_oracle_preds), c("f_hat", "f_hat_train"))
  expect_equal(dim(output$small_oracle_preds$f_hat[[1]]), c(25, length(landmark_times)))
  expect_equal(dim(output$small_oracle_preds$f_hat_train[[1]]), c(75, length(landmark_times)))
})

output <- vim(type = "accuracy",
              time = y,
              event = delta,
              X = X,
              landmark_times = landmark_times,
              large_feature_vector = 1:2,
              small_feature_vector = 1,
              conditional_surv_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                        bin_size = 0.1,
                                                        V = 2),
              large_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                    V = 2),
              small_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                    V = 2),
              cf_fold_num = 2,
              sample_split = TRUE,
              scale_est = TRUE)
test_that("vim(). accuracy, xfit, sample split", {
  expect_equal(dim(output$result)[1], 3)
  expect_equal(dim(output$result)[2], 12)
  expect_equal(names(output$result), c("landmark_time", "est", "var_est", "cil", "ciu",
                                       "cil_1sided", "p", "large_predictiveness", "small_predictiveness",
                                       "vim", "large_feature_vector", "small_feature_vector"))
  expect_equal(sum(is.na(output$result)), 0)
  expect_equal(names(output$folds), c("cf_folds", "ss_folds"))
  expect_equal(length(output$folds$cf_folds), 100)
  expect_equal(length(output$folds$ss_folds), 100)
  expect_equal(sort(unique(output$folds$cf_folds)), c(1,2,3,4))
  expect_equal(sort(unique(output$folds$ss_folds)), c(0,1))
  expect_equal(output$approx_times, sort(unique(c(quantile(y[delta == 1 & y <= max(landmark_times)],
                                                           probs = seq(0, 1, by = 0.01)),
                                                         landmark_times))))
  expect_equal(names(output$conditional_surv_preds), c("S_hat", "S_hat_train", "G_hat", "G_hat_train"))
  expect_equal(dim(output$conditional_surv_preds$S_hat[[1]]), c(25, length(output$approx_times)))
  expect_equal(dim(output$conditional_surv_preds$S_hat_train[[1]]), c(75, length(output$approx_times)))
  expect_equal(dim(output$conditional_surv_preds$G_hat[[1]]), c(25, length(output$approx_times)))
  expect_equal(dim(output$conditional_surv_preds$G_hat_train[[1]]), c(75, length(output$approx_times)))
  expect_equal(names(output$large_oracle_preds), c("f_hat", "f_hat_train"))
  expect_equal(dim(output$large_oracle_preds$f_hat[[1]]), c(25, length(landmark_times)))
  expect_equal(dim(output$large_oracle_preds$f_hat_train[[1]]), c(75, length(landmark_times)))
  expect_equal(names(output$small_oracle_preds), c("f_hat", "f_hat_train"))
  expect_equal(dim(output$small_oracle_preds$f_hat[[1]]), c(25, length(landmark_times)))
  expect_equal(dim(output$small_oracle_preds$f_hat_train[[1]]), c(75, length(landmark_times)))
})

output <- vim(type = "Brier",
              time = y,
              event = delta,
              X = X,
              landmark_times = landmark_times,
              large_feature_vector = 1:2,
              small_feature_vector = 1,
              conditional_surv_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                        bin_size = 0.1,
                                                        V = 2),
              large_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                    V = 2),
              small_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                    V = 2),
              cf_fold_num = 2,
              sample_split = TRUE,
              scale_est = TRUE)
test_that("vim(). Brier, xfit, sample split", {
  expect_equal(dim(output$result)[1], 3)
  expect_equal(dim(output$result)[2], 12)
  expect_equal(names(output$result), c("landmark_time", "est", "var_est", "cil", "ciu",
                                       "cil_1sided", "p", "large_predictiveness", "small_predictiveness",
                                       "vim", "large_feature_vector", "small_feature_vector"))
  expect_equal(sum(is.na(output$result)), 0)
  expect_equal(names(output$folds), c("cf_folds", "ss_folds"))
  expect_equal(length(output$folds$cf_folds), 100)
  expect_equal(length(output$folds$ss_folds), 100)
  expect_equal(sort(unique(output$folds$cf_folds)), c(1,2,3,4))
  expect_equal(sort(unique(output$folds$ss_folds)), c(0,1))
  expect_equal(output$approx_times, sort(unique(c(quantile(y[delta == 1 & y <= max(landmark_times)],
                                                           probs = seq(0, 1, by = 0.01)),
                                                         landmark_times))))
  expect_equal(names(output$conditional_surv_preds), c("S_hat", "S_hat_train", "G_hat", "G_hat_train"))
  expect_equal(dim(output$conditional_surv_preds$S_hat[[1]]), c(25, length(output$approx_times)))
  expect_equal(dim(output$conditional_surv_preds$S_hat_train[[1]]), c(75, length(output$approx_times)))
  expect_equal(dim(output$conditional_surv_preds$G_hat[[1]]), c(25, length(output$approx_times)))
  expect_equal(dim(output$conditional_surv_preds$G_hat_train[[1]]), c(75, length(output$approx_times)))
  expect_equal(names(output$large_oracle_preds), c("f_hat", "f_hat_train"))
  expect_equal(dim(output$large_oracle_preds$f_hat[[1]]), c(25, length(landmark_times)))
  expect_equal(dim(output$large_oracle_preds$f_hat_train[[1]]), c(75, length(landmark_times)))
  expect_equal(names(output$small_oracle_preds), c("f_hat", "f_hat_train"))
  expect_equal(dim(output$small_oracle_preds$f_hat[[1]]), c(25, length(landmark_times)))
  expect_equal(dim(output$small_oracle_preds$f_hat_train[[1]]), c(75, length(landmark_times)))
})

output <- vim(type = "R-squared",
              time = y,
              event = delta,
              X = X,
              landmark_times = landmark_times,
              large_feature_vector = 1:2,
              small_feature_vector = 1,
              conditional_surv_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                        bin_size = 0.1,
                                                        V = 2),
              large_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                    V = 2),
              small_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                    V = 2),
              cf_fold_num = 2,
              sample_split = TRUE,
              scale_est = TRUE)
test_that("vim(). R-squared, xfit, sample split", {
  expect_equal(dim(output$result)[1], 3)
  expect_equal(dim(output$result)[2], 12)
  expect_equal(names(output$result), c("landmark_time", "est", "var_est", "cil", "ciu",
                                       "cil_1sided", "p", "large_predictiveness", "small_predictiveness",
                                       "vim", "large_feature_vector", "small_feature_vector"))
  expect_equal(sum(is.na(output$result)), 0)
  expect_equal(names(output$folds), c("cf_folds", "ss_folds"))
  expect_equal(length(output$folds$cf_folds), 100)
  expect_equal(length(output$folds$ss_folds), 100)
  expect_equal(sort(unique(output$folds$cf_folds)), c(1,2,3,4))
  expect_equal(sort(unique(output$folds$ss_folds)), c(0,1))
  expect_equal(output$approx_times, sort(unique(c(quantile(y[delta == 1 & y <= max(landmark_times)],
                                                           probs = seq(0, 1, by = 0.01)),
                                                         landmark_times))))
  expect_equal(names(output$conditional_surv_preds), c("S_hat", "S_hat_train", "G_hat", "G_hat_train"))
  expect_equal(dim(output$conditional_surv_preds$S_hat[[1]]), c(25, length(output$approx_times)))
  expect_equal(dim(output$conditional_surv_preds$S_hat_train[[1]]), c(75, length(output$approx_times)))
  expect_equal(dim(output$conditional_surv_preds$G_hat[[1]]), c(25, length(output$approx_times)))
  expect_equal(dim(output$conditional_surv_preds$G_hat_train[[1]]), c(75, length(output$approx_times)))
  expect_equal(names(output$large_oracle_preds), c("f_hat", "f_hat_train"))
  expect_equal(dim(output$large_oracle_preds$f_hat[[1]]), c(25, length(landmark_times)))
  expect_equal(dim(output$large_oracle_preds$f_hat_train[[1]]), c(75, length(landmark_times)))
  expect_equal(names(output$small_oracle_preds), c("f_hat", "f_hat_train"))
  expect_equal(dim(output$small_oracle_preds$f_hat[[1]]), c(25, length(landmark_times)))
  expect_equal(dim(output$small_oracle_preds$f_hat_train[[1]]), c(75, length(landmark_times)))
})

output <- vim(type = "survival_time_MSE",
              time = y,
              event = delta,
              X = X,
              restriction_time = landmark_times[3],
              large_feature_vector = 1:2,
              small_feature_vector = 1,
              conditional_surv_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                        bin_size = 0.1,
                                                        V = 2),
              large_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                    V = 2),
              small_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                    V = 2),
              cf_fold_num = 2,
              sample_split = TRUE,
              scale_est = TRUE)
test_that("vim(). survival time MSE, xfit, sample split", {
  expect_equal(dim(output$result)[1], 1)
  expect_equal(dim(output$result)[2], 12)
  expect_equal(names(output$result), c("restriction_time", "est", "var_est", "cil", "ciu",
                                       "cil_1sided", "p", "large_predictiveness", "small_predictiveness",
                                       "vim", "large_feature_vector", "small_feature_vector"))
  expect_equal(sum(is.na(output$result)), 0)
  expect_equal(names(output$folds), c("cf_folds", "ss_folds"))
  expect_equal(length(output$folds$cf_folds), 100)
  expect_equal(length(output$folds$ss_folds), 100)
  expect_equal(sort(unique(output$folds$cf_folds)), c(1,2,3,4))
  expect_equal(sort(unique(output$folds$ss_folds)), c(0,1))
  expect_equal(output$approx_times, sort(unique(c(quantile(y[delta == 1 & y <= landmark_times[3]],
                                                           probs = seq(0, 1, by = 0.01)),
                                                         landmark_times[3]))))
  expect_equal(names(output$conditional_surv_preds), c("S_hat", "S_hat_train", "G_hat", "G_hat_train"))
  expect_equal(dim(output$conditional_surv_preds$S_hat[[1]]), c(25, length(output$approx_times)))
  expect_equal(dim(output$conditional_surv_preds$S_hat_train[[1]]), c(75, length(output$approx_times)))
  expect_equal(dim(output$conditional_surv_preds$G_hat[[1]]), c(25, length(output$approx_times)))
  expect_equal(dim(output$conditional_surv_preds$G_hat_train[[1]]), c(75, length(output$approx_times)))
  expect_equal(names(output$large_oracle_preds), c("f_hat", "f_hat_train"))
  expect_equal(length(output$large_oracle_preds$f_hat[[1]]), 25)
  expect_equal(length(output$large_oracle_preds$f_hat_train[[1]]), 75)
  expect_equal(names(output$small_oracle_preds), c("f_hat", "f_hat_train"))
  expect_equal(length(output$small_oracle_preds$f_hat[[1]]), 25)
  expect_equal(length(output$small_oracle_preds$f_hat_train[[1]]), 75)
})

output <- vim(type = "C-index",
              time = y,
              event = delta,
              X = X,
              restriction_time = landmark_times[3],
              large_feature_vector = 1:2,
              small_feature_vector = 1,
              conditional_surv_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                        bin_size = 0.1,
                                                        V = 2),
              large_oracle_generator_control = list(V = 2,
                                                    params = list(mstop = c(20),
                                                                  nu = c(0.1),
                                                                  sigma = c(0.01),
                                                                  learner = c("glm"))),
              small_oracle_generator_control = list(V = 2,
                                                    params = list(mstop = c(20),
                                                                  nu = c(0.1),
                                                                  sigma = c(0.01),
                                                                  learner = c("glm"))),
              cf_fold_num = 2,
              sample_split = TRUE,
              scale_est = TRUE)
test_that("vim(). C-index, xfit, sample split", {
  expect_equal(dim(output$result)[1], 1)
  expect_equal(dim(output$result)[2], 12)
  expect_equal(names(output$result), c("restriction_time", "est", "var_est", "cil", "ciu",
                                       "cil_1sided", "p", "large_predictiveness", "small_predictiveness",
                                       "vim", "large_feature_vector", "small_feature_vector"))
  expect_equal(sum(is.na(output$result)), 0)
  expect_equal(names(output$folds), c("cf_folds", "ss_folds"))
  expect_equal(length(output$folds$cf_folds), 100)
  expect_equal(length(output$folds$ss_folds), 100)
  expect_equal(sort(unique(output$folds$cf_folds)), c(1,2,3,4))
  expect_equal(sort(unique(output$folds$ss_folds)), c(0,1))
  expect_equal(output$approx_times, sort(unique(c(quantile(y[delta == 1 & y <= landmark_times[3]],
                                                           probs = seq(0, 1, by = 0.01)),
                                                         landmark_times[3]))))
  expect_equal(names(output$conditional_surv_preds), c("S_hat", "S_hat_train", "G_hat", "G_hat_train"))
  expect_equal(dim(output$conditional_surv_preds$S_hat[[1]]), c(25, length(output$approx_times)))
  expect_equal(dim(output$conditional_surv_preds$S_hat_train[[1]]), c(75, length(output$approx_times)))
  expect_equal(dim(output$conditional_surv_preds$G_hat[[1]]), c(25, length(output$approx_times)))
  expect_equal(dim(output$conditional_surv_preds$G_hat_train[[1]]), c(75, length(output$approx_times)))
  expect_equal(names(output$large_oracle_preds), c("f_hat", "f_hat_train"))
  expect_equal(length(output$large_oracle_preds$f_hat[[1]]), 25)
  expect_equal(length(output$large_oracle_preds$f_hat_train[[1]]), 75)
  expect_equal(names(output$small_oracle_preds), c("f_hat", "f_hat_train"))
  expect_equal(length(output$small_oracle_preds$f_hat[[1]]), 25)
  expect_equal(length(output$small_oracle_preds$f_hat_train[[1]]), 75)
})

output <- vim(type = "C-index",
              time = y,
              event = delta,
              X = X,
              restriction_time = landmark_times[3],
              large_feature_vector = 1:2,
              small_feature_vector = 1,
              conditional_surv_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                        bin_size = 0.1,
                                                        V = 2),
              large_oracle_generator_control = list(V = 2,
                                                    tuning = "CV",
                                                    params = list(mstop = c(20),
                                                                  nu = c(0.1),
                                                                  sigma = c(0.01),
                                                                  learner = c("glm"))),
              small_oracle_generator_control = list(V = 2,
                                                    tuning = "CV",
                                                    params = list(mstop = c(20),
                                                                  nu = c(0.1),
                                                                  sigma = c(0.01),
                                                                  learner = c("glm"))),
              cf_fold_num = 2,
              sample_split = TRUE,
              scale_est = TRUE)
test_that("vim(). C-index with CV, xfit, sample split", {
  expect_equal(dim(output$result)[1], 1)
  expect_equal(dim(output$result)[2], 12)
  expect_equal(names(output$result), c("restriction_time", "est", "var_est", "cil", "ciu",
                                       "cil_1sided", "p", "large_predictiveness", "small_predictiveness",
                                       "vim", "large_feature_vector", "small_feature_vector"))
  expect_equal(sum(is.na(output$result)), 0)
  expect_equal(names(output$folds), c("cf_folds", "ss_folds"))
  expect_equal(length(output$folds$cf_folds), 100)
  expect_equal(length(output$folds$ss_folds), 100)
  expect_equal(sort(unique(output$folds$cf_folds)), c(1,2,3,4))
  expect_equal(sort(unique(output$folds$ss_folds)), c(0,1))
  expect_equal(output$approx_times, sort(unique(c(quantile(y[delta == 1 & y <= landmark_times[3]],
                                                           probs = seq(0, 1, by = 0.01)),
                                                         landmark_times[3]))))
  expect_equal(names(output$conditional_surv_preds), c("S_hat", "S_hat_train", "G_hat", "G_hat_train"))
  expect_equal(dim(output$conditional_surv_preds$S_hat[[1]]), c(25, length(output$approx_times)))
  expect_equal(dim(output$conditional_surv_preds$S_hat_train[[1]]), c(75, length(output$approx_times)))
  expect_equal(dim(output$conditional_surv_preds$G_hat[[1]]), c(25, length(output$approx_times)))
  expect_equal(dim(output$conditional_surv_preds$G_hat_train[[1]]), c(75, length(output$approx_times)))
  expect_equal(names(output$large_oracle_preds), c("f_hat", "f_hat_train"))
  expect_equal(length(output$large_oracle_preds$f_hat[[1]]), 25)
  expect_equal(length(output$large_oracle_preds$f_hat_train[[1]]), 75)
  expect_equal(names(output$small_oracle_preds), c("f_hat", "f_hat_train"))
  expect_equal(length(output$small_oracle_preds$f_hat[[1]]), 25)
  expect_equal(length(output$small_oracle_preds$f_hat_train[[1]]), 75)
})

