set.seed(1)
y <- rexp(100, 1)
delta <- rbinom(100, size = 1, prob = 0.5)
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
                       folds = folds,
                       ss_folds = ss_folds,
                       sample_split = FALSE)
test_that("vim_accuracy(). no xfit, no sample split", {
  expect_equal(dim(output)[1], 3)
  expect_equal(dim(output)[2], 11)
  expect_equal(names(output), c("tau", "full_one_step", "reduced_one_step", "one_step",
                                "full_plug_in", "reduced_plug_in", "var_est", "cil", "ciu",
                                "cil_1sided", "p"))
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
                  folds = folds,
                  ss_folds = ss_folds,
                  sample_split = FALSE)
test_that("vim_AUC(). no xfit, no sample split", {
  expect_equal(dim(output)[1], 3)
  expect_equal(dim(output)[2], 11)
  expect_equal(names(output), c("tau", "full_one_step", "reduced_one_step", "one_step",
                                "full_plug_in", "reduced_plug_in", "var_est", "cil", "ciu",
                                "cil_1sided", "p"))
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
                    folds = folds,
                    ss_folds = ss_folds,
                    sample_split = FALSE)
test_that("vim_brier(). no xfit, no sample split", {
  expect_equal(dim(output)[1], 3)
  expect_equal(dim(output)[2], 11)
  expect_equal(names(output), c("tau", "full_one_step", "reduced_one_step", "one_step",
                                "full_plug_in", "reduced_plug_in", "var_est", "cil", "ciu",
                                "cil_1sided", "p"))
  expect_equal(sum(is.na(output)), 0)
})

# RMST MSE
output <- vim_rmst_mse(time = y,
                       event = delta,
                       approx_times = approx_times,
                       tau = landmark_times[3],
                       f_hat = lapply(f_hat, function(x) x[,1]),
                       fs_hat = lapply(fs_hat, function(x) x[,1]),
                       S_hat = S_hat,
                       G_hat = G_hat,
                       folds = folds,
                       ss_folds = ss_folds,
                       sample_split = FALSE)
test_that("vim_rmst_mse(). no xfit, no sample split", {
  expect_equal(dim(output)[1], 1)
  expect_equal(dim(output)[2], 11)
  expect_equal(names(output), c("tau", "full_one_step", "reduced_one_step", "one_step",
                                "full_plug_in", "reduced_plug_in", "var_est", "cil", "ciu",
                                "cil_1sided", "p"))
  expect_equal(sum(is.na(output)), 0)
})

# C-index
output <- vim_cindex(time = y,
                     event = delta,
                     approx_times = approx_times,
                     tau = landmark_times[3],
                     f_hat = lapply(f_hat, function(x) x[,1]),
                     fs_hat = lapply(fs_hat, function(x) x[,1]),
                     S_hat = S_hat,
                     G_hat = G_hat,
                     folds = folds,
                     ss_folds = ss_folds,
                     sample_split = FALSE)
test_that("vim_cindex(). no xfit, no sample split", {
  expect_equal(dim(output)[1], 1)
  expect_equal(dim(output)[2], 11)
  expect_equal(names(output), c("tau", "full_one_step", "reduced_one_step", "one_step",
                                "full_plug_in", "reduced_plug_in", "var_est", "cil", "ciu",
                                "cil_1sided", "p"))
  expect_equal(sum(is.na(output)), 0)
})


################################
### no crossfit, sample split
################################
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
                       folds = folds,
                       ss_folds = ss_folds,
                       sample_split = TRUE)
test_that("vim_accuracy(). no xfit, sample split", {
  expect_equal(dim(output)[1], 3)
  expect_equal(dim(output)[2], 11)
  expect_equal(names(output), c("tau", "full_one_step", "reduced_one_step", "one_step",
                                "full_plug_in", "reduced_plug_in", "var_est", "cil", "ciu",
                                "cil_1sided", "p"))
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
                  folds = folds,
                  ss_folds = ss_folds,
                  sample_split = TRUE)
test_that("vim_AUC(). no xfit, sample split", {
  expect_equal(dim(output)[1], 3)
  expect_equal(dim(output)[2], 11)
  expect_equal(names(output), c("tau", "full_one_step", "reduced_one_step", "one_step",
                                "full_plug_in", "reduced_plug_in", "var_est", "cil", "ciu",
                                "cil_1sided", "p"))
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
                    folds = folds,
                    ss_folds = ss_folds,
                    sample_split = TRUE)
test_that("vim_brier(). no xfit, sample split", {
  expect_equal(dim(output)[1], 3)
  expect_equal(dim(output)[2], 11)
  expect_equal(names(output), c("tau", "full_one_step", "reduced_one_step", "one_step",
                                "full_plug_in", "reduced_plug_in", "var_est", "cil", "ciu",
                                "cil_1sided", "p"))
  expect_equal(sum(is.na(output)), 0)
})

# RMST MSE
output <- vim_rmst_mse(time = y,
                       event = delta,
                       approx_times = approx_times,
                       tau = landmark_times[3],
                       f_hat = lapply(f_hat, function(x) x[,1]),
                       fs_hat = lapply(fs_hat, function(x) x[,1]),
                       S_hat = S_hat,
                       G_hat = G_hat,
                       folds = folds,
                       ss_folds = ss_folds,
                       sample_split = TRUE)
test_that("vim_rmst_mse(). no xfit, sample split", {
  expect_equal(dim(output)[1], 1)
  expect_equal(dim(output)[2], 11)
  expect_equal(names(output), c("tau", "full_one_step", "reduced_one_step", "one_step",
                                "full_plug_in", "reduced_plug_in", "var_est", "cil", "ciu",
                                "cil_1sided", "p"))
  expect_equal(sum(is.na(output)), 0)
})

# C-index
output <- vim_cindex(time = y,
                     event = delta,
                     approx_times = approx_times,
                     tau = landmark_times[3],
                     f_hat = lapply(f_hat, function(x) x[,1]),
                     fs_hat = lapply(fs_hat, function(x) x[,1]),
                     S_hat = S_hat,
                     G_hat = G_hat,
                     folds = folds,
                     ss_folds = ss_folds,
                     sample_split = TRUE)
test_that("vim_cindex(). no xfit, sample split", {
  expect_equal(dim(output)[1], 1)
  expect_equal(dim(output)[2], 11)
  expect_equal(names(output), c("tau", "full_one_step", "reduced_one_step", "one_step",
                                "full_plug_in", "reduced_plug_in", "var_est", "cil", "ciu",
                                "cil_1sided", "p"))
  expect_equal(sum(is.na(output)), 0)
})
