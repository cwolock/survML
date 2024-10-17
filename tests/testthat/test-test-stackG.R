#####################################
### basic functioning with truncation
#####################################
set.seed(1)
n <- 100
X <- data.frame(X1 = rnorm(n), X2 = rbinom(n, size = 1, prob = 0.5))
T <- rexp(n, rate = exp(-2 + X[,1] - X[,2] + .5 *  X[,1] * X[,2]))
C <- rexp(n, exp(-2 -.5 * X[,1] - .25 * X[,2] + .5 * X[,1] * X[,2]))
C[C > 15] <- 15
entry <- runif(n, 0, 15)

time <- pmin(T, C)
event <- as.numeric(T <= C)

sampled <- which(time >= entry)
X <- X[sampled,]
time <- time[sampled]
event <- event[sampled]
entry <- entry[sampled]

SL.library <- c("SL.mean")

# exponential form
fit <- stackG(time = time,
              event = event,
              entry = entry,
              X = X,
              newX = X[c(1,2,3),],
              newtimes = seq(0, 15, 3),
              direction = "prospective",
              bin_size = 0.05,
              time_basis = "continuous",
              time_grid_approx = sort(unique(time)),
              SL_control = list(SL.library = SL.library,
                                V = 5,
                                method = "method.NNLS"),
              surv_form = "exp")

true_S_T_preds <- rbind(c(1, 0.767, 0.767, 0.767, 0.767, 0.767),
                        c(1, 0.767, 0.767, 0.767, 0.767, 0.767),
                        c(1, 0.767, 0.767, 0.767, 0.767, 0.767))

estimated_S_T_preds <- round(fit$S_T_preds, digits = 3)

diffs <- (abs(true_S_T_preds - estimated_S_T_preds) > 0.01)

test_that("stackG() returns expected output (truncation, exponential)", {
  expect_equal(sum(diffs), 0)
})

preds <- predict(fit,
                 newX = X[c(1,2,3),],
                 newtimes = seq(0, 15, 3))

estimated_S_T_preds <- round(preds$S_T_preds, digits = 3)

diffs <- (abs(true_S_T_preds - estimated_S_T_preds) > 0.01)

test_that("stackG() predict() method is not broken (truncation, exponential)", {
  expect_equal(sum(diffs), 0)
})

# PI form
fit <- stackG(time = time,
              event = event,
              entry = entry,
              X = X,
              newX = X[c(1,2,3),],
              newtimes = seq(0, 15, 3),
              direction = "prospective",
              bin_size = 0.05,
              time_basis = "continuous",
              time_grid_approx = sort(unique(time)),
              SL_control = list(SL.library = SL.library,
                                V = 5,
                                method = "method.NNLS"),
              surv_form = "PI")


true_S_T_preds <- rbind(c(1, 0.735, 0.735, 0.735, 0.735, 0.735),
                        c(1, 0.735, 0.735, 0.735, 0.735, 0.735),
                        c(1, 0.735, 0.735, 0.735, 0.735, 0.735))

estimated_S_T_preds <- round(fit$S_T_preds, digits = 3)

diffs <- (abs(true_S_T_preds - estimated_S_T_preds) > 0.01)

test_that("stackG() returns expected output (truncation, PI)", {
  expect_equal(sum(diffs), 0)
})

preds <- predict(fit,
                 newX = X[c(1,2,3),],
                 newtimes = seq(0, 15, 3),
                 surv_form = "PI")

estimated_S_T_preds <- round(preds$S_T_preds, digits = 3)

diffs <- (abs(true_S_T_preds - estimated_S_T_preds) > 0.01)

test_that("stackG() predict() method is not broken (truncation, PI)", {
  expect_equal(sum(diffs), 0)
})

# dummy time
set.seed(1)
suppressWarnings({
  fit <- stackG(time = time,
                event = event,
                entry = entry,
                X = X,
                newX = X[c(1,2,3),],
                newtimes = seq(0, 15, 3),
                direction = "prospective",
                bin_size = 0.05,
                time_basis = "dummy",
                time_grid_approx = sort(unique(time)),
                SL_control = list(SL.library = SL.library,
                                  V = 5,
                                  method = "method.NNLS"),
                surv_form = "PI")
})

true_S_T_preds <- rbind(c(1,0.735, 0.735, 0.735, 0.735, 0.735),
                        c(1,  0.735, 0.735, 0.735, 0.735, 0.735),
                        c(1, 0.735, 0.735, 0.735, 0.735, 0.735))

estimated_S_T_preds <- round(fit$S_T_preds, digits = 3)

diffs <- (abs(true_S_T_preds - estimated_S_T_preds) > 0.01)

test_that("stackG() returns expected output (truncation, PI, dummy)", {
  expect_equal(sum(diffs), 0)
})


########################################
### basic functioning without truncation
########################################
set.seed(1)
n <- 100
X <- data.frame(X1 = rnorm(n), X2 = rbinom(n, size = 1, prob = 0.5))
T <- rexp(n, rate = exp(-2 + X[,1] - X[,2] + .5 *  X[,1] * X[,2]))
C <- rexp(n, exp(-2 -.5 * X[,1] - .25 * X[,2] + .5 * X[,1] * X[,2]))
C[C > 15] <- 15

time <- pmin(T, C)
event <- as.numeric(T <= C)

SL.library <- c("SL.mean")

# exponential form
fit <- stackG(time = time,
              event = event,
              X = X,
              newX = X[c(1,2,3),],
              newtimes = seq(0, 15, 3),
              direction = "prospective",
              bin_size = 0.05,
              time_basis = "continuous",
              time_grid_approx = sort(unique(time)),
              SL_control = list(SL.library = SL.library,
                                V = 5,
                                method = "method.NNLS"),
              surv_form = "exp")

true_S_T_preds <- rbind(c(1, 0.799, 0.799, 0.799, 0.799, 0.799),
                        c(1, 0.799, 0.799, 0.799, 0.799, 0.799),
                        c(1, 0.799, 0.799, 0.799, 0.799, 0.799))

estimated_S_T_preds <- round(fit$S_T_preds, digits = 3)

diffs <- (abs(true_S_T_preds - estimated_S_T_preds) > 0.01)

test_that("stackG() returns expected output (no truncation, exponential)", {
  expect_equal(sum(diffs), 0)
})

preds <- predict(fit,
                 newX = X[c(1,2,3),],
                 newtimes = seq(0, 15, 3))

estimated_S_T_preds <- round(preds$S_T_preds, digits = 3)

diffs <- (abs(true_S_T_preds - estimated_S_T_preds) > 0.01)

test_that("stackG() predict() method is not broken (no truncation, exponential)", {
  expect_equal(sum(diffs), 0)
})

# PI form
fit <- stackG(time = time,
              event = event,
              X = X,
              newX = X[c(1,2,3),],
              newtimes = seq(0, 15, 3),
              direction = "prospective",
              bin_size = 0.05,
              time_basis = "continuous",
              time_grid_approx = sort(unique(time)),
              SL_control = list(SL.library = SL.library,
                                V = 5,
                                method = "method.NNLS"),
              surv_form = "PI")

true_S_T_preds <- rbind(c(1, 0.776, 0.776, 0.776, 0.776, 0.776),
                        c(1, 0.776, 0.776, 0.776, 0.776, 0.776),
                        c(1, 0.776, 0.776, 0.776, 0.776, 0.776))

estimated_S_T_preds <- round(fit$S_T_preds, digits = 3)

diffs <- (abs(true_S_T_preds - estimated_S_T_preds) > 0.01)

test_that("stackG() returns expected output (no truncation, PI)", {
  expect_equal(sum(diffs), 0)
})

preds <- predict(fit,
                 newX = X[c(1,2,3),],
                 newtimes = seq(0, 15, 3),
                 surv_form = "PI")

estimated_S_T_preds <- round(preds$S_T_preds, digits = 3)

diffs <- (abs(true_S_T_preds - estimated_S_T_preds) > 0.01)

test_that("stackG() predict() method is not broken (no truncation, PI)", {
  expect_equal(sum(diffs), 0)
})

# dummy time
set.seed(1)
suppressWarnings({
  fit <- stackG(time = time,
                event = event,
                X = X,
                newX = X[c(1,2,3),],
                newtimes = seq(0, 15, 3),
                direction = "prospective",
                bin_size = 0.05,
                time_basis = "dummy",
                time_grid_approx = sort(unique(time)),
                SL_control = list(SL.library = SL.library,
                                  V = 5,
                                  method = "method.NNLS"),
                surv_form = "PI")
})

true_S_T_preds <- rbind(c(1, 0.776, 0.776, 0.776, 0.776, 0.776),
                        c(1, 0.776, 0.776, 0.776, 0.776, 0.776),
                        c(1, 0.776, 0.776, 0.776, 0.776, 0.776))

estimated_S_T_preds <- round(fit$S_T_preds, digits = 3)

diffs <- (abs(true_S_T_preds - estimated_S_T_preds) > 0.01)

test_that("stackG() returns expected output (no truncation, PI, dummy)", {
  expect_equal(sum(diffs), 0)
})

####################################################
### basic functioning with truncation, retrospective
####################################################
set.seed(1)
n <- 100
X <- data.frame(X1 = rnorm(n), X2 = rbinom(n, size = 1, prob = 0.5))
T <- rexp(n, rate = exp(-2 + X[,1] - X[,2] + .5 *  X[,1] * X[,2]))
entry <- runif(n, 0, 15)

time <- T

sampled <- which(time <= entry)
X <- X[sampled,]
time <- time[sampled]
entry <- entry[sampled]

SL.library <- c("SL.mean")

# exponential form
fit <- stackG(time = time,
              entry = entry,
              X = X,
              newX = X[c(1,2,3),],
              newtimes = seq(0, 15, 3),
              direction = "retrospective",
              bin_size = 0.05,
              time_basis = "continuous",
              time_grid_approx = sort(unique(time)),
              SL_control = list(SL.library = SL.library,
                                V = 5,
                                method = "method.NNLS"),
              surv_form = "exp")

true_S_T_preds <- rbind(c(0.426, 0.426, 0.426, 0.426, 0, 0),
                        c(0.426, 0.426, 0.426, 0.426, 0, 0),
                        c(0.426, 0.426, 0.426, 0.426, 0, 0))

estimated_S_T_preds <- round(fit$S_T_preds, digits = 3)

diffs <- (abs(true_S_T_preds - estimated_S_T_preds) > 0.01)

test_that("stackG() returns expected output (right truncation, exponential)", {
  expect_equal(sum(diffs), 0)
})
