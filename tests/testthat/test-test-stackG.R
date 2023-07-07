# ensure basic functioning doesn't change

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

SL.library <- c("SL.mean", "SL.gam")

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

true_S_T_preds <- rbind(c(1, 0.796, 0.575, 0.532, 0.529, 0.528),
                        c(1, 0.849, 0.840, 0.745, 0.646, 0.546),
                        c(1, 0.795, 0.744, 0.505, 0.420, 0.377))

estimated_S_T_preds <- round(fit$S_T_preds, digits = 3)

diffs <- (abs(true_S_T_preds - estimated_S_T_preds) > 0.01)

test_that("stackG() returns expected output (truncation)", {
  expect_equal(sum(diffs), 0)
})

preds <- predict(fit,
                 newX = X[c(1,2,3),],
                 newtimes = seq(0, 15, 3))

estimated_S_T_preds <- round(fit$S_T_preds, digits = 3)

diffs <- (abs(true_S_T_preds - estimated_S_T_preds) > 0.01)

test_that("stackG() predict() method is not broken (truncation)", {
  expect_equal(sum(diffs), 0)
})


# no truncation
set.seed(1)
n <- 100
X <- data.frame(X1 = rnorm(n), X2 = rbinom(n, size = 1, prob = 0.5))
T <- rexp(n, rate = exp(-2 + X[,1] - X[,2] + .5 *  X[,1] * X[,2]))
C <- rexp(n, exp(-2 -.5 * X[,1] - .25 * X[,2] + .5 * X[,1] * X[,2]))
C[C > 15] <- 15

time <- pmin(T, C)
event <- as.numeric(T <= C)

SL.library <- c("SL.mean", "SL.gam")

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

true_S_T_preds <- rbind(c(1, 0.889, 0.720, 0.560, 0.502, 0.451),
                        c(1, 0.668, 0.400, 0.243, 0.199, 0.165),
                        c(1, 0.985, 0.953, 0.911, 0.892, 0.872))

estimated_S_T_preds <- round(fit$S_T_preds, digits = 3)

diffs <- (abs(true_S_T_preds - estimated_S_T_preds) > 0.01)

test_that("stackG() returns expected output (no truncation)", {
  expect_equal(sum(diffs), 0)
})


preds <- predict(fit,
                 newX = X[c(1,2,3),],
                 newtimes = seq(0, 15, 3))

estimated_S_T_preds <- round(preds$S_T_preds, digits = 3)

diffs <- (abs(true_S_T_preds - estimated_S_T_preds) > 0.01)

test_that("stackG() predict() method is not broken", {
  expect_equal(sum(diffs), 0)
})
