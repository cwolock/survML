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
              bin_size = 0.02,
              time_basis = "continuous",
              time_grid_approx = sort(unique(time)),
              SL_control = list(SL.library = SL.library, V = 5),
              surv_form = "PI")

true_S_T_preds <- rbind(c(1, 0.791, 0.536, 0.491, 0.489, 0.487),
                        c(1, 0.863, 0.855, 0.761, 0.666, 0.530),
                        c(1, 0.799, 0.747, 0.492, 0.409, 0.355))

test_that("stackG is not broken", {
  expect_equal(round(fit$S_T_preds, digits = 3), true_S_T_preds)
})


preds <- predict(fit,
                 newX = X[c(1,2,3),],
                 newtimes = seq(0, 15, 3))

test_that("stackG predict() method is not broken", {
  expect_equal(round(preds$S_T_preds, digits = 3), true_S_T_preds)
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
              bin_size = 0.02,
              time_basis = "continuous",
              time_grid_approx = sort(unique(time)),
              SL_control = list(SL.library = SL.library, V = 5),
              surv_form = "PI")

true_S_T_preds <- rbind(c(1, 0.791, 0.536, 0.491, 0.489, 0.487),
                        c(1, 0.863, 0.855, 0.761, 0.666, 0.530),
                        c(1, 0.799, 0.747, 0.492, 0.409, 0.355))

test_that("stackG is not broken", {
  expect_equal(round(fit$S_T_preds, digits = 3), true_S_T_preds)
})


preds <- predict(fit,
                 newX = X[c(1,2,3),],
                 newtimes = seq(0, 15, 3))

test_that("stackG predict() method is not broken", {
  expect_equal(round(preds$S_T_preds, digits = 3), true_S_T_preds)
})

