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

SL.library <- c("SL.gam", "SL.glm")

fit <- stackL(time = time,
              event = event,
              entry = entry,
              X = X,
              newX = X[c(1,2,3),],
              newtimes = seq(0, 15, 3),
              direction = "prospective",
              bin_size = 0.05,
              time_basis = "continuous",
              SL_control = list(SL.library = SL.library,
                                V = 5))

true_S_T_preds <- rbind(c(0.974, 0.974, 0.630, 0.385, 0.319, 0.283),
                        c(0.985, 0.985, 0.769, 0.581, 0.522, 0.487),
                        c(0.962, 0.962, 0.508, 0.247, 0.188, 0.158))

estimated_S_T_preds <- round(fit$S_T_preds, digits = 3)

diffs <- (abs(true_S_T_preds - estimated_S_T_preds) > 0.01)

test_that("stackL() returns expected output (truncation)", {
  expect_equal(sum(diffs), 0)
})

preds <- predict(fit,
                 newX = X[c(1,2,3),],
                 newtimes = seq(0, 15, 3))

estimated_S_T_preds <- round(preds$S_T_preds, digits = 3)

diffs <- (abs(true_S_T_preds - estimated_S_T_preds) > 0.01)


test_that("stackL() predict() method is not broken", {
  expect_equal(sum(diffs), 0)
})
