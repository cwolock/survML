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
               bin_size = 0.02,
               time_basis = "continuous",
              SL_control = list(SL.library = SL.library,
                                V = 5))

true_S_T_preds <- rbind(c(0.999, 0.999, 0.812, 0.557, 0.355, 0.309),
                        c(0.999, 0.999, 0.841, 0.613, 0.420, 0.373),
                        c(0.997, 0.997, 0.637, 0.282, 0.111, 0.083))

test_that("survMLc is not broken", {
  expect_equal(round(fit$S_T_preds, digits = 3), true_S_T_preds)
})

preds <- predict(fit,
                 newX = X[c(1,2,3),],
                 newtimes = seq(0, 15, 3))


test_that("survMLs predict() method is not broken", {
  expect_equal(round(preds$S_T_preds, digits = 3), true_S_T_preds)
})
