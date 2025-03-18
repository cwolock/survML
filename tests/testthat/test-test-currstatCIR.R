set.seed(123)
n <- 200
x <- cbind(2*rbinom(n, size = 1, prob = 0.5)-1,
           2*rbinom(n, size = 1, prob = 0.5)-1)
t <- rweibull(n,
              shape = 0.75,
              scale = exp(0.4*x[,1] - 0.2*x[,2]))
y <- rweibull(n,
              shape = 0.75,
              scale = exp(0.4*x[,1] - 0.2*x[,2]))

# round y to nearest quantile of y, just so there aren't so many unique values
quants <- quantile(y, probs = seq(0, 1, by = 0.1), type = 1)
for (i in 1:length(y)){
  y[i] <- quants[which.min(abs(y[i] - quants))]
}
delta <- as.numeric(t <= y)

dat <- data.frame(y = y, delta = delta, x1 = x[,1], x2 = x[,2])

dat$delta[dat$y > 1.8] <- NA
dat$y[dat$y > 1.8] <- NA
eval_region <- c(0.05, 1.5)
res <- survML::currstatCIR(time = dat$y,
                           event = dat$delta,
                           X = dat[,3:4],
                           SL_control = list(SL.library = c("SL.mean", "SL.glm"),
                                             V = 2),
                           HAL_control = list(n_bins = c(5),
                                              grid_type = c("equal_mass"),
                                              V = 2),
                           eval_region = eval_region)$primary_results

test_that("currstatCIR()", {
  expect_equal(dim(res), c(100, 4))
  expect_equal(names(res), c("t", "S_hat_est", "S_hat_cil", "S_hat_ciu"))
})
