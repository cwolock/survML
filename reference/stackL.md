# Estimate a conditional survival function via local survival stacking

Estimate a conditional survival function via local survival stacking

## Usage

``` r
stackL(
  time,
  event = rep(1, length(time)),
  entry = NULL,
  X,
  newX,
  newtimes,
  direction = "prospective",
  bin_size = NULL,
  time_basis = "continuous",
  learner = "SuperLearner",
  SL_control = list(SL.library = c("SL.mean"), V = 10, method = "method.NNLS", stratifyCV
    = FALSE),
  tau = NULL
)
```

## Arguments

- time:

  `n x 1` numeric vector of observed follow-up times If there is
  censoring, these are the minimum of the event and censoring times.

- event:

  `n x 1` numeric vector of status indicators of whether an event was
  observed. Defaults to a vector of 1s, i.e. no censoring.

- entry:

  Study entry variable, if applicable. Defaults to `NULL`, indicating
  that there is no truncation.

- X:

  `n x p` data.frame of observed covariate values on which to train the
  estimator.

- newX:

  `m x p` data.frame of new observed covariate values at which to obtain
  `m` predictions for the estimated algorithm. Must have the same names
  and structure as `X`.

- newtimes:

  `k x 1` numeric vector of times at which to obtain `k` predicted
  conditional survivals.

- direction:

  Whether the data come from a prospective or retrospective study. This
  determines whether the data are treated as subject to left truncation
  and right censoring (`"prospective"`) or right truncation alone
  (`"retrospective"`).

- bin_size:

  Size of bins for the discretization of time. A value between 0 and 1
  indicating the size of observed event time quantiles on which to grid
  times (e.g. 0.02 creates a grid of 50 times evenly spaced on the
  quantile scaled). If NULL, defaults to every observed event time.

- time_basis:

  How to treat time for training the binary classifier. Options are
  `"continuous"` and `"dummy"`, meaning an indicator variable is
  included for each time in the time grid.

- learner:

  Which binary regression algorithm to use. Currently, only
  `SuperLearner` is supported, but more learners will be added. See
  below for algorithm-specific arguments.

- SL_control:

  Named list of parameters controlling the Super Learner fitting
  process. These parameters are passed directly to the `SuperLearner`
  function. Parameters include `SL.library` (library of algorithms to
  include in the binary classification Super Learner), `V` (Number of
  cross validation folds on which to train the Super Learner classifier,
  defaults to 10), `method` (Method for estimating coefficients for the
  Super Learner, defaults to `"method.NNLS"`), `stratifyCV` (logical
  indicating whether to stratify by outcome in `SuperLearner`'s
  cross-validation scheme), and `obsWeights` (observation weights,
  passed directly to prediction algorithms by `SuperLearner`).

- tau:

  The maximum time of interest in a study, used for retrospective
  conditional survival estimation. Rather than dealing with right
  truncation separately than left truncation, it is simpler to estimate
  the survival function of `tau - time`. Defaults to `NULL`, in which
  case the maximum study entry time is chosen as the reference point.

## Value

A named list of class `stackL`.

- S_T_preds:

  An `m x k` matrix of estimated event time survival probabilities at
  the `m` covariate vector values and `k` times provided by the user in
  `newX` and `newtimes`, respectively.

- fit:

  The Super Learner fit for binary classification on the stacked
  dataset.

## References

Polley E.C. and van der Laan M.J. (2011). "Super Learning for
Right-Censored Data" in Targeted Learning.

Craig E., Zhong C., and Tibshirani R. (2021). "Survival stacking:
casting survival analysis as a classification problem."

## See also

[predict.stackL](https://cwolock.github.io/survML/reference/predict.stackL.md)
for `stackL` prediction method.

## Examples

``` r
# This is a small simulation example
set.seed(123)
n <- 500
X <- data.frame(X1 = rnorm(n), X2 = rbinom(n, size = 1, prob = 0.5))

S0 <- function(t, x){
  pexp(t, rate = exp(-2 + x[,1] - x[,2] + .5 * x[,1] * x[,2]), lower.tail = FALSE)
}
T <- rexp(n, rate = exp(-2 + X[,1] - X[,2] + .5 *  X[,1] * X[,2]))

G0 <- function(t, x) {
  as.numeric(t < 15) *.9*pexp(t,
                              rate = exp(-2 -.5*x[,1]-.25*x[,2]+.5*x[,1]*x[,2]),
                              lower.tail=FALSE)
}
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

# Note that this a very small Super Learner library, for computational purposes.
SL.library <- c("SL.mean", "SL.glm")

fit <- stackL(time = time,
               event = event,
               entry = entry,
               X = X,
               newX = X,
               newtimes = seq(0, 15, .1),
               direction = "prospective",
               bin_size = 0.1,
               time_basis = "continuous",
               SL_control = list(SL.library = SL.library,
                                 V = 5))

plot(fit$S_T_preds[1,], S0(t =  seq(0, 15, .1), X[1,]))
abline(0,1,col='red')

```
