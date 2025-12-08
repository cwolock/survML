# Estimate conditional survival function nuisance parameters using survival stacking

Estimate conditional survival function nuisance parameters using
survival stacking

## Usage

``` r
generate_nuisance_predictions_stackG(
  time,
  event,
  X,
  X_holdout,
  newtimes,
  SL.library = c("SL.mean", "SL.glm", "SL.earth", "SL.gam", "SL.ranger"),
  V = 5,
  bin_size = 0.05,
  approx_times
)
```

## Arguments

- time:

  `n x 1` numeric vector of observed follow-up times. If there is
  censoring, these are the minimum of the event and censoring times.

- event:

  `n x 1` numeric vector of status indicators of whether an event was
  observed.

- X:

  `n x p` data.frame of observed covariate values

- X_holdout:

  `m x p` data.frame of new observed covariate values at which to obtain
  `m` predictions for the estimated algorithm. Must have the same names
  and structure as `X`.

- newtimes:

  `k x 1` numeric vector of times at which to obtain `k` predicted
  conditional survivals.

- SL.library:

  Super Learner library

- V:

  Number of cross-validation folds, to be passed to `SuperLearner`

- bin_size:

  Size of time bin on which to discretize for estimation of cumulative
  probability functions. Can be a number between 0 and 1, indicating the
  size of quantile grid (e.g. `0.1` estimates the cumulative probability
  functions on a grid based on deciles of observed `time`s). If `NULL`,
  creates a grid of all observed `time`s. See
  [stackG](https://cwolock.github.io/survML/reference/stackG.md)
  documentation.

- approx_times:

  Numeric vector of times at which to approximate product integral or
  cumulative hazard interval. See
  [stackG](https://cwolock.github.io/survML/reference/stackG.md)
  documentation.

## Value

A list containing elements `S_hat` (conditional event survival function,
corresponding to `X_holdout` and `newtimes`), `S_hat_train` (conditional
event survival function, corresponding to `X` and `newtimes`), `G_hat`
(conditional censoring survival function, corresponding to `X_holdout`
and `newtimes`), and `G_hat_train` (conditional censoring survival
function, corresponding to `X` and `newtimes`)

## See also

[stackG](https://cwolock.github.io/survML/reference/stackG.md)
