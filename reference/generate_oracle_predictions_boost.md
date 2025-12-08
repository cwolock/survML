# Estimate oracle prediction function using DR gradient boosting

Estimate oracle prediction function using DR gradient boosting

## Usage

``` r
generate_oracle_predictions_boost(
  time,
  event,
  X,
  X_holdout,
  nuisance_preds,
  restriction_time,
  approx_times,
  V = 5,
  indx,
  tuning = FALSE,
  subsample_n = length(time),
  boosting_params = list(mstop = c(100), nu = c(0.1), sigma = c(0.01), learner =
    c("glm"))
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

- nuisance_preds:

  Named list of conditional survival function predictions with elements
  `"S_hat"`, `"S_hat_train"`, `"G_hat"`, and `"G_hat_train"`. This
  should match the output of `conditional_surv_generator`.

- restriction_time:

  Maximum follow-up time for calculation of the C-index. Essentially,
  this time should be chosen such that the conditional survival function
  is identified at this time for all covariate values `X` present in the
  data. Choosing the restriction time such that roughly 10% of
  individuals remain at-risk at that time has been shown to work
  reasonably well in simulations.

- approx_times:

  Numeric vector of length J2 giving times at which to approximate
  C-index integral.

- V:

  Number of cross-validation folds for selection of tuning parameters

- indx:

  Numeric index of column(s) of `X` to be removed, i.e., not used in the
  oracle prediction function.

- tuning:

  Logical, whether or not to use cross-validation to select tuning
  parameters

- subsample_n:

  Number of samples to use for boosting procedure. Using a subsample of
  the full sample can greatly reduce runtime

- boosting_params:

  Named list of parameter values for the boosting procedure. Elements of
  this list include `mstop` (number of boosting iterations), `nu`
  (learning rate), `sigma` (smoothness parameter for sigmoid
  approximation, with smaller meaning less smoothing), and `learner`
  (base learner, can take values `"glm"`, `"gam"`, or `"tree"`)

## Value

A list containing elements `f0_hat` and `f0_hat_train`, the estimated
oracle prediction functions for `X_holdout` and `X`, respectively.

## See also

[boost_c_index](https://cwolock.github.io/survML/reference/boost_c_index.md)
