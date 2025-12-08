# Estimate full oracle prediction function using DR pseudo-outcome regression

Estimate full oracle prediction function using DR pseudo-outcome
regression

## Usage

``` r
generate_oracle_predictions_DR(
  time,
  event,
  X,
  X_holdout,
  nuisance_preds,
  outcome,
  landmark_times,
  restriction_time,
  approx_times,
  SL.library = c("SL.mean", "SL.glm", "SL.earth", "SL.gam", "SL.ranger"),
  V = 5,
  indx
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

- outcome:

  Outcome type, either `"survival_probability"` or
  `"restricted_survival_time"`

- landmark_times:

  Numeric vector of length J1 giving landmark times at which to estimate
  VIM (`"accuracy"`, `"AUC"`, `"Brier"`, `"R-squared"`).

- restriction_time:

  Maximum follow-up time for calculation of `"survival_time_MSE"`.
  Essentially, this time should be chosen such that the conditional
  survival function is identified at this time for all covariate values
  `X` present in the data. Choosing the restriction time such that
  roughly 10% of individuals remain at-risk at that time has been shown
  to work reasonably well in simulations.

- approx_times:

  Numeric vector of length J2 giving times at which to approximate
  integral appearing in the pseudo-outcomes

- SL.library:

  Super Learner library

- V:

  Number of cross-validation folds, to be passed to `SuperLearner`

- indx:

  Numeric index of column(s) of `X` to be removed, i.e., not used in the
  oracle prediction function.

## Value

A list containing elements `f0_hat` and `f0_hat_train`, the estimated
oracle prediction functions for `X_holdout` and `X`, respectively.

## See also

[DR_pseudo_outcome_regression](https://cwolock.github.io/survML/reference/DR_pseudo_outcome_regression.md)
