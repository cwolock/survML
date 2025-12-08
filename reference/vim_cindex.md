# Estimate concordance index VIM

Estimate concordance index VIM

## Usage

``` r
vim_cindex(
  time,
  event,
  approx_times,
  restriction_time,
  f_hat,
  fs_hat,
  S_hat,
  G_hat,
  cf_folds,
  sample_split,
  ss_folds,
  scale_est = FALSE,
  alpha = 0.05
)
```

## Arguments

- time:

  `n x 1` numeric vector of observed follow-up times If there is
  censoring, these are the minimum of the event and censoring times.

- event:

  `n x 1` numeric vector of status indicators of whether an event was
  observed. Defaults to a vector of 1s, i.e. no censoring.

- approx_times:

  Numeric vector of length J1 giving times at which to approximate
  integrals.

- restriction_time:

  Restriction time (upper bound for event times to be compared in
  computing the C-index)

- f_hat:

  Full oracle predictions (n x J1 matrix)

- fs_hat:

  Residual oracle predictions (n x J1 matrix)

- S_hat:

  Estimates of conditional event time survival function (n x J2 matrix)

- G_hat:

  Estimate of conditional censoring time survival function (n x J2
  matrix)

- cf_folds:

  Numeric vector of length n giving cross-fitting folds

- sample_split:

  Logical indicating whether or not to sample split

- ss_folds:

  Numeric vector of length n giving sample-splitting folds

- scale_est:

  Logical, whether or not to force the VIM estimate to be nonnegative

- alpha:

  The level at which to compute confidence intervals and hypothesis
  tests. Defaults to 0.05

## Value

A data frame giving results, with the following columns:

- restriction_time:

  Restriction time (upper bound for event times to be compared in
  computing the C-index).

- est:

  VIM point estimate.

- var_est:

  Estimated variance of the VIM estimate.

- cil:

  Lower bound of the VIM confidence interval.

- ciu:

  Upper bound of the VIM confidence interval.

- cil_1sided:

  Lower bound of a one-sided confidence interval.

- p:

  p-value corresponding to a hypothesis test of null importance.

- large_predictiveness:

  Estimated predictiveness of the large oracle prediction function.

- small_predictiveness:

  Estimated predictiveness of the small oracle prediction function.

- vim:

  VIM type.

- large_feature_vector:

  Group of features available for the large oracle prediction function.

- small_feature_vector:

  Group of features available for the small oracle prediction function.

## See also

[vim](https://cwolock.github.io/survML/reference/vim.md) for example
usage
