# Estimate variable importance

Compute estimates of and confidence intervals for nonparametric variable
importance based on the difference predictiveness obtained with and
without the feature of interest. Designed for use with time-to-event
outcomes subject to right censoring that may be informed by measured
covariates.

## Usage

``` r
vim(
  type,
  time,
  event,
  X,
  landmark_times = stats::quantile(time[event == 1], probs = c(0.25, 0.5, 0.75)),
  restriction_time = max(time[event == 1]),
  approx_times = NULL,
  large_feature_vector,
  small_feature_vector,
  conditional_surv_generator = NULL,
  conditional_surv_generator_control = NULL,
  large_oracle_generator = NULL,
  large_oracle_generator_control = NULL,
  small_oracle_generator = NULL,
  small_oracle_generator_control = NULL,
  conditional_surv_preds = NULL,
  large_oracle_preds = NULL,
  small_oracle_preds = NULL,
  cf_folds = NULL,
  cf_fold_num = 5,
  sample_split = TRUE,
  ss_folds = NULL,
  robust = TRUE,
  scale_est = FALSE,
  alpha = 0.05,
  verbose = FALSE
)
```

## Arguments

- type:

  Type of VIM to compute. Options include `"accuracy"`, `"AUC"`,
  `"Brier"`, `"R-squared"` `"C-index"`, and `"survival_time_MSE"`.

- time:

  `n x 1` numeric vector of observed follow-up times. If there is
  censoring, these are the minimum of the event and censoring times.

- event:

  `n x 1` numeric vector of status indicators of whether an event was
  observed.

- X:

  `n x p` data.frame of observed covariate values

- landmark_times:

  Numeric vector of length J1 giving landmark times at which to estimate
  VIM (`"accuracy"`, `"AUC"`, `"Brier"`, `"R-squared"`).

- restriction_time:

  Maximum follow-up time for calculation of `"C-index"` and
  `"survival_time_MSE"`. Essentially, this time should be chosen such
  that the conditional survival function is identified at this time for
  all covariate values `X` present in the data. Choosing the restriction
  time such that roughly 10% of individuals remain at-risk at that time
  has been shown to work reasonably well in simulations.

- approx_times:

  Numeric vector of length J2 giving times at which to approximate
  integrals. Defaults to a grid of 100 timepoints, evenly spaced on the
  quantile scale of the distribution of observed event times.

- large_feature_vector:

  Numeric vector giving indices of features to include in the 'large'
  prediction model.

- small_feature_vector:

  Numeric vector giving indices of features to include in the 'small'
  prediction model. Must be a subset of `large_feature_vector`.

- conditional_surv_generator:

  A function to estimate the conditional survival functions of the event
  and censoring variables. Must take arguments (`time`, `event`, `X`)
  (for training purposes) and (`X_holdout` and `newtimes`) (covariate
  values and times at which to generate predictions). Defaults to
  [generate_nuisance_predictions_stackG](https://cwolock.github.io/survML/reference/generate_nuisance_predictions_stackG.md),
  a pre-built generator function based on the
  [stackG](https://cwolock.github.io/survML/reference/stackG.md)
  function. Alternatively, the user can provide their own function for
  this argument, or provide pre-computed estimates to
  `conditional_surv_preds` in lieu of this argument.

- conditional_surv_generator_control:

  A list of arguments to pass to `conditional_surv_generator`.

- large_oracle_generator:

  A function to estimate the oracle prediction function using
  `large_feature_vector`. Must take arguments `time`, `event`, `X`,
  `X_holdout`, and `nuisance_preds`. For all VIM types except for
  `"C-index"`, defaults to
  [generate_oracle_predictions_DR](https://cwolock.github.io/survML/reference/generate_oracle_predictions_DR.md),
  a pre-built generator function using doubly-robust pseudo-outcome
  regression. For `"C-index"`, defaults to
  [generate_oracle_predictions_boost](https://cwolock.github.io/survML/reference/generate_oracle_predictions_boost.md),
  a pre-built generator function using doubly-robust gradient boosting.
  Alternatively, the user can provide their own function, or provide
  pre-computed estimates to `large_oracle_preds` in lieu of this
  argument.

- large_oracle_generator_control:

  A list of arguments to pass to `large_oracle_generator`.

- small_oracle_generator:

  A function to estimate the oracle prediction function using
  `small_feature_vector`. Must take arguments `time`, `event`, `X`,
  `X_holdout`, and `nuisance_preds`. For all VIM types except for
  `"C-index"`, defaults to
  [generate_oracle_predictions_SL](https://cwolock.github.io/survML/reference/generate_oracle_predictions_SL.md),
  a pre-built generator function based on regression the large oracle
  predictions on the small feature vector. For `"C-index"`, defaults to
  [generate_oracle_predictions_boost](https://cwolock.github.io/survML/reference/generate_oracle_predictions_boost.md),
  a pre-built generator function using doubly-robust gradient boosting.
  Alternatively, the user can provide their own function, or provide
  pre-computed estimates to `small_oracle_preds` in lieu of this
  argument.

- small_oracle_generator_control:

  A list of arguments to pass to `small_oracle_generator`.

- conditional_surv_preds:

  User-provided estimates of the conditional survival functions of the
  event and censoring variables given the full covariate vector (if not
  using the `conditional_surv_generator` functionality to compute these
  nuisance estimates). Must be a named list of lists with elements
  `S_hat`, `S_hat_train`, `G_hat`, and `G_hat_train`. If using sample
  splitting, each of these is itself a list of length `2K`, where `K` is
  the number of cross-fitting folds (if not using sample splitting, each
  is a list of length `K`). Each element of these lists is a matrix with
  J2 columns and number of rows equal to either the number of samples in
  the `k`th fold (for `S_hat` and `G_hat`) or the number of samples used
  to compute the nuisance estimates for the `k`th fold (for
  `S_hat_train` and `G_hat_train`).

- large_oracle_preds:

  User-provided estimates of the oracle prediction function using
  `large_feature_vector` (if not using the `large_oracle_generator`
  functionality to compute these nuisance estimates). Must be a named
  list of lists with elements `f0_hat` and `f0_hat_train`. If using
  sample splitting, each of these is itself a list of length `2K` (if
  not using sample splitting, each is a list of length `K`). Each
  element of these lists is a matrix with J1 columns (for landmark time
  VIMs) or 1 column (for `"C-index"` and `"survival_time_MSE"`) and
  number of rows equal to either the number of samples in the `k`th fold
  (for `f0_hat`) or the number of samples used to compute the nuisance
  estimates for the `k`th fold (for `f0_hat_train`).

- small_oracle_preds:

  User-provided estimates of the oracle prediction function using
  `small_feature_vector` (if not using the `small_oracle_generator`
  functionality to compute these nuisance estimates). Must be a named
  list of lists with elements `f0_hat` and `f0_hat_train`. If using
  sample splitting, each of these is itself a list of length `2K` (if
  not using sample splitting, each is a list of length `K`). Each
  element of these lists is a matrix with J1 columns (for landmark time
  VIMs) or 1 column (for `"C-index"` and `"survival_time_MSE"`) and
  number of rows equal to either the number of samples in the `k`th fold
  (for `f0_hat`) or the number of samples used to compute the nuisance
  estimates for the `k`th fold (for `f0_hat_train`).

- cf_folds:

  Numeric vector of length `n` giving cross-fitting folds, if specifying
  the folds explicitly. This is required if you are providing
  pre-computed nuisance estimations — if providing a nuisance generator
  function, the `vim()` will assign folds.

- cf_fold_num:

  The number of cross-fitting folds, if not providing `cf_folds`. Note
  that with samples-splitting, the data will be split into
  `2 x cf_fold_num` folds (i.e., there will be `cf_fold_num` folds
  within each half of the data).

- sample_split:

  Logical indicating whether or not to sample split. Sample-splitting is
  required for valid hypothesis testing of null importance and is
  generally recommended. Defaults to `TRUE`.

- ss_folds:

  Numeric vector of length `n` giving sample-splitting folds, if
  specifying the folds explicitly. This is required if you are providing
  pre-computed nuisance estimations — if providing a nuisance generator
  function, the `vim()` will assign folds.

- robust:

  Logical, whether or not to use the doubly-robust debiasing approach.
  This option is meant for illustration purposes only — it should be
  left as `TRUE`.

- scale_est:

  Logical, whether or not to force the VIM estimate to be nonnegative.

- alpha:

  The level at which to compute confidence intervals and hypothesis
  tests. Defaults to 0.05.

- verbose:

  Whether to print progress messages.

## Value

Named list with the following elements:

- result:

  Data frame giving results. See the documentation of the individual
  `vim_*` functions for details.

- folds:

  A named list giving the cross-fitting fold IDs (`cf_folds`) and
  sample-splitting fold IDs (`ss_folds`).

- approx_times:

  A vector of times used to approximate integrals appearing in the form
  of the VIM estimator.

- conditional_surv_preds:

  A named list containing the estimated conditional event and censoring
  survival functions.

- large_oracle_preds:

  A named list containing the estimated large oracle prediction
  function.

- small_oracle_preds:

  A named list containing the estimated small oracle prediction
  function.

## Details

For nuisance estimation, it is generally advisable to use the pre-built
nuisance generator functions provided by `survML`. See the ”Variable
importance in survival analysis” vignette, or the [package
website](https://cwolock.github.io/survML/) for an illustration.

## References

Wolock C.J., Gilbert P.B., Simon N., and Carone, M. (2025). "Assessing
variable importance in survival analysis using machine learning."

## See also

[vim_accuracy](https://cwolock.github.io/survML/reference/vim_accuracy.md)
[vim_AUC](https://cwolock.github.io/survML/reference/vim_AUC.md)
[vim_brier](https://cwolock.github.io/survML/reference/vim_brier.md)
[vim_cindex](https://cwolock.github.io/survML/reference/vim_cindex.md)
[vim_rsquared](https://cwolock.github.io/survML/reference/vim_rsquared.md)
[vim_survival_time_mse](https://cwolock.github.io/survML/reference/vim_survival_time_mse.md)

## Examples

``` r
# This is a small simulation example
set.seed(123)
n <- 100
X <- data.frame(X1 = rnorm(n), X2 = rbinom(n, size = 1, prob = 0.5))

T <- rexp(n, rate = exp(-2 + X[,1] - X[,2] + .5 *  X[,1] * X[,2]))

C <- rexp(n, exp(-2 -.5 * X[,1] - .25 * X[,2] + .5 * X[,1] * X[,2]))
C[C > 15] <- 15

time <- pmin(T, C)
event <- as.numeric(T <= C)

# landmark times for AUC
landmark_times <- c(3)

output <- vim(type = "AUC",
              time = time,
              event = event,
              X = X,
              landmark_times = landmark_times,
              large_feature_vector = 1:2,
              small_feature_vector = 2,
              conditional_surv_generator_control = list(SL.library = c("SL.mean", "SL.glm")),
              large_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm")),
              small_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm")),
              cf_fold_num = 2,
              sample_split = FALSE,
              scale_est = TRUE)

print(output$result)
#>   landmark_time       est  var_est        cil       ciu cil_1sided  p
#> 1             3 0.2823303 1.407984 0.04976388 0.5148967 0.08715441 NA
#>   large_predictiveness small_predictiveness vim large_feature_vector
#> 1            0.8209323             0.538602 AUC                  1,2
#>   small_feature_vector
#> 1                    2
```
