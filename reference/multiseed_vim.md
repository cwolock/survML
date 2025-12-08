# Estimate variable importance with multiple seeds

Repeat the VIM estimation procedure multiple times and aggregate the
results, mitigating the additional randomness introduced by
sample-splitting and cross-fitting.

## Usage

``` r
multiseed_vim(
  n_seed,
  agg_method = "compound_bg",
  ci_grid,
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
  cf_fold_num = 5,
  sample_split = TRUE,
  scale_est = FALSE,
  alpha = 0.05,
  verbose = FALSE
)
```

## Arguments

- n_seed:

  Number of iterations (seeds) to perform the VIM estimation procedure.
  These will be aggregated into a single result.

- agg_method:

  P-value aggregation method use to combine results from different
  seeds. Current options are `"bonferroni"` (Bonferroni's method),
  `"hommel"` (Hommel's method), `"arithmetic"` (arithmetic mean),
  `"geometric"` (geometric mean), `"harmonic"` (harmonic mean),
  `"compound_bg"` (compound Bonferroni and geometric mean), and
  `"compound_ba"` (compound Bonferroni and arithmetic mean). These
  approaches are discussed at length in Vovk and Wang (2020). Defaults
  to `"compound_bg"`, which has been shown to work well in many
  settings.

- ci_grid:

  Grid of VIM values over which to construct a confidence interval by
  inverting a hypothesis test. The aggregation works by constructing
  hypothesis tests (at level `alpha`) of the null corresponding to each
  value in `ci_grid`, and then inverting these tests to yield a 1 -
  `alpha` confidence interval. For example, for `"AUC"` importance, the
  VIM takes values in (0,1), so a grid of values between 0 and 1 would
  be a reasonable choice.

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

- cf_fold_num:

  The number of cross-fitting folds, if not providing `cf_folds`. Note
  that with samples-splitting, the data will be split into
  `2 x cf_fold_num` folds (i.e., there will be `cf_fold_num` folds
  within each half of the data).

- sample_split:

  Logical indicating whether or not to sample split. Sample-splitting is
  required for valid hypothesis testing of null importance and is
  generally recommended. Defaults to `TRUE`.

- scale_est:

  Logical, whether or not to force the VIM estimate to be nonnegative.

- alpha:

  The level at which to compute confidence intervals and hypothesis
  tests. Defaults to 0.05.

- verbose:

  Whether to print progress messages.

## Value

Named list with the following elements:

- agg_result:

  Data frame giving results aggregated over seeds.

- agg_method:

  P-value aggregation method used.

- n_seed:

  Number of iterations (seeds) used to perform the VIM estimation
  procedure.

- vim_objects:

  A list of `vim` return objects, each corresponding to a different
  seed.

## Details

Using a larger value of `n_seed` will result in more stable results, at
a greater computational cost.

## References

Vovk V. and Wang R. (2020). "Combining p-values via averaging."

Wolock C.J., Gilbert P.B., Simon N., and Carone, M. (2025). "Assessing
variable importance in survival analysis using machine learning."

## See also

[vim](https://cwolock.github.io/survML/reference/vim.md)

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

output <- multiseed_vim(n_seed = 2,
              agg_method = "compound_bg",
              ci_grid = seq(0, 1, by = 0.01),
              type = "AUC",
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
              sample_split = TRUE,
              scale_est = TRUE)

print(output$result)
#> NULL
```
