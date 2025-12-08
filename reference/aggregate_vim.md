# Aggregate multiseed VIM results

Aggregate multiseed VIM results

## Usage

``` r
aggregate_vim(result_list, agg_method, ci_grid, n_eff, alpha = 0.05)
```

## Arguments

- result_list:

  List of result data frames return by the `vim` function.

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

- n_eff:

  The effective sample size. Without sample-splitting, this is simply
  the sample size. With sample-splitting, this is the sample size
  divided by two (i.e., the size of each of the two halves of the data).

- alpha:

  The level at which to compute confidence intervals and hypothesis
  tests. Defaults to 0.05.

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
