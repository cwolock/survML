# Generate cross-fitted conditional survival predictions

Generate cross-fitted conditional survival predictions

## Usage

``` r
crossfit_surv_preds(time, event, X, newtimes, folds, pred_generator, ...)
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

- newtimes:

  Numeric vector of times on which to estimate the conditional survival
  functions

- folds:

  `n x 1` numeric vector of folds identifiers for cross-fitting

- pred_generator:

  Function to be used to estimate conditional survival function.

- ...:

  Additional arguments to be passed to `pred_generator`.

## Value

Named list of cross-fitted conditional survival predictions
