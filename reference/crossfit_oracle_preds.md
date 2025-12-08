# Generate cross-fitted oracle prediction function estimates

Generate cross-fitted oracle prediction function estimates

## Usage

``` r
crossfit_oracle_preds(
  time,
  event,
  X,
  folds,
  nuisance_preds,
  pred_generator,
  ...
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

- folds:

  `n x 1` numeric vector of folds identifiers for cross-fitting

- nuisance_preds:

  Named list of conditional event and censoring survival functions that
  will be used to estimate the oracle prediction function.

- pred_generator:

  Function to be used to estimate oracle prediction function.

- ...:

  Additional arguments to be passed to `pred_generator`.

## Value

Named list of cross-fitted oracle prediction estimates
