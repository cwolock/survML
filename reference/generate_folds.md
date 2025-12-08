# Generate cross-fitting and sample-splitting folds

Generate cross-fitting and sample-splitting folds

## Usage

``` r
generate_folds(n, V, sample_split, cf_folds = NULL)
```

## Arguments

- n:

  Total sample size

- V:

  Number of cross-fitting folds to use

- sample_split:

  Logical, whether or not sample-splitting is being used

- cf_folds:

  Cross-fitting folds, if already provided by user.

## Value

Named list of cross-fitting and sample-splitting folds
