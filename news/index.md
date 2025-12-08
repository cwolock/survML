# Changelog

## survML (development version)

## survML 1.2.0

CRAN release: 2024-10-30

## survML 1.2.0

CRAN release: 2024-10-30

- Added variable importance and current status isotonic regression
  functions; these changes do not affect the existing `stackG` and
  `stackL` functions.
- Fixed bug where `stackG` threw an error if there was no censoring.
- Added cumulative hazard estimates to `stackG` output.

## survML 1.1.0

CRAN release: 2024-03-17

- Added `gam` to `SUGGESTS` in order to allow `SuperLearner` package to
  make corresponding change without breaking vignettes.
- Added `time_grid_fit` option to main `stackG` function in order to
  allow more flexibility in choosing time grids.
- Minor bug fixes.

## survML 1.0.0

CRAN release: 2023-07-08

- Initial CRAN submission.
