
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start  -->

![GitHub R package
version](https://img.shields.io/github/r-package/v/cwolock/survML)
[![CRAN
status](https://www.r-pkg.org/badges/version/survML)](https://CRAN.R-project.org/package=survML)
![GitHub](https://img.shields.io/github/license/cwolock/survML)
[![R-CMD-check](https://github.com/cwolock/survML/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/cwolock/survML/actions/workflows/R-CMD-check.yml)
[![Codecov test
coverage](https://codecov.io/gh/cwolock/survML/branch/main/graph/badge.svg)](https://app.codecov.io/gh/cwolock/survML?branch=main)

<!-- badges: end -->

# `survML`: Tools for Flexible Survival Analysis Using Machine Learning

**Note: The current development version of `survML` now has
functionality for estimating variable importance and for estimating a
covariate-adjusted survival curve under current status sampling, in
addition to the original survival stacking functionality that was
included in versions 1.1.0 and earlier. A new version on CRAN is
forthcoming.**

The `survML` package contains a variety of functions for analyzing
survival data using machine learning. These include:

1.  *Global and local survival stacking*: Use off-the-shelf machine
    learning tools to estimate conditional survival functions.

2.  *Algorithm-agnostic variable importance*: Use debiased machine
    learning to estimate and make inference on variable importance for
    prediction of time-to-event outcomes.

3.  *Current-status isotonic regression*: Use isotonic regression to
    estimate the covariate-adjusted survival function of a time-to-event
    outcome under current status sampling.

See the package vignettes and function reference for more details.

## Installing `survML`

You can install a stable version of `survML` from CRAN using

``` r
install.packages("survML")
```

Alternatively, the development version of `survML` is available on
GitHub. You can install it using the `devtools` package as follows:

``` r
## install.packages("devtools") # run only if necessary
install_github(repo = "cwolock/survML")
```

## Documentation

Full documentation can be found on the `survML` website at
<https://cwolock.github.io/survML/>.

## Bugs reports and feature requests

To submit a bug report or request a new feature, please submit a new
[GitHub Issue](https://github.com/cwolock/survML/issues).

## References

For details of global survival stacking, please see the following paper:

- Charles J. Wolock, Peter B. Gilbert, Noah Simon and Marco Carone. [“A
  framework for leveraging machine learning tools to estimate
  personalized survival
  curves.”](https://www.tandfonline.com/doi/full/10.1080/10618600.2024.2304070)
  *Journal of Computational and Graphical Statistics* (2024).

The following preprints contain details on assessing variable importance
in survival analysis, and using isotonic regression to estimate survival
curves from current status data, respectively:

- Charles J. Wolock, Peter B. Gilbert, Noah Simon and Marco Carone.
  “Nonparametric variable importance for time-to-event outcomes with
  application to prediction of HIV infection.”
  [arXiv:2311.12726.](https://arxiv.org/abs/2311.12726)

- Charles J. Wolock, Susan Jacob, Julia C. Bennett, Anna Elias-Warren,
  Jessica O’Hanlon, Avi Kenny, Nicholas P. Jewell, Andrea Rotnitzky,
  Ana A. Weil, Helen Y. Chu and Marco Carone. “Investigating symptom
  duration using current status data: a case study of post-acute
  COVID-19 syndrome.”
  [arXiv:2407.04214.](https://arxiv.org/abs/2407.04214)

Local survival stacking is described in:

- Eric C. Polley and Mark J. van der Laan. “Super Learning for
  Right-Censored Data” in *Targeted Learning* (2011).

- Erin Craig, Chenyang Zhong, and Robert Tibshirani. “Survival stacking:
  casting survival analysis as a classification problem.”
  [arXiv:2107.13480.](https://arxiv.org/abs/2107.13480)

## Citation

After using the `survML` package for conditional survival estimation,
please cite the following:

    @article{wolock2024framework,
            title={A framework for leveraging machine learning tools to estimate personalized survival curves},
            author={Wolock, Charles J and Gilbert, Peter B and Simon, Noah and Carone, Marco},
            journal={Journal of Computational and Graphical Statistics},
            year={2024},
            publisher={Taylor \& Francis},
            doi={10.1080/10618600.2024.2304070}
    }

After using the variable importance functions, please cite the
following:

    @article{wolock2023assessing,
             title={Assessing variable importance in survival analysis using machine learning},
             author={Wolock, Charles J and Gilbert, Peter B and Simon, Noah and Carone, Marco},
             journal={arXiv preprint arXiv:2311.12726},
             year={2023}
    }

After using the functionality for current status data, please cite the
following:

    @article{wolock2024investigating,
      title={Investigating symptom duration using current status data: a case study of post-acute COVID-19 syndrome},
      author={Wolock, Charles J and Jacob, Susan and Bennett, Julia C and Elias-Warren, Anna and O'Hanlon, Jessica and Kenny, Avi and Jewell, Nicholas P and Rotnitzky, Andrea and Weil, Ana A and Chu, Helen Y and Carone, Marco},
      journal={arXiv preprint arXiv:2407.04214},
      year={2024}
    }
