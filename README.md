
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `survML`: Conditional survival function estimation using machine learning

The `survML` package implements two methods for estimating a conditional
survival function using off-the-shelf machine learning. The first,
called cumulative probability estimation and performed using the
`survMLc()` function, involves decomposing the conditional hazard
function into regression functions depending only on the observed data.
The second, called survival stacking or discrete-time hazard estimation,
involves discretizing time and estimating the probability of an event of
interest within each discrete time period. More details on each method,
as well as examples, follow.

# Installing `survML`

**Note: `survML` is currently under development.**

The development version of `survML` is available on GitHub. You can
install it using the `devtools` package as follows:

``` r
## install.packages("devtools") # run only if necessary
install_github(repo = "cwolock/survML")
```

# Cumulative probability estimation

In a basic survival analysis setting with right-censored data (but no
truncation), the cumulative probability estimator requires three pieces:
(1) the conditional probability of censoring given covariates, (2) the
CDF of the observed time given covariates, among censored subjects, and
(3) the CDF of the observed time given covariates, among uncensored
subjects. All three of these can be estimated using standard binary
regression or classification methods.

For maximum flexibility, `survML` uses Super Learner for binary
regression. Estimating (1) is a standard binary regression problem. We
use pooled binary regression to estimate (2) and (3). In essence, at
time
![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t "t")
each on a user-specified grid, the CDF is a binary regression using the
outcome
![I(Y \\leq t)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;I%28Y%20%5Cleq%20t%29 "I(Y \leq t)").
The data sets for each
![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t "t")
are combined into a single, pooled data set, including
![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t "t")
as a covariate.

The `survMLc` function performs cumulative probability estimation. The
most important user-specified arguments are described below:

-   `bin_size`: This is the size of time grid used for estimating (2)
    and (3). In most cases, a finer grid performs better than a coarser
    grid, at increased computational cost. We recommend using as fine a
    grid as computational resources and time allow. In simulations, a
    grid of 40 time points performed similarly to a grid of every
    observed time. Bin size is given in quantile terms;
    `bin_size = 0.025` will use times corresponding to quantiles
    ![\\{0, 0.025, 0.05, \\dots, 0.975, 1\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%7B0%2C%200.025%2C%200.05%2C%20%5Cdots%2C%200.975%2C%201%5C%7D "\{0, 0.025, 0.05, \dots, 0.975, 1\}").
    If `NULL`, a grid of every observed time is used.
-   `time_basis`: This is how the time variable
    ![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t "t")
    is included in the pooled data set. The default is `continuous`. It
    is also possible to include a dummy variable for each time in the
    grid using option `dummy`.
-   `SL.library`

# Survival stacking

# Example

## Subtitle for example

# References

For details of this method, please see the following preprint:

Charles J. Wolock, Peter B. Gilbert, Noah Simon, Marco Carone. “Flexible
estimation of the conditional survival function via observable
regression models.” (2022).

Super Learner citation.
