% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiseed_vim.R
\name{aggregate_vim}
\alias{aggregate_vim}
\title{Aggregate multiseed VIM results}
\usage{
aggregate_vim(result_list, agg_method, ci_grid, n_eff, alpha = 0.05)
}
\arguments{
\item{result_list}{List of result data frames return by the \code{vim} function.}

\item{agg_method}{P-value aggregation method use to combine results from different seeds. Current options are \code{"bonferroni"}
(Bonferroni's method), \code{"hommel"} (Hommel's method), \code{"arithmetic"} (arithmetic mean), \code{"geometric"} (geometric mean),
\code{"harmonic"} (harmonic mean), \code{"compound_bg"} (compound Bonferroni and geometric mean), and \code{"compound_ba"}
(compound Bonferroni and arithmetic mean). These approaches are discussed at length in Vovk and Wang (2020). Defaults to \code{"compound_bg"}, which has been shown to work well in many settings.}

\item{ci_grid}{Grid of VIM values over which to construct a confidence interval by inverting a hypothesis test. The aggregation works by constructing
hypothesis tests (at level \code{alpha}) of the null corresponding to each value in \code{ci_grid}, and then inverting these tests to yield a
1 - \code{alpha} confidence interval. For example, for \code{"AUC"} importance, the VIM takes values in (0,1), so a grid of values between 0 and 1
would be a reasonable choice.}

\item{n_eff}{The effective sample size. Without sample-splitting, this is simply the sample size. With sample-splitting, this is the sample size divided by
two (i.e., the size of each of the two halves of the data).}

\item{alpha}{The level at which to compute confidence intervals and hypothesis tests. Defaults to 0.05.}
}
\value{
Named list with the following elements:
\item{agg_result}{Data frame giving results aggregated over seeds.}
\item{agg_method}{P-value aggregation method used.}
\item{n_seed}{Number of iterations (seeds) used to perform the VIM estimation procedure.}
\item{vim_objects}{A list of \code{vim} return objects, each corresponding to a different seed.}
}
\description{
Aggregate multiseed VIM results
}
