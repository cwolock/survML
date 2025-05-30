% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crossfit_surv_preds.R
\name{crossfit_surv_preds}
\alias{crossfit_surv_preds}
\title{Generate cross-fitted conditional survival predictions}
\usage{
crossfit_surv_preds(time, event, X, newtimes, folds, pred_generator, ...)
}
\arguments{
\item{time}{\code{n x 1} numeric vector of observed
follow-up times. If there is censoring, these are the minimum of the
event and censoring times.}

\item{event}{\code{n x 1} numeric vector of status indicators of
whether an event was observed.}

\item{X}{\code{n x p} data.frame of observed covariate values}

\item{newtimes}{Numeric vector of times on which to estimate the conditional survival functions}

\item{folds}{\code{n x 1} numeric vector of folds identifiers for cross-fitting}

\item{pred_generator}{Function to be used to estimate conditional survival function.}

\item{...}{Additional arguments to be passed to \code{pred_generator}.}
}
\value{
Named list of cross-fitted conditional survival predictions
}
\description{
Generate cross-fitted conditional survival predictions
}
