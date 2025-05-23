% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stackG.R
\name{predict.stackG}
\alias{predict.stackG}
\title{Obtain predicted conditional survival and cumulative hazard functions from a global survival stacking object}
\usage{
\method{predict}{stackG}(
  object,
  newX,
  newtimes,
  surv_form = object$surv_form,
  time_grid_approx = object$time_grid_approx,
  ...
)
}
\arguments{
\item{object}{Object of class \code{stackG}}

\item{newX}{\code{m x p} data.frame of new observed covariate
values at which to obtain \code{m} predictions for the estimated algorithm.
Must have the same names and structure as \code{X}.}

\item{newtimes}{\code{k x 1} numeric vector of times at which to obtain \code{k}
predicted conditional survivals.}

\item{surv_form}{Mapping from hazard estimate to survival estimate.
Can be either \code{"PI"} (product integral mapping) or \code{"exp"}
(exponentiated cumulative hazard estimate). Defaults to the value
saved in \code{object}.}

\item{time_grid_approx}{Numeric vector of times at which to
approximate product integral or cumulative hazard interval. Defaults to the value
saved in \code{object}.}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
A named list with the following components:
\item{S_T_preds}{An \code{m x k} matrix of estimated event time survival probabilities at the
\code{m} covariate vector values and \code{k} times provided by the user in
\code{newX} and \code{newtimes}, respectively.}
\item{S_C_preds}{An \code{m x k} matrix of estimated censoring time survival probabilities at the
\code{m} covariate vector values and \code{k} times provided by the user in
\code{newX} and \code{newtimes}, respectively.}
\item{Lambda_T_preds}{An \code{m x k} matrix of estimated event time cumulative hazard function values at the
\code{m} covariate vector values and \code{k} times provided by the user in
\code{newX} and \code{newtimes}, respectively.}
\item{Lambda_C_preds}{An \code{m x k} matrix of estimated censoring time cumulative hazard function values at the
\code{m} covariate vector values and \code{k} times provided by the user in
\code{newX} and \code{newtimes}, respectively.}
\item{time_grid_approx}{The approximation grid for the product integral or cumulative hazard integral,
(user-specified).}
\item{surv_form}{Exponential or product-integral form (user-specified).}
}
\description{
Obtain predicted conditional survival and cumulative hazard functions from a global survival stacking object
}
\examples{

# This is a small simulation example
set.seed(123)
n <- 250
X <- data.frame(X1 = rnorm(n), X2 = rbinom(n, size = 1, prob = 0.5))

S0 <- function(t, x){
  pexp(t, rate = exp(-2 + x[,1] - x[,2] + .5 * x[,1] * x[,2]), lower.tail = FALSE)
}
T <- rexp(n, rate = exp(-2 + X[,1] - X[,2] + .5 *  X[,1] * X[,2]))

G0 <- function(t, x) {
  as.numeric(t < 15) *.9*pexp(t,
                              rate = exp(-2 -.5*x[,1]-.25*x[,2]+.5*x[,1]*x[,2]),
                              lower.tail=FALSE)
}
C <- rexp(n, exp(-2 -.5 * X[,1] - .25 * X[,2] + .5 * X[,1] * X[,2]))
C[C > 15] <- 15

entry <- runif(n, 0, 15)

time <- pmin(T, C)
event <- as.numeric(T <= C)

sampled <- which(time >= entry)
X <- X[sampled,]
time <- time[sampled]
event <- event[sampled]
entry <- entry[sampled]

# Note that this a very small Super Learner library, for computational purposes.
SL.library <- c("SL.mean", "SL.glm")

fit <- stackG(time = time,
              event = event,
              entry = entry,
              X = X,
              newX = X,
              newtimes = seq(0, 15, .1),
              direction = "prospective",
              bin_size = 0.1,
              time_basis = "continuous",
              time_grid_approx = sort(unique(time)),
              surv_form = "exp",
              learner = "SuperLearner",
              SL_control = list(SL.library = SL.library,
                                V = 5))

preds <- predict(object = fit,
                 newX = X,
                 newtimes = seq(0, 15, 0.1))

plot(preds$S_T_preds[1,], S0(t =  seq(0, 15, .1), X[1,]))
abline(0,1,col='red')
}
\seealso{
\link{stackG}
}
