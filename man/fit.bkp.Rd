% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_BKP.R
\name{fit.BKP}
\alias{fit.BKP}
\title{Fit a Beta Kernel Process (BKP) Model}
\usage{
fit.BKP(
  X,
  y,
  m,
  Xbounds = NULL,
  prior = c("noninformative", "fixed", "adaptive"),
  r0 = 2,
  p0 = 0.5,
  kernel = c("gaussian", "matern52", "matern32"),
  loss = c("brier", "log_loss"),
  n_multi_start = NULL
)
}
\arguments{
\item{X}{A numeric input matrix of size \eqn{n \times d}, where each row
represents a covariate vector.}

\item{y}{A numeric vector of observed successes (length \code{n}).}

\item{m}{A numeric vector of total binomial trials (length \code{n}),
corresponding to each \code{y}.}

\item{Xbounds}{Optional \eqn{d \times 2} matrix specifying the lower and
upper bounds of each input dimension. Used to normalize inputs to
\eqn{[0,1]^d}. If \code{Xbounds} is \code{NULL}, the input is assumed to
have already been normalized, and the default bounds are set to \eqn{[0,1]^d}.}

\item{prior}{Type of prior to use. One of \code{"noninformative"},
\code{"fixed"}, or \code{"adaptive"}.}

\item{r0}{Global prior precision (only used when \code{prior = "fixed"} or
\code{"adaptive"}).}

\item{p0}{Global prior mean (only used when \code{prior = "fixed"}).}

\item{kernel}{Kernel function for local weighting. Choose from
\code{"gaussian"}, \code{"matern52"}, or \code{"matern32"}.}

\item{loss}{Loss function for kernel hyperparameter tuning. One of
\code{"brier"} (default) or \code{"log_loss"}.}

\item{n_multi_start}{Number of random initializations for multi-start
optimization. Default is \code{10 × d}.}
}
\value{
A list of class \code{"BKP"} containing the fitted BKP model, with
the following elements:
\describe{
\item{\code{theta_opt}}{Optimized kernel hyperparameters (lengthscales).}
\item{\code{kernel}}{Kernel function used, as a string.}
\item{\code{loss}}{Loss function used for hyperparameter tuning.}
\item{\code{loss_min}}{Minimum loss value achieved during optimization.}
\item{\code{X}}{Original (unnormalized) input matrix of size \code{n × d}.}
\item{\code{Xnorm}}{Normalized input matrix scaled to \eqn{[0,1]^d}.}
\item{\code{Xbounds}}{Matrix specifying normalization bounds for each input dimension.}
\item{\code{y}}{Observed success counts.}
\item{\code{m}}{Observed binomial trial counts.}
\item{\code{prior}}{Type of prior used.}
\item{\code{r0}}{Prior precision parameter.}
\item{\code{p0}}{Prior mean (for fixed priors).}
\item{\code{alpha0}}{Prior shape parameter \eqn{\alpha_0(\mathbf{x})}, either a scalar or vector.}
\item{\code{beta0}}{Prior shape parameter \eqn{\beta_0(\mathbf{x})}, either a scalar or vector.}
\item{\code{alpha_n}}{Posterior shape parameter \eqn{\alpha_n(\mathbf{x})}.}
\item{\code{beta_n}}{Posterior shape parameter \eqn{\beta_n(\mathbf{x})}.}
}
}
\description{
Fits a BKP model to binomial or binary response data via local
kernel smoothing. The model constructs a flexible latent probability
surface by updating Beta priors using kernel-weighted observations.
}
\examples{
#-------------------------- 1D Example ---------------------------
set.seed(123)

# Define true success probability function
true_pi_fun <- function(x) {
  (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
}

n <- 30
Xbounds <- matrix(c(-2,2), nrow=1)
X <- tgp::lhs(n = n, rect = Xbounds)
true_pi <- true_pi_fun(X)
m <- sample(100, n, replace = TRUE)
y <- rbinom(n, size = m, prob = true_pi)

# Fit BKP model
model1 <- fit.BKP(X, y, m, Xbounds=Xbounds)
print(model1)


#-------------------------- 2D Example ---------------------------
set.seed(123)

# Define 2D latent function and probability transformation
true_pi_fun <- function(X) {
  if(is.null(nrow(X))) X <- matrix(X, nrow=1)
  m <- 8.6928
  s <- 2.4269
  x1 <- 4*X[,1]- 2
  x2 <- 4*X[,2]- 2
  a <- 1 + (x1 + x2 + 1)^2 *
    (19- 14*x1 + 3*x1^2- 14*x2 + 6*x1*x2 + 3*x2^2)
  b <- 30 + (2*x1- 3*x2)^2 *
    (18- 32*x1 + 12*x1^2 + 48*x2- 36*x1*x2 + 27*x2^2)
  f <- log(a*b)
  f <- (f- m)/s
  return(pnorm(f))  # Transform to probability
}

n <- 100
Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
X <- tgp::lhs(n = n, rect = Xbounds)
true_pi <- true_pi_fun(X)
m <- sample(100, n, replace = TRUE)
y <- rbinom(n, size = m, prob = true_pi)

# Fit BKP model
model2 <- fit.BKP(X, y, m, Xbounds=Xbounds)
print(model2)

}
\seealso{
\code{\link{fit.DKP}} for modeling multinomial responses using the
Dirichlet Kernel Process. \code{\link{predict.BKP}},
\code{\link{plot.BKP}}, \code{\link{simulate.BKP}} for making predictions,
visualizing results, and generating simulations from a fitted BKP model.
\code{\link{summary.BKP}}, \code{\link{print.BKP}} for inspecting model
details.
}
