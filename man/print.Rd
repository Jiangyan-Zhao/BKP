% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_BKP.R, R/print_DKP.R
\name{print}
\alias{print}
\alias{print.BKP}
\alias{print.DKP}
\title{Print Summary of a Fitted BKP or DKP Model}
\usage{
\method{print}{BKP}(x, ...)

\method{print}{DKP}(x, ...)
}
\arguments{
\item{x}{An object of class \code{"BKP"} (from \code{\link{fit.BKP}}) or
\code{"DKP"} (from \code{\link{fit.DKP}}).}

\item{...}{Additional arguments passed to the generic \code{print} method
(currently not used).}
}
\value{
Invisibly returns the input object (of class \code{"BKP"} or
\code{"DKP"}). The function is called for its side effect of printing a
summary to the console.
}
\description{
Displays a concise summary of a fitted BKP or DKP model. The
output includes key characteristics such as sample size, input
dimensionality, kernel type, loss function, optimized kernel
hyperparameters, and minimum loss.
}
\examples{
# ============================================================== #
# ========================= BKP Examples ======================= #
# ============================================================== #

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

# ============================================================== #
# ========================= DKP Examples ======================= #
# ============================================================== #

#-------------------------- 1D Example ---------------------------
set.seed(123)

# Define true class probability function (3-class)
true_pi_fun <- function(X) {
  p <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
  return(matrix(c(p/2, p/2, 1 - p), nrow = length(p)))
}

n <- 30
Xbounds <- matrix(c(-2, 2), nrow = 1)
X <- tgp::lhs(n = n, rect = Xbounds)
true_pi <- true_pi_fun(X)
m <- sample(100, n, replace = TRUE)

# Generate multinomial responses
Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))

# Fit DKP model
model1 <- fit.DKP(X, Y, Xbounds = Xbounds)
print(model1)


#-------------------------- 2D Example ---------------------------
set.seed(123)

# Define latent function and transform to 3-class probabilities
true_pi_fun <- function(X) {
  if (is.null(nrow(X))) X <- matrix(X, nrow = 1)
  m <- 8.6928; s <- 2.4269
  x1 <- 4 * X[,1] - 2
  x2 <- 4 * X[,2] - 2
  a <- 1 + (x1 + x2 + 1)^2 *
    (19 - 14*x1 + 3*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2)
  b <- 30 + (2*x1 - 3*x2)^2 *
    (18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2)
  f <- (log(a * b) - m) / s
  p <- pnorm(f)
  return(matrix(c(p/2, p/2, 1 - p), nrow = length(p)))
}

n <- 100
Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
X <- tgp::lhs(n = n, rect = Xbounds)
true_pi <- true_pi_fun(X)
m <- sample(100, n, replace = TRUE)

# Generate multinomial responses
Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))

# Fit DKP model
model2 <- fit.DKP(X, Y, Xbounds = Xbounds)
print(model2)

}
\seealso{
\code{\link{fit.BKP}}, \code{\link{fit.DKP}},
\code{\link{summary.BKP}}, \code{\link{summary.DKP}}.
}
\keyword{BKP}
\keyword{DKP}
