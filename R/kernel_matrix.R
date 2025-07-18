#' @title Compute Kernel Matrix Between Input Locations
#'
#' @description Computes the kernel matrix between two sets of input locations
#'   using a specified kernel function. Supports both isotropic and anisotropic
#'   lengthscales. Available kernels include the Gaussian, Matérn 5/2, and
#'   Matérn 3/2.
#'
#' @param X A numeric matrix (or vector) of input locations with shape \eqn{n
#'   \times d}.
#' @param Xprime An optional numeric matrix of input locations with shape \eqn{m
#'   \times d}. If \code{NULL} (default), it is set to \code{X}, resulting in a
#'   symmetric matrix.
#' @param theta A positive numeric value or vector specifying the kernel
#'   lengthscale(s). If a scalar, the same lengthscale is applied to all input
#'   dimensions. If a vector, it must be of length \code{d}, corresponding to
#'   anisotropic scaling.
#' @param kernel A character string specifying the kernel function. Must be one
#'   of \code{"gaussian"}, \code{"matern32"}, or \code{"matern52"}.
#' @param anisotropic Logical. If \code{TRUE} (default), \code{theta} is
#'   interpreted as a vector of per-dimension lengthscales. If \code{FALSE},
#'   \code{theta} is treated as a scalar.
#'
#' @return A numeric matrix of size \eqn{n \times m}, where each element
#'   \eqn{K_{ij}} gives the kernel similarity between input \eqn{X_i} and
#'   \eqn{X'_j}.
#'
#' @details Let \eqn{\mathbf{x}} and \eqn{\mathbf{x}'} denote two input points.
#'   The scaled distance is defined as
#'   \deqn{
#'      r = \left\| \frac{\mathbf{x} - \mathbf{x}'}{\boldsymbol{\theta}} \right\|_2.
#'   }
#'
#'   The available kernels are defined as:
#'   \itemize{
#'      \item \strong{Gaussian:}
#'      \deqn{
#'        k(\mathbf{x}, \mathbf{x}') = \exp(-r^2)
#'      }
#'      \item \strong{Matérn 5/2:}
#'      \deqn{
#'        k(\mathbf{x}, \mathbf{x}') = \left(1 + \sqrt{5} r + \frac{5}{3} r^2 \right) \exp(-\sqrt{5} r)
#'      }
#'      \item \strong{Matérn 3/2:}
#'      \deqn{
#'        k(\mathbf{x}, \mathbf{x}') = \left(1 + \sqrt{3} r \right) \exp(-\sqrt{3} r)
#'      }
#'   }
#'
#'   The function performs consistency checks on input dimensions and
#'   automatically broadcasts \code{theta} when it is a scalar.
#'
#' @references Rasmussen, C. E., & Williams, C. K. I. (2006). \emph{Gaussian
#'   Processes for Machine Learning}. MIT Press.
#'
#' @examples
#' # Basic usage with default Xprime = X
#' X <- matrix(runif(20), ncol = 2)
#' K1 <- kernel_matrix(X, theta = 0.2, kernel = "gaussian")
#'
#' # Anisotropic lengthscales with Matérn 5/2
#' K2 <- kernel_matrix(X, theta = c(0.1, 0.3), kernel = "matern52")
#'
#' # Isotropic Matérn 3/2
#' K3 <- kernel_matrix(X, theta = 1, kernel = "matern32", anisotropic = FALSE)
#'
#' # Use Xprime different from X
#' Xprime <- matrix(runif(10), ncol = 2)
#' K4 <- kernel_matrix(X, Xprime, theta = 0.2, kernel = "gaussian")
#'
#' @export

kernel_matrix <- function(X, Xprime = NULL, theta = 0.1,
                          kernel = c("gaussian", "matern52", "matern32"),
                          anisotropic = TRUE) {
  # Match the kernel argument explicitly
  kernel <- match.arg(kernel)

  if (is.null(Xprime)) Xprime <- X

  # Convert vectors to matrices
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  if (is.null(Xprime)) {
    Xprime <- X
    symmetric <- TRUE
  } else {
    if (is.vector(Xprime)) Xprime <- matrix(Xprime, ncol = 1)
    symmetric <- identical(X, Xprime)
  }

  # Check that input dimensions match
  if (ncol(X) != ncol(Xprime)) {
    stop("X and Xprime must have the same number of columns (i.e., input dimensions).")
  }

  # Expand or validate theta
  d <- ncol(X)
  if (anisotropic) {
    if (length(theta) == 1) theta <- rep(theta, d)
    if (length(theta) != d) stop("For anisotropic=TRUE, theta must be of length 1 or equal to number of columns.")
    # Rescale inputs by theta
    X_scaled <- sweep(X, 2, theta, "/")
    Xp_scaled <- sweep(Xprime, 2, theta, "/")
  } else {
    if (length(theta) != 1) stop("For anisotropic=FALSE, theta must be a scalar.")
    X_scaled <- X / theta
    Xp_scaled <- Xprime / theta
  }

  # Compute pairwise distances
  if (symmetric) {
    # Efficient computation when X == Xprime using symmetry
    dist <- as.matrix(dist(X_scaled))
    dist[dist < 0] <- 0  # numerical stability
    dist_sq <- dist^2

  } else {
    # General case: X and Xprime are different
    X_sq <- rowSums(X_scaled^2)
    Xp_sq <- rowSums(Xp_scaled^2)
    # Use identity: ||x - x'||^2 = ||x||^2 + ||x'||^2 - 2*x^T*x'
    dist_sq <- outer(X_sq, Xp_sq, "+") - 2 * tcrossprod(X_scaled, Xp_scaled)
    dist_sq[dist_sq < 0] <- 0  # numerical stability
    dist <- sqrt(dist_sq)
  }

  # Evaluate kernel function
  if (kernel == "gaussian") {
    K <- exp(-dist_sq)
  } else if (kernel == "matern52") {
    sqrt5 <- sqrt(5)
    K <- (1 + sqrt5 * dist + (5/3) * dist_sq) * exp(-sqrt5 * dist)
  } else if (kernel == "matern32") {
    sqrt3 <- sqrt(3)
    K <- (1 + sqrt3 * dist) * exp(-sqrt3 * dist)
  } else {
    stop(paste("Unsupported kernel type:", kernel))
  }

  return(K)
}
