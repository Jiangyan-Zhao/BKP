#' @title Compute Kernel Matrix Between Input Locations
#'
#' @description Computes the kernel matrix between two sets of input locations
#'   using a specified kernel function. Supports both isotropic and anisotropic
#'   lengthscales. Available kernels include the Gaussian, Matérn 5/2,
#'   Matérn 3/2, and Wendland compactly supported kernel.
#'
#' @param X A numeric matrix (or vector) of input locations with shape \eqn{n
#'   \times d}.
#' @param Xprime An optional numeric matrix of input locations with shape \eqn{m
#'   \times d}. If \code{NULL} (default), it is set to \code{X}, resulting in a
#'   symmetric matrix.
#' @param theta A positive numeric value or vector specifying the kernel
#'   lengthscale(s). If \code{isotropic = TRUE} (default), this must be a
#'   scalar shared by all input dimensions. If \code{isotropic = FALSE}, this
#'   can be a scalar (broadcasted) or a vector of length \code{d} for
#'   per-dimension scaling.
#' @param kernel A character string specifying the kernel function. Must be one
#'   of \code{"gaussian"}, \code{"matern32"}, \code{"matern52"}, or
#'   \code{"wendland"}.
#' @param isotropic Logical. If \code{TRUE} (default), use a single shared
#'   lengthscale across dimensions. If \code{FALSE}, use per-dimension
#'   lengthscales.
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
#'      \item \strong{Wendland:}
#'      \deqn{
#'        k(\mathbf{x}, \mathbf{x}') = (q r + 1)\max(0,1-r)^q,
#'        \quad q = \lfloor d/2 \rfloor + 3
#'      }
#'   }
#'
#'   The function performs consistency checks on input dimensions and
#'   automatically broadcasts \code{theta} when it is a scalar.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}.  arXiv.
#'   https://doi.org/10.48550/arXiv.2508.10447.
#'
#'   Rasmussen, C. E., & Williams, C. K. I. (2006). \emph{Gaussian
#'   Processes for Machine Learning}. MIT Press.
#'
#' @examples
#' # Basic usage with default Xprime = X
#' X <- matrix(runif(20), ncol = 2)
#' K1 <- kernel_matrix(X, theta = 0.2, kernel = "gaussian")
#'
#' # Anisotropic lengthscales with Matérn 5/2
#' K2 <- kernel_matrix(X, theta = c(0.1, 0.3), kernel = "matern52", isotropic = FALSE)
#'
#' # Isotropic Matérn 3/2
#' K3 <- kernel_matrix(X, theta = 1, kernel = "matern32")
#'
#' # Use Xprime different from X
#' Xprime <- matrix(runif(10), ncol = 2)
#' K4 <- kernel_matrix(X, Xprime, theta = 0.2, kernel = "gaussian")
#'
#' @export

kernel_matrix <- function(X, Xprime = NULL, theta = 0.1,
                          kernel = c("gaussian", "matern52", "matern32", "wendland"),
                          isotropic = TRUE) {
  # ---- Argument checking ----
  if (!is.numeric(X)) stop("'X' must be numeric or a numeric matrix.")
  if (anyNA(X)) stop("'X' contains NA values.")

  if (!is.null(Xprime)) {
    if (!is.numeric(Xprime)) stop("'Xprime' must be numeric or a numeric matrix.")
    if (anyNA(Xprime)) stop("'Xprime' contains NA values.")
  }

  if (!is.numeric(theta) || length(theta) < 1 || anyNA(theta) || any(theta <= 0)) {
    stop("'theta' must be numeric and strictly positive.")
  }
  if (!is.logical(isotropic) || length(isotropic) != 1) {
    stop("'isotropic' must be a single logical value.")
  }

  kernel <- match.arg(kernel)

  # Convert vector -> matrix (n x 1)
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  if (!is.null(Xprime) && is.vector(Xprime)) Xprime <- matrix(Xprime, ncol = 1)

  # Dimension checks
  if (!is.null(Xprime) && ncol(X) != ncol(Xprime)) {
    stop("'X' and 'Xprime' must have the same number of columns (input dimensions).")
  }

  d <- ncol(X)

  # theta checks aligned with doc
  if (isotropic) {
    if (length(theta) != 1) {
      stop("For isotropic=TRUE, 'theta' must be a scalar.")
    }
  } else {
    if (length(theta) == 1) {
      theta <- rep(theta, d)
    } else if (length(theta) != d) {
      stop("For isotropic=FALSE, 'theta' must be scalar or of length equal to ncol(X).")
    }
  }

  kernel_matrix_rcpp(
    X = X,
    Xprime = Xprime,
    theta = as.numeric(theta),
    kernel = kernel,
    isotropic = isTRUE(isotropic)
  )
}
