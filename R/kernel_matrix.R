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
#'   of \code{"gaussian"}, \code{"matern52"}, \code{"matern32"}, or
#'   \code{"wendland"}.
#' @param isotropic Logical. If \code{TRUE} (default), use a single shared
#'   lengthscale across dimensions. If \code{FALSE}, use per-dimension
#'   lengthscales.
#'
#' @return A numeric matrix of size \eqn{n \times m}, where each element
#'   \eqn{K_{ij}} gives the kernel similarity between input \eqn{X_i} and
#'   \eqn{X'_j}.
#'
#' @details Let \eqn{\mathbf{x}} and \eqn{\mathbf{x}^{\prime}} denote two input
#'   points. The scaled Euclidean distance is
#'   \deqn{
#'     h(\mathbf{x}, \mathbf{x}^{\prime}; \boldsymbol{\theta}) =
#'     \left\|
#'       \frac{\mathbf{x} - \mathbf{x}^{\prime}}{\boldsymbol{\theta}}
#'     \right\|_2.
#'   }
#'   For isotropic kernels, \eqn{\boldsymbol{\theta}} is a scalar shared by all
#'   input dimensions. For anisotropic kernels, \eqn{\boldsymbol{\theta}} is a
#'   vector of dimension-specific lengthscales.
#'
#'   The available kernels are
#'   \itemize{
#'     \item \strong{Gaussian:}
#'     \deqn{
#'       k(\mathbf{x}, \mathbf{x}^{\prime}) = \exp(-h^2).
#'     }
#'     \item \strong{Matérn 5/2:}
#'     \deqn{
#'       k(\mathbf{x}, \mathbf{x}^{\prime}) =
#'       \left(1 + \sqrt{5}h + \frac{5}{3}h^2\right)
#'       \exp(-\sqrt{5}h).
#'     }
#'     \item \strong{Matérn 3/2:}
#'     \deqn{
#'       k(\mathbf{x}, \mathbf{x}^{\prime}) =
#'       \left(1 + \sqrt{3}h\right)\exp(-\sqrt{3}h).
#'     }
#'     \item \strong{Wendland:}
#'     \deqn{
#'       k(\mathbf{x}, \mathbf{x}^{\prime}) =
#'       (\zeta h + 1)\max(0, 1 - h)^{\zeta},
#'       \qquad \zeta = \lfloor d/2 \rfloor + 3.
#'     }
#'   }
#'
#'   The returned matrix has one row for each row of \code{X} and one column for
#'   each row of \code{Xprime}. If \code{Xprime = NULL}, the function returns the
#'   symmetric kernel matrix for \code{X}.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}. arXiv. <doi:10.48550/arXiv.2508.10447>.
#'
#'   Rasmussen, C. E., & Williams, C. K. I. (2006). \emph{Gaussian
#'   Processes for Machine Learning}. MIT Press.
#'
#'   Wendland, H. (1995). Piecewise polynomial, positive definite and compactly
#'   supported radial functions of minimal degree. \emph{Advances in
#'   Computational Mathematics}, 4(1), 389--396.
#'
#' @examples
#' set.seed(123)
#'
#' X <- matrix(runif(20), ncol = 2)
#'
#' # Compute kernel matrices for all implemented kernels
#' kernels <- c("gaussian", "matern52", "matern32", "wendland")
#' K_list <- lapply(kernels, function(k) {
#'   kernel_matrix(X, theta = 0.2, kernel = k)
#' })
#' names(K_list) <- kernels
#'
#' # Inspect part of the Gaussian and Wendland kernel matrices
#' K_list$gaussian[1:3, 1:3]
#' K_list$wendland[1:3, 1:3]
#'
#' # Anisotropic lengthscales
#' K_aniso <- kernel_matrix(
#'   X,
#'   theta = c(0.1, 0.3),
#'   kernel = "matern52",
#'   isotropic = FALSE
#' )
#'
#' # Cross-kernel matrix between two input sets
#' Xprime <- matrix(runif(10), ncol = 2)
#' K_cross <- kernel_matrix(X, Xprime, theta = 0.2, kernel = "gaussian")
#' dim(K_cross)
#'
#' @export

kernel_matrix <- function(X, Xprime = NULL, theta = 0.1,
                          kernel = c("gaussian", "matern52", "matern32", "wendland"),
                          isotropic = TRUE) {
  # ---- Argument checking ----
  if (!is.numeric(X)) stop("'X' must be numeric or a numeric matrix.")
  if (any(!is.finite(X))) stop("'X' must contain only finite numeric values.")

  if (!is.null(Xprime)) {
    if (!is.numeric(Xprime)) stop("'Xprime' must be numeric or a numeric matrix.")
    if (any(!is.finite(Xprime))) stop("'Xprime' must contain only finite numeric values.")
  }

  if (!is.numeric(theta) || length(theta) < 1 || any(!is.finite(theta)) || any(theta <= 0)) {
    stop("'theta' must contain only finite numeric and strictly positive values.")
  }
  if (!is.logical(isotropic) || length(isotropic) != 1) {
    stop("'isotropic' must be a single logical value.")
  }

  kernel <- match.arg(kernel)

  # Convert vector -> matrix (n x 1)
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  if (!is.null(Xprime) && is.vector(Xprime)) Xprime <- matrix(Xprime, ncol = 1)

  # Empty input checks
  if (nrow(X) < 1L || ncol(X) < 1L) {
    stop("'X' must have at least one row and one column.")
  }

  if (!is.null(Xprime) && (nrow(Xprime) < 1L || ncol(Xprime) < 1L)) {
    stop("'Xprime' must have at least one row and one column.")
  }

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
