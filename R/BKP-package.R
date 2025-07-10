"_PACKAGE"

#' Beta Kernel Process Modeling
#'
#' @name BKP-package
#'
#' @description \pkg{BKP} provides tools for modeling binary or binomial
#'   response data using the Beta Kernel Process (BKP), a flexible nonparametric
#'   method that estimates latent probability surfaces by local kernel smoothing
#'   with beta-binomial likelihood. The package supports efficient model
#'   fitting, prediction with uncertainty quantification, and visualization
#'   tools for both 1D and 2D covariates.
#'
#' @section Functions:
#'
#'   Key functions included in the package:
#' \describe{
#'   \item{\link{fit.BKP}}{
#'   Fit a BKP model to binomial or binary response data.
#'   }
#'   \item{\link{predict.BKP}}{
#'   Predict response probabilities and compute credible/confidence intervals at new input locations.
#'   }
#'   \item{\link{plot.BKP}}{
#'   Visualize fitted BKP models for 1D or 2D inputs.
#'   }
#'   \item{\link{print.BKP}, \link{summary.BKP}}{
#'   Print or summarize the BKP model fit results.
#'   }
#' }
#'
#' @importFrom lattice levelplot panel.levelplot panel.contourplot
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices hcl.colors
#' @importFrom graphics legend lines points polygon
#' @importFrom stats as.formula qbeta
#' @importFrom tgp lhs
#' @importFrom optimx multistart
NULL
