"_PACKAGE"

#' @title Beta Kernel Process Modeling
#'
#' @name BKP-package
#'
#' @description
#' \pkg{BKP} provides tools for modeling binary or binomial
#' response data using the Beta Kernel Process (BKP), a flexible nonparametric
#' method that estimates latent probability surfaces via local kernel smoothing
#' under a beta-binomial framework. The package supports efficient model
#' fitting, probabilistic prediction with uncertainty quantification, and
#' visualization tools for 1D and 2D input spaces. It also provides simulation
#' utilities for generating synthetic data from a fitted model.
#'
#' @section Functions:
#' Key functions included in the package:
#' \describe{
#'   \item{\link{fit.BKP}, \link{fit.DKP}}{
#'     Fit a BKP model to binomial or binary response data.
#'   }
#'   \item{\link{predict.BKP}, \link{predict.DKP}}{
#'     Predict success probabilities and confidence/credible intervals at new input locations.
#'     Classification labels are automatically returned when \code{m = 1} for all training data.
#'   }
#'   \item{\link{simulate.BKP}, \link{simulate.DKP}}{
#'     Generate simulated binary or binomial responses from a fitted BKP model.
#'   }
#'   \item{\link{plot.BKP}, \link{plot.DKP}}{
#'     Visualize fitted BKP models for 1D or 2D inputs.
#'   }
#'   \item{\link{print.BKP}, \link{print.DKP}, \link{summary.BKP}, \link{summary.DKP}}{
#'     Print or summarize model fit results.
#'   }
#' }
#'
#' @importFrom lattice levelplot panel.levelplot panel.contourplot
#' @importFrom grid gpar textGrob
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices hcl.colors
#' @importFrom graphics legend lines points polygon
#' @importFrom stats as.formula qbeta rbeta rgamma
#' @importFrom tgp lhs
#' @importFrom optimx multistart
NULL
