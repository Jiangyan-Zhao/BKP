"_PACKAGE"

#' @name BKP-package
#'
#' @title Beta Kernel Process Modeling
#'
#' @description The \pkg{BKP} package provides tools for Bayesian-inspired
#'   nonparametric modeling of binary/binomial and categorical/multinomial
#'   response data using the Beta Kernel Process (BKP), the Dirichlet Kernel
#'   Process (DKP), and scalable twinning-based approximations such as
#'   TwinBKP and TwinDKP. These methods estimate latent probability surfaces through
#'   localized kernel smoothing with conjugate posterior updates.
#'
#'   The package offers functionality for model fitting, posterior inference
#'   with uncertainty quantification, simulation of posterior draws, and
#'   visualization in one- and two-dimensional input spaces. It also supports
#'   flexible prior specification, kernel selection, hyperparameter tuning, and
#'   global-local approximation for scalable BKP and DKP modeling.
#'
#' @section Main Functions: Core functionality is organized as follows:
#' \describe{
#'   \item{\code{\link{fit_BKP}}, \code{\link{fit_DKP}}, \code{\link{fit_TwinBKP}}, \code{\link{fit_TwinDKP}}}{
#'     Fit BKP, DKP, TwinBKP, or TwinDKP models to binomial or multinomial response data.
#'     TwinBKP and TwinDKP provide scalable global-local approximations using
#'     twinning-selected global subsets and local nearest-neighbour updates.
#'   }
#'   \item{\code{\link{predict.BKP}}, \code{\link{predict.DKP}}, \code{\link{predict.TwinBKP}}, \code{\link{predict.TwinDKP}}}{
#'     Perform posterior inference at new input locations, including posterior
#'     means, variances, and credible intervals. Decision labels, when needed,
#'     can be obtained from posterior means or method-specific outputs.
#'   }
#'   \item{\code{\link{simulate.BKP}}, \code{\link{simulate.DKP}}, \code{\link{simulate.TwinBKP}}, \code{\link{simulate.TwinDKP}}}{
#'     Generate posterior draws of latent success probabilities or class
#'     probability vectors from fitted models.
#'   }
#'   \item{\code{\link{plot.BKP}}, \code{\link{plot.DKP}}, \code{\link{plot.TwinBKP}}, \code{\link{plot.TwinDKP}}}{
#'     Visualize model predictions and associated uncertainty in one- and
#'     two-dimensional input spaces. For inputs with more than two dimensions,
#'     users can select one or two dimensions to display via the \code{dims}
#'     argument. TwinBKP and TwinDKP plots can optionally highlight the selected global
#'     subset.
#'   }
#'   \item{\code{\link{summary.BKP}}, \code{\link{summary.DKP}}, \code{\link{summary.TwinBKP}}, \code{\link{summary.TwinDKP}},
#'         \code{\link{print.BKP}}, \code{\link{print.DKP}}, \code{\link{print.TwinBKP}}, \code{\link{print.TwinDKP}}}{
#'     Summarize or print fitted model objects and associated results.
#'   }
#'   \item{\code{\link{fitted.BKP}}, \code{\link{fitted.DKP}}, \code{\link{fitted.TwinBKP}}, \code{\link{fitted.TwinDKP}},
#'         \code{\link{parameter}}, \code{\link{quantile.BKP}},
#'         \code{\link{quantile.DKP}}, \code{\link{quantile.TwinBKP}}, \code{\link{quantile.TwinDKP}}}{
#'     Extract fitted posterior means, model parameters, and posterior quantiles.
#'   }
#' }
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}.  arXiv. \doi{10.48550/arXiv.2508.10447}
#'
#'   Vakayil A, Joseph VR (2024). \emph{A Global-Local Approximation Framework
#'   for Large-Scale Gaussian Process Modeling}. Technometrics, 66(2), 295--305.
#'
#'   Rolland P, Kavis A, Singla A, Cevher V (2019). \emph{Efficient learning of
#'   smooth probability functions from Bernoulli tests with guarantees}. In
#'   Proceedings of the 36th International Conference on Machine Learning, ICML
#'   2019, 9-15 June 2019, Long Beach, California, USA, volume 97 of Proceedings
#'   of Machine Learning Research, pp. 5459-5467. PMLR.
#'
#'   MacKenzie CA, Trafalis TB, Barker K (2014). \emph{A Bayesian Beta Kernel Model
#'   for Binary Classification and Online Learning Problems}. Statistical
#'   Analysis and Data Mining: The ASA Data Science Journal, 7(6), 434-449.
#'
#'   Goetschalckx R, Poupart P, Hoey J (2011). \emph{Continuous
#'   Correlated Beta Processes}. In Proceedings of the Twenty-Second
#'   International Joint Conference on Artificial Intelligence - Volume Volume
#'   Two, IJCAI’11, p. 1269-1274. AAAI Press.
#'
#'
#' @importFrom dirmult rdirichlet
#' @importFrom graphics abline legend lines par points polygon text
#' @importFrom grDevices hcl.colors rainbow
#' @importFrom grid gpar textGrob unit
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point geom_raster
#' @importFrom ggplot2 geom_contour scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 scale_fill_viridis_c scale_color_manual scale_fill_manual
#' @importFrom ggplot2 theme_bw theme element_blank element_rect element_text
#' @importFrom ggplot2 guide_legend guide_colorbar coord_cartesian annotate geom_hline
#' @importFrom ggplot2 labs scale_color_discrete guides scale_fill_brewer scale_shape_manual
#' @importFrom lattice levelplot panel.contourplot panel.levelplot panel.points
#' @importFrom rlang .data
#' @importFrom stats as.formula median qbeta quantile rbeta rgamma sd simulate
#' @importFrom tgp lhs
#' @importFrom utils head
#' @import nloptr
#' @importFrom Rcpp evalCpp
#' @useDynLib BKP, .registration = TRUE
NULL
