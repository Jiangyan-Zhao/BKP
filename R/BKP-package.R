"_PACKAGE"

#' @name BKP-package
#'
#' @title Beta Kernel Process Modeling
#'
#' @description The \pkg{BKP} package provides tools for Bayesian
#'   nonparametric modeling of binary/binomial and categorical/multinomial
#'   response data using the Beta Kernel Process (BKP), the Dirichlet Kernel
#'   Process (DKP), and their scalable global-local approximations, TwinBKP and
#'   TwinDKP. These methods estimate latent probability surfaces by combining
#'   kernel-based smoothing with conjugate posterior updates.
#'
#'   The package supports model fitting, posterior prediction with uncertainty
#'   quantification, posterior simulation, and visualization in one- and
#'   two-dimensional input spaces. It also provides flexible prior
#'   specification, multiple kernel choices, hyperparameter tuning, and
#'   twinning-based global-local computation for scalable BKP and DKP modeling.
#'
#' @section Main Functions:
#' Core functionality is organized as follows:
#' \describe{
#'   \item{Model fitting}{
#'     Use \code{\link{fit_BKP}} and \code{\link{fit_DKP}} for full BKP and DKP
#'     models, and \code{\link{fit_TwinBKP}} and \code{\link{fit_TwinDKP}} for
#'     scalable global-local approximations using twinning-selected global
#'     subsets and local nearest-neighbour updates.
#'   }
#'   \item{Prediction}{
#'     Use \code{\link{predict}} with fitted \code{BKP}, \code{DKP},
#'     \code{TwinBKP}, or \code{TwinDKP} objects for posterior inference at new
#'     input locations, including posterior means, variances, and credible
#'     intervals. Decision labels, when needed, can be obtained from posterior
#'     means or method-specific outputs.
#'   }
#'   \item{Simulation}{
#'     Use \code{\link{simulate}} with fitted model objects to generate posterior
#'     draws of latent success probabilities or class-probability vectors.
#'   }
#'   \item{Visualization}{
#'     Use \code{\link{plot}} with fitted model objects to visualize predictions and
#'     associated uncertainty in one- and two-dimensional input spaces. For
#'     inputs with more than two dimensions, users can select one or two
#'     dimensions to display via the \code{dims} argument. TwinBKP and TwinDKP
#'     plots can optionally highlight the selected global subset.
#'   }
#'   \item{Summaries and printing}{
#'     Use \code{\link{summary}} and \code{\link{print}} to summarize or display fitted
#'     model objects and associated results.
#'   }
#'   \item{Extraction}{
#'     Use \code{\link{fitted}}, \code{\link{parameter}}, and \code{\link{quantile}} to
#'     extract fitted posterior means, model parameters, and posterior
#'     quantiles.
#'   }
#'   \item{Utilities}{
#'     Use \code{\link{kernel_matrix}} to construct kernel matrices,
#'     \code{\link{get_prior}} to construct BKP or DKP prior parameters, and
#'     \code{\link{loss_fun}} to evaluate leave-one-out tuning losses.
#'   }
#' }
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}.  arXiv. \doi{10.48550/arXiv.2508.10447}.
#'
#'   Vakayil A, Joseph VR (2024). \emph{A Global-Local Approximation Framework
#'   for Large-Scale Gaussian Process Modeling}. Technometrics, 66(2), 295--305.
#'   \doi{10.1080/00401706.2023.2296451}.
#'
#'   Rolland P, Kavis A, Immer A, Singla A, Cevher V (2019). \emph{Efficient learning of
#'   smooth probability functions from Bernoulli tests with guarantees}. In
#'   Proceedings of the 36th International Conference on Machine Learning, ICML
#'   2019, 9-15 June 2019, Long Beach, California, USA, volume 97 of Proceedings
#'   of Machine Learning Research, pp. 5459--5467. PMLR.
#'   \url{https://proceedings.mlr.press/v97/rolland19a.html}.
#'
#'   MacKenzie CA, Trafalis TB, Barker K (2014). \emph{A Bayesian Beta Kernel Model
#'   for Binary Classification and Online Learning Problems}. Statistical
#'   Analysis and Data Mining: The ASA Data Science Journal, 7(6), 434--449.
#'   \doi{10.1002/sam.11241}.
#'
#'   Goetschalckx R, Poupart P, Hoey J (2011). \emph{Continuous
#'   Correlated Beta Processes}. In Proceedings of the Twenty-Second
#'   International Joint Conference on Artificial Intelligence,
#'   IJCAI-11, p. 1269--1274. AAAI Press.
#'   \url{https://www.ijcai.org/Proceedings/11/Papers/215.pdf}.
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
