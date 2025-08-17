"_PACKAGE"

#' @name BKP-package
#'
#' @title Beta Kernel Process Modeling
#'
#' @description The \pkg{BKP} package provides tools for nonparametric modeling
#'   of binary/binomial or categorical/multinomial response data using the Beta
#'   Kernel Process (BKP) and its extension, the Dirichlet Kernel Process (DKP).
#'   These methods estimate latent probability surfaces through localized kernel
#'   smoothing under a Bayesian framework.
#'
#'   The package includes functionality for model fitting, probabilistic
#'   prediction with uncertainty quantification, posterior simulation, and
#'   visualization in both one- and two-dimensional input spaces. It also
#'   supports hyperparameter tuning and flexible prior specification.
#'
#' @section Main Functions: Core functionality is organized into the following
#'   groups:
#' \describe{
#'   \item{\code{\link{fit.BKP}}, \code{\link{fit.DKP}}}{
#'     Fit a BKP or DKP model to (multi)binomial response data.
#'   }
#'   \item{\code{\link{predict.BKP}}, \code{\link{predict.DKP}}}{
#'     Perform posterior predictive inference at new input locations, including
#'     predictive means, variances, and credible intervals.
#'     Classification labels are returned automatically
#'     when observations represent single trials (i.e., binary outcomes).
#'   }
#'   \item{\code{\link{simulate.BKP}}, \code{\link{simulate.DKP}}}{
#'     Draw simulated responses from the posterior predictive distribution of a fitted model.
#'   }
#'   \item{\code{\link{plot.BKP}}, \code{\link{plot.DKP}}}{
#'     Visualize model predictions and uncertainty bands in 1D and 2D input spaces.
#'   }
#'   \item{\code{\link{summary.BKP}}, \code{\link{summary.DKP}}, \code{\link{print.BKP}}, \code{\link{print.DKP}}}{
#'     Summarize or print details of a fitted BKP or DKP model.
#'   }
#' }
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}.  arXiv.
#'   https://doi.org/10.48550/arXiv.2508.10447.
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
#' @importFrom graphics abline legend lines par points polygon text
#' @importFrom grDevices hcl.colors rainbow
#' @importFrom grid gpar textGrob
#' @importFrom gridExtra grid.arrange
#' @importFrom lattice levelplot panel.contourplot panel.levelplot panel.points
#' @importFrom optimx multistart
#' @importFrom stats as.formula qbeta rbeta rgamma
#' @importFrom tgp lhs
NULL
