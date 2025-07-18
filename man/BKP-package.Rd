% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BKP-package.R
\name{BKP-package}
\alias{BKP-package}
\title{Beta and Dirichlet Kernel Process Modeling}
\description{
The \pkg{BKP} package provides tools for nonparametric modeling
of binary, binomial, or multinomial response data using the Beta Kernel
Process (BKP) and its extension, the Dirichlet Kernel Process (DKP). These
methods estimate latent probability surfaces through localized kernel
smoothing under a Bayesian framework.

The package includes functionality for model fitting, probabilistic
prediction with uncertainty quantification, posterior simulation, and
visualization in both one- and two-dimensional input spaces. It also
supports hyperparameter tuning and flexible prior specification.
}
\section{Main Functions}{
 Core functionality is organized into the following
groups:
\describe{
\item{\code{\link{fit.BKP}}, \code{\link{fit.DKP}}}{
Fit a BKP or DKP model to (multi)binomial response data.
}
\item{\code{\link{predict.BKP}}, \code{\link{predict.DKP}}}{
Perform posterior predictive inference at new input locations, including
predictive means, variances, and credible intervals.
Classification labels are returned automatically
when observations represent single trials (i.e., binary outcomes).
}
\item{\code{\link{simulate.BKP}}, \code{\link{simulate.DKP}}}{
Draw simulated responses from the posterior predictive distribution of a fitted model.
}
\item{\code{\link{plot.BKP}}, \code{\link{plot.DKP}}}{
Visualize model predictions and uncertainty bands in 1D and 2D input spaces.
}
\item{\code{\link{summary.BKP}}, \code{\link{summary.DKP}}, \code{\link{print.BKP}}, \code{\link{print.DKP}}}{
Summarize or print details of a fitted BKP or DKP model.
}
}
}

