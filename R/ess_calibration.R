# Internal helpers for optional BKP effective-sample-size calibration.

.bkp_check_unique_locations <- function(Xnorm) {
  if (anyDuplicated(as.data.frame(Xnorm)) > 0L) {
    stop(
      paste0(
        "ESS calibration with ess = 'shepard' requires unique input locations; ",
        "duplicated rows in 'X' are not supported because strict Shepard ",
        "interpolation is not well-defined for duplicates."
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.bkp_shepard_m <- function(Xquery_norm, Xtrain_norm, m, power = 2) {
  Xquery_norm <- as.matrix(Xquery_norm)
  Xtrain_norm <- as.matrix(Xtrain_norm)
  m <- as.numeric(m)

  if (ncol(Xquery_norm) != ncol(Xtrain_norm)) {
    stop("'Xquery_norm' and 'Xtrain_norm' must have the same number of columns.", call. = FALSE)
  }

  out <- numeric(nrow(Xquery_norm))
  for (i in seq_len(nrow(Xquery_norm))) {
    diff <- sweep(Xtrain_norm, 2L, Xquery_norm[i, ], "-")
    dist_sq <- rowSums(diff^2)
    exact <- which(dist_sq == 0)

    if (length(exact) > 0L) {
      out[i] <- m[exact[1L]]
    } else {
      weights <- dist_sq^(-power / 2)
      out[i] <- sum(weights * m) / sum(weights)
    }
  }
  out
}

.bkp_shepard_m_loo <- function(Xnorm, m, power = 2) {
  Xnorm <- as.matrix(Xnorm)
  m <- as.numeric(m)

  .bkp_check_unique_locations(Xnorm)

  n <- nrow(Xnorm)
  if (n < 2L) {
    stop(
      "ESS calibration with ess = 'shepard' requires at least two input locations for leave-one-out calibration.",
      call. = FALSE
    )
  }

  out <- numeric(n)
  for (i in seq_len(n)) {
    diff <- sweep(Xnorm[-i, , drop = FALSE], 2L, Xnorm[i, ], "-")
    dist_sq <- rowSums(diff^2)
    weights <- dist_sq^(-power / 2)
    out[i] <- sum(weights * m[-i]) / sum(weights)
  }
  out
}

.bkp_ess_calibration <- function(Xquery_norm, Xtrain_norm, m, K) {
  .bkp_check_unique_locations(Xtrain_norm)

  m_kernel <- as.vector(K %*% as.numeric(m))
  m_shepard <- .bkp_shepard_m(Xquery_norm, Xtrain_norm, m, power = 2)
  rho <- apply(K, 1L, max)
  m_target <- rho * m_shepard

  scale <- numeric(length(m_kernel))
  positive_kernel_mass <- m_kernel > 0
  scale[positive_kernel_mass] <- m_target[positive_kernel_mass] / m_kernel[positive_kernel_mass]

  list(
    scale = scale,
    m_kernel = m_kernel,
    m_shepard = m_shepard,
    m_target = m_target,
    rho = rho
  )
}

.bkp_ess_none_info <- function(K, m) {
  list(
    scale = rep(1, nrow(K)),
    m_kernel = as.vector(K %*% as.numeric(m)),
    m_shepard = NULL,
    m_target = NULL,
    rho = apply(K, 1L, max)
  )
}
