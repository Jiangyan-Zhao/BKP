// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;


enum KernelType {
  GAUSSIAN = 0,
  MATERN52 = 1,
  MATERN32 = 2,
  WENDLAND = 3
};

static inline KernelType parse_kernel_type(const std::string& kernel) {
  if (kernel == "gaussian") {
    return GAUSSIAN;
  }
  if (kernel == "matern52") {
    return MATERN52;
  }
  if (kernel == "matern32") {
    return MATERN32;
  }
  return WENDLAND;
}

static inline double scaled_dist_sq_scalar(
    const arma::mat& Xm,
    const arma::mat& Xpm,
    const double inv_th,
    const arma::uword i,
    const arma::uword j
) {
  const arma::uword d = Xm.n_cols;
  double out = 0.0;

  for (arma::uword k = 0; k < d; ++k) {
    const double diff = (Xm(i, k) - Xpm(j, k)) * inv_th;
    out += diff * diff;
  }

  return out;
}

static inline double scaled_dist_sq_aniso(
    const arma::mat& Xm,
    const arma::mat& Xpm,
    const arma::vec& theta,
    const arma::uword i,
    const arma::uword j
) {
  const arma::uword d = Xm.n_cols;
  double out = 0.0;

  for (arma::uword k = 0; k < d; ++k) {
    const double diff = (Xm(i, k) - Xpm(j, k)) / theta[k];
    out += diff * diff;
  }

  return out;
}

static inline double kernel_from_dist_sq_fast(
    const double dist_sq,
    const KernelType kernel_type,
    const double sqrt3,
    const double sqrt5,
    const double q_w
) {
  switch (kernel_type) {
  case GAUSSIAN:
    return std::exp(-dist_sq);

  case MATERN52: {
    const double dist = std::sqrt(dist_sq);
    return (1.0 + sqrt5 * dist + (5.0 / 3.0) * dist_sq) *
      std::exp(-sqrt5 * dist);
  }

  case MATERN32: {
    const double dist = std::sqrt(dist_sq);
    return (1.0 + sqrt3 * dist) * std::exp(-sqrt3 * dist);
  }

  case WENDLAND: {
    const double dist = std::sqrt(dist_sq);
    const double one_minus = std::max(0.0, 1.0 - dist);
    return (q_w * dist + 1.0) * std::pow(one_minus, q_w);
  }
  }

  return R_NaReal;
}

arma::mat kernel_matrix_arma_loop(
    const arma::mat& Xm,
    const arma::mat& Xpm,
    const arma::vec& theta,
    const std::string& kernel,
    const bool isotropic,
    const bool symmetric
) {
  const arma::uword n = Xm.n_rows;
  const arma::uword m = symmetric ? Xm.n_rows : Xpm.n_rows;
  const arma::uword d = Xm.n_cols;
  const KernelType kernel_type = parse_kernel_type(kernel);
  const double sqrt3 = std::sqrt(3.0);
  const double sqrt5 = std::sqrt(5.0);
  const double q_w = std::floor(static_cast<double>(d) / 2.0) + 3.0;
  const bool scalar_lengthscale = isotropic || theta.n_elem == 1;
  arma::mat K(n, m, arma::fill::none);

  if (scalar_lengthscale) {
    const double inv_th = 1.0 / theta[0];

    if (symmetric) {
      for (arma::uword i = 0; i < n; ++i) {
        K(i, i) = 1.0;
        for (arma::uword j = i + 1; j < n; ++j) {
          const double dist_sq = scaled_dist_sq_scalar(Xm, Xm, inv_th, i, j);
          const double val = kernel_from_dist_sq_fast(
            dist_sq, kernel_type, sqrt3, sqrt5, q_w
          );
          K(i, j) = val;
          K(j, i) = val;
        }
      }
    } else {
      for (arma::uword i = 0; i < n; ++i) {
        for (arma::uword j = 0; j < m; ++j) {
          const double dist_sq = scaled_dist_sq_scalar(Xm, Xpm, inv_th, i, j);
          K(i, j) = kernel_from_dist_sq_fast(
            dist_sq, kernel_type, sqrt3, sqrt5, q_w
          );
        }
      }
    }
  } else {
    if (symmetric) {
      for (arma::uword i = 0; i < n; ++i) {
        K(i, i) = 1.0;
        for (arma::uword j = i + 1; j < n; ++j) {
          const double dist_sq = scaled_dist_sq_aniso(Xm, Xm, theta, i, j);
          const double val = kernel_from_dist_sq_fast(
            dist_sq, kernel_type, sqrt3, sqrt5, q_w
          );
          K(i, j) = val;
          K(j, i) = val;
        }
      }
    } else {
      for (arma::uword i = 0; i < n; ++i) {
        for (arma::uword j = 0; j < m; ++j) {
          const double dist_sq = scaled_dist_sq_aniso(Xm, Xpm, theta, i, j);
          K(i, j) = kernel_from_dist_sq_fast(
            dist_sq, kernel_type, sqrt3, sqrt5, q_w
          );
        }
      }
    }
  }

  return K;
}

// -----------------------------------------------------------------------------
// Internal C++ kernel matrix engine.
//
// This function assumes that user-facing validation has already been handled
// by the R wrapper kernel_matrix(). It intentionally avoids scanning Xm/Xpm
// for NA/Inf, because this function may be called repeatedly during kernel
// hyperparameter optimization.
// -----------------------------------------------------------------------------

arma::mat kernel_matrix_arma(
    const arma::mat& Xm,
    const arma::mat& Xpm,
    const arma::vec& theta,
    const std::string& kernel,
    const bool isotropic,
    const bool symmetric
) {
  const arma::uword d = Xm.n_cols;
  const double n_pairs = static_cast<double>(Xm.n_rows) *
    static_cast<double>(symmetric ? Xm.n_rows : Xpm.n_rows);

  // GEMM is usually faster for moderate kernel matrices.  A 2000 x 2000
  // double matrix is about 32 MB; even with several temporaries this is
  // acceptable on ordinary machines.  For larger kernels, switch to the
  // loop engine to cap peak memory.
  const double loop_cutoff_pairs = 4e6;

  if (n_pairs > loop_cutoff_pairs) {
    return kernel_matrix_arma_loop(Xm, Xpm, theta, kernel, isotropic, symmetric);
  }

  // ---- Scale inputs by lengthscale(s) ----
  arma::mat X_scaled;
  arma::mat Xp_scaled;

  if (isotropic) {
    const double inv_th = 1.0 / theta[0];

    X_scaled  = Xm  * inv_th;
    Xp_scaled = Xpm * inv_th;

  } else {
    arma::rowvec th(d);

    if (theta.n_elem == 1) {
      th.fill(theta[0]);
    } else {
      th = theta.t();
    }

    X_scaled  = Xm.each_row()  / th;
    Xp_scaled = Xpm.each_row() / th;
  }

  // ---- Pairwise squared distances ----
  arma::mat dist_sq;

  if (symmetric) {
    arma::vec g = arma::sum(arma::square(X_scaled), 1); // n x 1
    arma::mat G = X_scaled * X_scaled.t();              // n x n

    dist_sq = arma::repmat(g, 1, g.n_elem) + arma::repmat(g.t(), g.n_elem, 1) - 2.0 * G;

  } else {
    arma::vec g  = arma::sum(arma::square(X_scaled), 1);  // n x 1
    arma::vec gp = arma::sum(arma::square(Xp_scaled), 1); // m x 1
    arma::mat G  = X_scaled * Xp_scaled.t();              // n x m

    dist_sq = arma::repmat(g, 1, gp.n_elem) + arma::repmat(gp.t(), g.n_elem, 1) - 2.0 * G;
  }

  // Guard against tiny negative values caused by floating-point roundoff.
  dist_sq.transform([](double v) {
    return (v < 0.0) ? 0.0 : v;
  });

  // ---- Kernel evaluation ----
  arma::mat K;

  if (kernel == "gaussian") {
    K = arma::exp(-dist_sq);

  } else if (kernel == "matern52") {
    const double sqrt5 = std::sqrt(5.0);
    arma::mat dist = arma::sqrt(dist_sq);

    K = (1.0 + sqrt5 * dist + (5.0 / 3.0) * dist_sq) % arma::exp(-sqrt5 * dist);

  } else if (kernel == "matern32") {
    const double sqrt3 = std::sqrt(3.0);
    arma::mat dist = arma::sqrt(dist_sq);

    K = (1.0 + sqrt3 * dist) % arma::exp(-sqrt3 * dist);

  } else {
    // Wendland compactly supported kernel:
    // K(r) = (q * r + 1) * max(0, 1 - r)^q,
    // q = floor(d / 2) + 3.
    const double q_w = std::floor(static_cast<double>(d) / 2.0) + 3.0;

    arma::mat dist = arma::sqrt(dist_sq);
    arma::mat one_minus = 1.0 - dist;

    one_minus.transform([](double v) {
      return (v > 0.0) ? v : 0.0;
    });

    K = (q_w * dist + 1.0) % arma::pow(one_minus, q_w);
  }

  return K;
}

// -----------------------------------------------------------------------------
// R-facing wrapper.
//
// The R wrapper kernel_matrix() is responsible for:
// - checking numeric inputs;
// - checking finite values;
// - checking theta;
// - matching kernel names;
// - converting vector inputs to n x 1 matrices;
// - checking input dimensions.
//
// Therefore this C++ wrapper only converts R objects to Armadillo views and
// delegates computation to kernel_matrix_arma().
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
arma::mat kernel_matrix_rcpp(
    SEXP X,
    SEXP Xprime = R_NilValue,
    NumericVector theta = NumericVector::create(0.1),
    std::string kernel = "gaussian",
    bool isotropic = true
) {
  NumericMatrix X_r(X);

  arma::mat Xm(
      X_r.begin(),
      X_r.nrow(),
      X_r.ncol(),
      false
  );

  arma::vec theta_arma(
      theta.begin(),
      static_cast<arma::uword>(theta.size()),
      false
  );

  if (Rf_isNull(Xprime)) {
    return kernel_matrix_arma(
      Xm,
      Xm,
      theta_arma,
      kernel,
      isotropic,
      true
    );
  }

  NumericMatrix Xp_r(Xprime);

  arma::mat Xpm(
      Xp_r.begin(),
      Xp_r.nrow(),
      Xp_r.ncol(),
      false
  );

  return kernel_matrix_arma(
    Xm,
    Xpm,
    theta_arma,
    kernel,
    isotropic,
    false
  );
}
