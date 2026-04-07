// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// ---------- helpers ----------
static inline void check_no_na_numeric_vector(const NumericVector& v, const std::string& name) {
  for (R_xlen_t i = 0; i < v.size(); ++i) {
    if (NumericVector::is_na(v[i])) {
      stop("'" + name + "' contains NA values.");
    }
  }
}

static inline arma::mat as_matrix_1d_or_2d(SEXP x, const std::string& name) {
  if (!Rf_isNumeric(x)) {
    stop("'" + name + "' must be numeric or a numeric matrix.");
  }

  SEXP dim = Rf_getAttrib(x, R_DimSymbol);

  // numeric vector -> n x 1 matrix
  if (dim == R_NilValue) {
    NumericVector v(x);
    check_no_na_numeric_vector(v, name);

    arma::mat out(v.size(), 1);
    for (R_xlen_t i = 0; i < v.size(); ++i) out(i, 0) = v[i];
    return out;
  }

  // numeric matrix
  NumericMatrix m(x);
  NumericVector flat = as<NumericVector>(x);
  check_no_na_numeric_vector(flat, name);

  arma::mat out(m.begin(), m.nrow(), m.ncol(), false); // view
  return arma::mat(out); // deep copy
}

// [[Rcpp::export]]
arma::mat kernel_matrix_rcpp(
    SEXP X,
    SEXP Xprime = R_NilValue,
    NumericVector theta = NumericVector::create(0.1),
    std::string kernel = "gaussian",
    bool isotropic = true
) {
  // ---- argument checks ----
  if (theta.size() < 1) {
    stop("'theta' must be numeric and strictly positive.");
  }
  for (R_xlen_t i = 0; i < theta.size(); ++i) {
    if (NumericVector::is_na(theta[i]) || theta[i] <= 0.0) {
      stop("'theta' must be numeric and strictly positive.");
    }
  }

  if (!(kernel == "gaussian" || kernel == "matern52" || kernel == "matern32" || kernel == "wendland")) {
    stop("'kernel' must be one of: 'gaussian', 'matern52', 'matern32', 'wendland'.");
  }

  arma::mat Xm = as_matrix_1d_or_2d(X, "X");
  arma::mat Xpm;
  bool symmetric = false;

  if (Rf_isNull(Xprime)) {
    Xpm = Xm;
    symmetric = true;
  } else {
    Xpm = as_matrix_1d_or_2d(Xprime, "Xprime");
    symmetric = false;
  }

  if (Xm.n_cols != Xpm.n_cols) {
    stop("'X' and 'Xprime' must have the same number of columns (input dimensions).");
  }

  const arma::uword d = Xm.n_cols;
  arma::mat X_scaled, Xp_scaled;

  if (isotropic) {
    if (theta.size() != 1) {
      stop("For isotropic=TRUE, 'theta' must be a scalar.");
    }
    const double th = theta[0];
    X_scaled  = Xm / th;
    Xp_scaled = Xpm / th;
  } else {
    NumericVector theta_use = clone(theta);
    if (theta_use.size() == 1) {
      theta_use = NumericVector(d, theta_use[0]);
    }
    if (static_cast<arma::uword>(theta_use.size()) != d) {
      stop("For isotropic=FALSE, 'theta' must be scalar or of length equal to ncol(X).");
    }

    arma::rowvec th(d);
    for (arma::uword j = 0; j < d; ++j) th[j] = theta_use[j];

    X_scaled  = Xm.each_row()  / th;
    Xp_scaled = Xpm.each_row() / th;
  }

  // ---- pairwise squared distances ----
  arma::mat dist_sq;
  if (symmetric) {
    arma::vec g = arma::sum(arma::square(X_scaled), 1); // n x 1
    arma::mat G = X_scaled * X_scaled.t();              // n x n
    dist_sq = arma::repmat(g, 1, g.n_elem) + arma::repmat(g.t(), g.n_elem, 1) - 2.0 * G;
  } else {
    arma::vec g  = arma::sum(arma::square(X_scaled), 1);   // n x 1
    arma::vec gp = arma::sum(arma::square(Xp_scaled), 1);  // m x 1
    arma::mat G  = X_scaled * Xp_scaled.t();               // n x m
    dist_sq = arma::repmat(g, 1, gp.n_elem) + arma::repmat(gp.t(), g.n_elem, 1) - 2.0 * G;
  }

  // numerical stability
  dist_sq.transform([](double v) { return (v < 0.0) ? 0.0 : v; });

  // ---- kernel evaluation ----
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
  } else { // wendland
    // L(x_a, x_b) = (q * r + 1) * max(0, 1-r)^q,
    // where r = ||(x_a - x_b)/theta||_2 and q = floor(d/2) + 3.
    const double q_w = std::floor(static_cast<double>(d) / 2.0) + 3.0;
    arma::mat dist = arma::sqrt(dist_sq);
    arma::mat one_minus = 1.0 - dist;
    one_minus.transform([](double v) { return (v > 0.0) ? v : 0.0; });
    K = (q_w * dist + 1.0) % arma::pow(one_minus, q_w);
  }

  return K;
}