// Input validation is handled in R/get_prior.R.
//
// This file provides:
// - get_prior_bkp_arma(): internal pure Armadillo BKP prior engine
// - get_prior_dkp_arma(): internal pure Armadillo DKP prior engine
// - get_prior_rcpp(): R-facing wrapper
//
// The *_arma functions intentionally avoid Rcpp objects so that they can be
// safely reused inside C++ optimization routines, including OpenMP-parallel
// multi-start optimization.
//
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// -----------------------------------------------------------------------------
// Internal C++ prior engines
// -----------------------------------------------------------------------------

void get_prior_bkp_arma(
    const std::string& prior,
    const double r0,
    const double p0,
    const arma::vec& y,
    const arma::vec& m,
    const arma::mat& K,
    arma::vec& alpha0,
    arma::vec& beta0
) {
  const arma::uword nrowK = K.n_rows;

  if (prior == "noninformative") {

    alpha0 = arma::vec(nrowK, arma::fill::ones);
    beta0  = arma::vec(nrowK, arma::fill::ones);

  } else if (prior == "fixed") {

    alpha0 = arma::vec(nrowK, arma::fill::value(r0 * p0));
    beta0  = arma::vec(nrowK, arma::fill::value(r0 * (1.0 - p0)));

  } else {  // prior == "adaptive"

    arma::vec row_sum_K = arma::sum(K, 1);

    // Equivalent to R's pmax(row_sum_K, 1e-6); no finite upper bound.
    arma::vec denom_K = arma::clamp(row_sum_K, 1e-6, arma::datum::inf);

    // Avoid materializing row-normalized weights W = K.each_col() / denom_K.
    // Algebraically equivalent to W * (y / m), with one fewer n x n
    // temporary matrix on the adaptive-prior hot path.
    arma::vec prop = y / m;
    arma::vec p_hat = (K * prop) / denom_K;
    arma::vec r_hat = r0 * arma::clamp(row_sum_K, 1e-3, arma::datum::inf);

    alpha0 = r_hat % p_hat;
    beta0  = r_hat % (1.0 - p_hat);

    alpha0 = arma::clamp(alpha0, 1e-2, arma::datum::inf);
    beta0  = arma::clamp(beta0,  1e-2, arma::datum::inf);
  }
}


arma::mat get_prior_dkp_arma(
    const std::string& prior,
    const double r0,
    const arma::vec& p0,
    const arma::mat& Y,
    const arma::mat& K
) {
  const arma::uword nrowK = K.n_rows;
  const arma::uword q = (Y.n_cols > 0) ? Y.n_cols : p0.n_elem;

  arma::mat alpha0;

  if (prior == "noninformative") {

    alpha0 = arma::mat(nrowK, q, arma::fill::ones);

  } else if (prior == "fixed") {

    alpha0 = arma::mat(nrowK, q, arma::fill::none);
    alpha0.each_row() = (r0 * p0).t();

  } else {  // prior == "adaptive"

    arma::vec row_sum_Y = arma::sum(Y, 1);
    arma::mat Pi = Y.each_col() / row_sum_Y;

    arma::vec row_sum_K = arma::sum(K, 1);
    arma::vec denom_K = arma::clamp(row_sum_K, 1e-6, arma::datum::inf);

    // Avoid materializing row-normalized weights W.  Compute K * Pi first,
    // then divide each row by the kernel row sum.
    arma::mat Pi_hat = K * Pi;
    Pi_hat.each_col() /= denom_K;
    arma::vec r_hat = r0 * arma::clamp(row_sum_K, 1e-3, arma::datum::inf);

    alpha0 = Pi_hat.each_col() % r_hat;
  }

  alpha0 = arma::clamp(alpha0, 1e-2, arma::datum::inf);

  return alpha0;
}


// -----------------------------------------------------------------------------
// R-facing wrapper
//
// The R wrapper get_prior() is responsible for input validation. This C++ wrapper
// only converts Rcpp objects to Armadillo objects and delegates computation to
// the internal *_arma engines.
//
// Important: K may be NULL for noninformative/fixed priors. In that case, the
// original implementation used nrowK = 1. We preserve that behavior here by
// using a 1 x 1 dummy kernel matrix.
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List get_prior_rcpp(
    std::string model,
    std::string prior,
    double r0,
    Nullable<NumericVector> p0 = R_NilValue,
    Nullable<NumericVector> y = R_NilValue,
    Nullable<NumericVector> m = R_NilValue,
    Nullable<NumericMatrix> Y = R_NilValue,
    Nullable<NumericMatrix> K = R_NilValue
) {
  if (model == "BKP") {

    NumericVector y_R;
    if (y.isNotNull()) y_R = NumericVector(y);
    arma::vec y_vec(y_R.begin(), y_R.size(), false);

    NumericVector m_R;
    if (m.isNotNull()) m_R = NumericVector(m);
    arma::vec m_vec(m_R.begin(), m_R.size(), false);

    double p0_scalar = 0.5;
    if (p0.isNotNull()) {
      NumericVector p0_R(p0);
      p0_scalar = p0_R[0];
    }

    arma::vec alpha0;
    arma::vec beta0;

    if (K.isNotNull()) {
      NumericMatrix K_R(K);
      arma::mat K_mat(K_R.begin(), K_R.nrow(), K_R.ncol(), false);
      get_prior_bkp_arma(prior, r0, p0_scalar, y_vec, m_vec, K_mat, alpha0, beta0);
    } else {
      // Preserve old behavior: when K is absent, prior length is 1.
      arma::mat K_mat(1, 1, arma::fill::ones);
      get_prior_bkp_arma(prior, r0, p0_scalar, y_vec, m_vec, K_mat, alpha0, beta0);
    }

    return List::create(
      Named("alpha0") = alpha0,
      Named("beta0")  = beta0
    );

  } else {  // model == "DKP"

    NumericMatrix Y_R;
    if (Y.isNotNull()) Y_R = NumericMatrix(Y);
    arma::mat Y_mat(Y_R.begin(), Y_R.nrow(), Y_R.ncol(), false);

    NumericVector p0_R;
    if (p0.isNotNull()) p0_R = NumericVector(p0);
    arma::vec p0_vec(p0_R.begin(), p0_R.size(), false);

    arma::mat alpha0;
    if (K.isNotNull()) {
      NumericMatrix K_R(K);
      arma::mat K_mat(K_R.begin(), K_R.nrow(), K_R.ncol(), false);
      alpha0 = get_prior_dkp_arma(prior, r0, p0_vec, Y_mat, K_mat);
    } else {
      // Preserve old behavior: when K is absent, prior has one row.
      arma::mat K_mat(1, 1, arma::fill::ones);
      alpha0 = get_prior_dkp_arma(prior, r0, p0_vec, Y_mat, K_mat);
    }

    return List::create(
      Named("alpha0") = alpha0
    );
  }
}
