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

    // Row-normalized kernel weights: W_ij = K_ij / rowSums(K)_i
    arma::mat W = K.each_col() / denom_K;

    arma::vec p_hat = W * (y / m);
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

    // Row-normalized kernel weights
    arma::mat W = K.each_col() / denom_K;

    arma::mat Pi_hat = W * Pi;
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

    arma::mat K_mat;

    if (K.isNotNull()) {
      NumericMatrix K_R(K);
      K_mat = as<arma::mat>(K_R);
    } else {
      // Preserve old behavior: when K is absent, prior length is 1.
      K_mat = arma::mat(1, 1, arma::fill::ones);
    }

    arma::vec y_vec;
    if (y.isNotNull()) {
      NumericVector y_R(y);
      y_vec = as<arma::vec>(y_R);
    }

    arma::vec m_vec;
    if (m.isNotNull()) {
      NumericVector m_R(m);
      m_vec = as<arma::vec>(m_R);
    }

    double p0_scalar = 0.5;
    if (p0.isNotNull()) {
      NumericVector p0_R(p0);
      p0_scalar = p0_R[0];
    }

    arma::vec alpha0;
    arma::vec beta0;

    get_prior_bkp_arma(
      prior,
      r0,
      p0_scalar,
      y_vec,
      m_vec,
      K_mat,
      alpha0,
      beta0
    );

    return List::create(
      Named("alpha0") = alpha0,
      Named("beta0")  = beta0
    );

  } else {  // model == "DKP"

    arma::mat K_mat;

    if (K.isNotNull()) {
      NumericMatrix K_R(K);
      K_mat = as<arma::mat>(K_R);
    } else {
      // Preserve old behavior: when K is absent, prior has one row.
      K_mat = arma::mat(1, 1, arma::fill::ones);
    }

    arma::mat Y_mat;
    if (Y.isNotNull()) {
      NumericMatrix Y_R(Y);
      Y_mat = as<arma::mat>(Y_R);
    }

    arma::vec p0_vec;
    if (p0.isNotNull()) {
      NumericVector p0_R(p0);
      p0_vec = as<arma::vec>(p0_R);
    }

    arma::mat alpha0 = get_prior_dkp_arma(
      prior,
      r0,
      p0_vec,
      Y_mat,
      K_mat
    );

    return List::create(
      Named("alpha0") = alpha0
    );
  }
}
