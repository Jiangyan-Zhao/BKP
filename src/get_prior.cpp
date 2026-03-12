// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// BKP noninformative prior
// [[Rcpp::export]]
List get_prior_bkp_noninformative_rcpp(int n) {
  arma::vec alpha0(n, arma::fill::ones);
  arma::vec beta0(n, arma::fill::ones);
  
  return List::create(
    Named("alpha0") = alpha0,
    Named("beta0") = beta0
  );
}

// BKP fixed prior
// [[Rcpp::export]]
List get_prior_bkp_fixed_rcpp(int n, double r0, double p0) {
  arma::vec alpha0(n);
  arma::vec beta0(n);
  alpha0.fill(r0 * p0);
  beta0.fill(r0 * (1.0 - p0));
  
  return List::create(
    Named("alpha0") = alpha0,
    Named("beta0") = beta0
  );
}

// BKP adaptive prior
// [[Rcpp::export]]
List get_prior_bkp_adaptive_rcpp(
    const arma::mat& K,
    const arma::vec& y,
    const arma::vec& m,
    double r0
) {
  arma::vec alpha0 = r0 * (K * (y / m));
  arma::vec beta0  = r0 * (K * ((m - y) / m));

  alpha0.elem(arma::find(alpha0 <= 1e-10)).zeros();
  beta0.elem(arma::find(beta0 <= 1e-10)).zeros();

  return List::create(
    Named("alpha0") = alpha0,
    Named("beta0") = beta0
  );
}

// DKP noninformative prior
// [[Rcpp::export]]
arma::mat get_prior_dkp_noninformative_rcpp(int n, int q) {
  arma::mat alpha0(n, q, arma::fill::ones);
  return alpha0;
}

// DKP fixed prior
// [[Rcpp::export]]
arma::mat get_prior_dkp_fixed_rcpp(int n, double r0, const arma::vec& p0) {
  const arma::uword q = p0.n_elem;
  
  arma::vec flat(n * q);
  arma::uword idx = 0;
  for (arma::uword j = 0; j < q; ++j) {
    for (int i = 0; i < n; ++i) {
      flat(idx++) = r0 * p0(j);
    }
  }

  arma::mat out(n, q);
  idx = 0;
  for (int i = 0; i < n; ++i) {
    for (arma::uword j = 0; j < q; ++j) {
      out(i, j) = flat(idx++);
    }
  }

  return out;
}

// DKP adaptive prior
// [[Rcpp::export]]
arma::mat get_prior_dkp_adaptive_rcpp(
    const arma::mat& K,
    const arma::mat& Y,
    double r0
) {
  // Row-normalized kernel weights
  arma::vec row_sum_K = arma::sum(K, 1);
  for (size_t i = 0; i < row_sum_K.n_elem; ++i) {
    if (row_sum_K(i) < 1e-10) row_sum_K(i) = 1.0;
  }
  arma::mat W = K.each_col() / row_sum_K;
  
  // Estimate local class proportions
  arma::vec row_sum_Y = arma::sum(Y, 1);
  arma::mat Pi_hat = W * (Y.each_col() / row_sum_Y);
  
  // Estimate local precision
  arma::vec r_hat = r0 * arma::sum(K, 1);
  r_hat = arma::clamp(r_hat, 1e-3, arma::datum::inf);
  
  // Compute prior parameters
  arma::mat alpha0 = Pi_hat.each_col() % r_hat;
  
  return alpha0;
}