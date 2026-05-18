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

  alpha0 = arma::clamp(alpha0, 1e-10, arma::datum::inf);
  beta0 = arma::clamp(beta0, 1e-10, arma::datum::inf);

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

  arma::mat alpha0(n, q, arma::fill::none);
  arma::rowvec row_val = (r0 * p0).t();
  alpha0.each_row() = row_val;
  return alpha0;
}

// DKP adaptive prior
// [[Rcpp::export]]
arma::mat get_prior_dkp_adaptive_rcpp(
    const arma::mat& K,
    const arma::mat& Y,
    double r0
) {
  arma::vec row_sum_Y = arma::sum(Y, 1);
  arma::mat P = Y.each_col() / row_sum_Y;

  arma::mat alpha0 = r0 * (K * P);
  alpha0 = arma::clamp(alpha0, 1e-10, arma::datum::inf);

  return alpha0;
}
