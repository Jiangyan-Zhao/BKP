// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List bkp_posterior_update_rcpp(
    const arma::mat& K,
    const arma::vec& y,
    const arma::vec& m,
    const arma::vec& alpha0,
    const arma::vec& beta0
) {
  arma::vec alpha_n = alpha0 + K * y;
  arma::vec beta_n  = beta0 + K * (m - y);

  return List::create(
    Named("alpha_n") = alpha_n,
    Named("beta_n")  = beta_n
  );
}