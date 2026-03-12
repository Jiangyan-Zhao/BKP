// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List dkp_posterior_update_rcpp(
    const arma::mat& K,
    const arma::mat& Y,
    const arma::mat& alpha0
) {
  arma::mat alpha_n = alpha0 + K * Y;

  return List::create(
    Named("alpha_n") = alpha_n
  );
}