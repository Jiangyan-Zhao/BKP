// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List predict_bkp_rcpp(
    const arma::mat& K,
    const arma::vec& alpha0,
    const arma::vec& beta0,
    const arma::vec& y,
    const arma::vec& m
) {
  // Posterior parameters
  arma::vec alpha_n = alpha0 + K * y;
  arma::vec beta_n  = beta0 + K * (m - y);

  // Predictive mean and variance
  double eps = 1e-10;
  arma::vec pred_mean = alpha_n / arma::max(alpha_n + beta_n, eps * arma::ones(alpha_n.n_elem));
  pred_mean = arma::clamp(pred_mean, eps, 1.0 - eps);
  arma::vec pred_var = pred_mean % (1.0 - pred_mean) / (alpha_n + beta_n + 1.0);

  return List::create(
    Named("alpha_n")  = alpha_n,
    Named("beta_n")   = beta_n,
    Named("mean")     = pred_mean,
    Named("variance") = pred_var
  );
}