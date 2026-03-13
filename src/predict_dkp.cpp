// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List predict_dkp_rcpp(
    const arma::mat& K,
    const arma::mat& alpha0,
    const arma::mat& Y
) {
  // Posterior parameters: alpha_n = alpha0 + K %*% Y
  arma::mat alpha_n = alpha0 + K * Y;

  // Row sums (total concentration)
  arma::vec row_sum = arma::sum(alpha_n, 1);

  // Predictive mean: each row divided by its row sum
  arma::mat pred_mean(alpha_n.n_rows, alpha_n.n_cols);
  for (size_t i = 0; i < alpha_n.n_rows; ++i) {
    pred_mean.row(i) = alpha_n.row(i) / std::max(row_sum(i), 1e-10);
  }
  pred_mean = arma::clamp(pred_mean, 1e-10, 1.0 - 1e-10);
  
  // Predictive variance
  arma::mat pred_var = pred_mean % (1.0 - pred_mean);
  for (size_t i = 0; i < alpha_n.n_rows; ++i) {
    pred_var.row(i) = pred_var.row(i) / (row_sum(i) + 1.0);
  }

  return List::create(
    Named("alpha_n")  = alpha_n,
    Named("row_sum")  = row_sum,
    Named("mean")     = pred_mean,
    Named("variance") = pred_var
  );
}