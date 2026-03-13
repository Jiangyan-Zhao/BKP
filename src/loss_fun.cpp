// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double loss_fun_brier_bkp_rcpp(
    const arma::mat& K,
    const arma::vec& y,
    const arma::vec& m,
    const arma::vec& alpha0,
    const arma::vec& beta0
) {
  arma::vec alpha_n = alpha0 + K * y;
  arma::vec beta_n  = beta0 + K * (m - y);
  arma::vec pi_hat = alpha_n / (alpha_n + beta_n);
  arma::vec pi_tilde = y / m;
  double brier = arma::as_scalar(arma::mean(arma::square(pi_hat - pi_tilde)));
  return brier;
}

// [[Rcpp::export]]
double loss_fun_logloss_bkp_rcpp(
    const arma::mat& K,
    const arma::vec& y,
    const arma::vec& m,
    const arma::vec& alpha0,
    const arma::vec& beta0
) {
  arma::vec alpha_n = alpha0 + K * y;
  arma::vec beta_n  = beta0 + K * (m - y);
  arma::vec pi_hat = alpha_n / (alpha_n + beta_n);
  pi_hat = arma::clamp(pi_hat, 1e-10, 1.0 - 1e-10);
  arma::vec pi_tilde = y / m;
  double log_loss = -arma::as_scalar(arma::mean(pi_tilde % arma::log(pi_hat)));
  return log_loss;
}

// [[Rcpp::export]]
double loss_fun_brier_dkp_rcpp(
    const arma::mat& K,
    const arma::mat& Y,
    const arma::mat& alpha0
) {
  arma::mat alpha_n = alpha0 + K * Y;
  arma::vec row_sums = arma::sum(alpha_n, 1);
  arma::mat pi_hat = alpha_n.each_col() / row_sums;
  arma::vec Y_row_sums = arma::sum(Y, 1);
  arma::mat pi_tilde = Y.each_col() / Y_row_sums;
  double brier = arma::as_scalar(arma::mean(arma::vectorise(arma::square(pi_hat - pi_tilde))));
  return brier;
}

// [[Rcpp::export]]
double loss_fun_logloss_dkp_rcpp(
    const arma::mat& K,
    const arma::mat& Y,
    const arma::mat& alpha0
) {
  arma::mat alpha_n = alpha0 + K * Y;
  arma::vec row_sums = arma::sum(alpha_n, 1);
  arma::mat pi_hat = alpha_n.each_col() / row_sums;
  pi_hat.transform([](double x) { return std::max(std::min(x, 1.0 - 1e-10), 1e-10); });
  arma::vec Y_row_sums = arma::sum(Y, 1);
  arma::mat pi_tilde = Y.each_col() / Y_row_sums;
  double log_loss = -arma::as_scalar(arma::mean(arma::vectorise(pi_tilde % arma::log(pi_hat))));
  return log_loss;
}