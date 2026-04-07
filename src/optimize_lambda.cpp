// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>

using namespace Rcpp;

double loss_fun_brier_bkp_rcpp(
  const arma::mat& K, const arma::vec& y, const arma::vec& m,
  const arma::vec& alpha0, const arma::vec& beta0
);

double loss_fun_logloss_bkp_rcpp(
  const arma::mat& K, const arma::vec& y, const arma::vec& m,
  const arma::vec& alpha0, const arma::vec& beta0
);

double loss_fun_brier_dkp_rcpp(
  const arma::mat& K, const arma::mat& Y, const arma::mat& alpha0
);

double loss_fun_logloss_dkp_rcpp(
  const arma::mat& K, const arma::mat& Y, const arma::mat& alpha0
);

static inline double eval_bkp_lambda_loss(
    const double lambda,
    const arma::mat& K_g,
    const arma::mat& K_l,
    const arma::vec& y,
    const arma::vec& m,
    const arma::vec& alpha0,
    const arma::vec& beta0,
    const std::string& loss
) {
  const double lam = std::max(0.0, std::min(1.0, lambda));
  const arma::mat K_mix = lam * K_g + (1.0 - lam) * K_l;

  if (loss == "brier") {
    return loss_fun_brier_bkp_rcpp(K_mix, y, m, alpha0, beta0);
  }
  if (loss == "log_loss") {
    return loss_fun_logloss_bkp_rcpp(K_mix, y, m, alpha0, beta0);
  }
  stop("Unsupported loss in optimize_lambda_bkp_rcpp: " + loss);
}

static inline double eval_dkp_lambda_loss(
    const double lambda,
    const arma::mat& K_g,
    const arma::mat& K_l,
    const arma::mat& Y,
    const arma::mat& alpha0,
    const std::string& loss
) {
  const double lam = std::max(0.0, std::min(1.0, lambda));
  const arma::mat K_mix = lam * K_g + (1.0 - lam) * K_l;

  if (loss == "brier") {
    return loss_fun_brier_dkp_rcpp(K_mix, Y, alpha0);
  }
  if (loss == "log_loss") {
    return loss_fun_logloss_dkp_rcpp(K_mix, Y, alpha0);
  }
  stop("Unsupported loss in optimize_lambda_dkp_rcpp: " + loss);
}

// [[Rcpp::export]]
Rcpp::List optimize_lambda_bkp_rcpp(
    const arma::mat& K_g,
    const arma::mat& K_l,
    const arma::vec& y,
    const arma::vec& m,
    const arma::vec& alpha0,
    const arma::vec& beta0,
    const std::string& loss,
    const int max_iter = 80,
    const double tol = 1e-8
) {
  if (K_g.n_rows != K_l.n_rows || K_g.n_cols != K_l.n_cols) {
    stop("'K_g' and 'K_l' must have the same dimensions.");
  }

  double a = 0.0, b = 1.0;
  const double phi = (std::sqrt(5.0) - 1.0) / 2.0;
  double x1 = b - phi * (b - a);
  double x2 = a + phi * (b - a);
  double f1 = eval_bkp_lambda_loss(x1, K_g, K_l, y, m, alpha0, beta0, loss);
  double f2 = eval_bkp_lambda_loss(x2, K_g, K_l, y, m, alpha0, beta0, loss);

  for (int it = 0; it < max_iter && (b - a) > tol; ++it) {
    if (f1 > f2) {
      a = x1;
      x1 = x2;
      f1 = f2;
      x2 = a + phi * (b - a);
      f2 = eval_bkp_lambda_loss(x2, K_g, K_l, y, m, alpha0, beta0, loss);
    } else {
      b = x2;
      x2 = x1;
      f2 = f1;
      x1 = b - phi * (b - a);
      f1 = eval_bkp_lambda_loss(x1, K_g, K_l, y, m, alpha0, beta0, loss);
    }
  }

  const double lambda_opt = 0.5 * (a + b);
  const double loss_opt = eval_bkp_lambda_loss(lambda_opt, K_g, K_l, y, m, alpha0, beta0, loss);

  return Rcpp::List::create(
    Rcpp::Named("lambda_opt") = lambda_opt,
    Rcpp::Named("loss_opt") = loss_opt
  );
}

// [[Rcpp::export]]
Rcpp::List optimize_lambda_dkp_rcpp(
    const arma::mat& K_g,
    const arma::mat& K_l,
    const arma::mat& Y,
    const arma::mat& alpha0,
    const std::string& loss,
    const int max_iter = 80,
    const double tol = 1e-8
) {
  if (K_g.n_rows != K_l.n_rows || K_g.n_cols != K_l.n_cols) {
    stop("'K_g' and 'K_l' must have the same dimensions.");
  }

  double a = 0.0, b = 1.0;
  const double phi = (std::sqrt(5.0) - 1.0) / 2.0;
  double x1 = b - phi * (b - a);
  double x2 = a + phi * (b - a);
  double f1 = eval_dkp_lambda_loss(x1, K_g, K_l, Y, alpha0, loss);
  double f2 = eval_dkp_lambda_loss(x2, K_g, K_l, Y, alpha0, loss);

  for (int it = 0; it < max_iter && (b - a) > tol; ++it) {
    if (f1 > f2) {
      a = x1;
      x1 = x2;
      f1 = f2;
      x2 = a + phi * (b - a);
      f2 = eval_dkp_lambda_loss(x2, K_g, K_l, Y, alpha0, loss);
    } else {
      b = x2;
      x2 = x1;
      f2 = f1;
      x1 = b - phi * (b - a);
      f1 = eval_dkp_lambda_loss(x1, K_g, K_l, Y, alpha0, loss);
    }
  }

  const double lambda_opt = 0.5 * (a + b);
  const double loss_opt = eval_dkp_lambda_loss(lambda_opt, K_g, K_l, Y, alpha0, loss);

  return Rcpp::List::create(
    Rcpp::Named("lambda_opt") = lambda_opt,
    Rcpp::Named("loss_opt") = loss_opt
  );
}
