// This file provides:
// - loss_bkp_arma(): internal pure Armadillo BKP loss engine
// - loss_dkp_arma(): internal pure Armadillo DKP loss engine
// - loss_fun_rcpp(): R-facing wrapper
//
// The *_arma functions intentionally avoid Rcpp objects so that they can be
// safely reused inside C++ optimization routines, including OpenMP-parallel
// multi-start optimization.
//
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// -----------------------------------------------------------------------------
// Internal C++ loss engines
// -----------------------------------------------------------------------------

double loss_bkp_arma(
    const std::string& loss,
    const arma::mat& K,
    const arma::vec& y,
    const arma::vec& m,
    const arma::vec& alpha0,
    const arma::vec& beta0,
    const arma::vec& data_scale
) {
  arma::vec Ky = K * y;
  arma::vec Km = K * m;

  arma::vec alpha_n = alpha0 + data_scale % Ky;
  arma::vec beta_n  = beta0 + data_scale % (Km - Ky);

  arma::vec pi_hat = alpha_n / (alpha_n + beta_n);
  arma::vec pi_tilde = y / m;

  if (loss == "brier") {

    return arma::as_scalar(
      arma::mean(arma::square(pi_hat - pi_tilde))
    );

  } else {  // loss == "log_loss"

    pi_hat = arma::clamp(pi_hat, 1e-10, 1.0 - 1e-10);

    return -arma::as_scalar(
        arma::mean(
          pi_tilde % arma::log(pi_hat) +
            (1.0 - pi_tilde) % arma::log(1.0 - pi_hat)
        )
    );
  }
}


double loss_dkp_arma(
    const std::string& loss,
    const arma::mat& K,
    const arma::mat& Y,
    const arma::mat& alpha0,
    const arma::vec& data_scale
) {
  arma::mat data_counts = K * Y;
  data_counts.each_col() %= data_scale;
  arma::mat alpha_n = alpha0 + data_counts;

  arma::vec alpha_row_sums = arma::sum(alpha_n, 1);
  arma::mat pi_hat = alpha_n.each_col() / alpha_row_sums;

  arma::vec Y_row_sums = arma::sum(Y, 1);
  arma::mat pi_tilde = Y.each_col() / Y_row_sums;

  if (loss == "brier") {

    // Keep current C++ definition: average over observations,
    // after summing squared errors over classes.
    return arma::accu(arma::square(pi_hat - pi_tilde)) / Y.n_rows;

  } else {  // loss == "log_loss"

    pi_hat = arma::clamp(pi_hat, 1e-10, 1.0 - 1e-10);

    return -arma::accu(pi_tilde % arma::log(pi_hat)) / Y.n_rows;
  }
}


// -----------------------------------------------------------------------------
// R-facing wrapper
//
// The R wrapper loss_fun() is responsible for input validation. This C++ wrapper
// only converts Rcpp objects to Armadillo objects and delegates computation to
// the internal *_arma engines.
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
double loss_fun_rcpp(
    std::string model,
    std::string loss,
    SEXP K,
    Nullable<NumericVector> y = R_NilValue,
    Nullable<NumericVector> m = R_NilValue,
    Nullable<NumericMatrix> Y = R_NilValue,
    Nullable<NumericVector> alpha0 = R_NilValue,
    Nullable<NumericVector> beta0 = R_NilValue,
    Nullable<NumericMatrix> alpha0_mat = R_NilValue,
    Nullable<NumericVector> data_scale = R_NilValue
) {
  NumericMatrix K_R(K);
  arma::mat K_mat(K_R.begin(), K_R.nrow(), K_R.ncol(), false);

  if (model == "BKP") {

    NumericVector y_R(y);
    NumericVector m_R(m);
    NumericVector alpha0_R(alpha0);
    NumericVector beta0_R(beta0);

    arma::vec y_vec(y_R.begin(), y_R.size(), false);
    arma::vec m_vec(m_R.begin(), m_R.size(), false);
    arma::vec alpha0_vec(alpha0_R.begin(), alpha0_R.size(), false);
    arma::vec beta0_vec(beta0_R.begin(), beta0_R.size(), false);

    if (data_scale.isNotNull()) {
      NumericVector data_scale_R(data_scale);
      arma::vec data_scale_vec(data_scale_R.begin(), data_scale_R.size(), false);
      return loss_bkp_arma(
        loss, K_mat, y_vec, m_vec, alpha0_vec, beta0_vec, data_scale_vec
      );
    }

    arma::vec data_scale_vec = arma::ones<arma::vec>(K_mat.n_rows);
    return loss_bkp_arma(
      loss, K_mat, y_vec, m_vec, alpha0_vec, beta0_vec, data_scale_vec
    );

  } else if (model == "DKP") {

    NumericMatrix Y_R(Y);
    NumericMatrix alpha0_R(alpha0_mat);

    arma::mat Y_mat(Y_R.begin(), Y_R.nrow(), Y_R.ncol(), false);
    arma::mat alpha0_dkp(alpha0_R.begin(), alpha0_R.nrow(), alpha0_R.ncol(), false);

    if (data_scale.isNotNull()) {
      NumericVector data_scale_R(data_scale);
      arma::vec data_scale_vec(data_scale_R.begin(), data_scale_R.size(), false);
      return loss_dkp_arma(loss, K_mat, Y_mat, alpha0_dkp, data_scale_vec);
    }

    arma::vec data_scale_vec = arma::ones<arma::vec>(K_mat.n_rows);
    return loss_dkp_arma(loss, K_mat, Y_mat, alpha0_dkp, data_scale_vec);

  } else {
    stop("Unsupported model: " + model);
  }
}
