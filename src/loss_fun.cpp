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
  arma::vec alpha_n = alpha0 + data_scale % (K * y);
  arma::vec beta_n  = beta0 + data_scale % (K * (m - y));

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
    const arma::mat& alpha0
) {
  arma::mat alpha_n = alpha0 + K * Y;

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
    const arma::mat& K,
    Nullable<NumericVector> y = R_NilValue,
    Nullable<NumericVector> m = R_NilValue,
    Nullable<NumericMatrix> Y = R_NilValue,
    Nullable<NumericVector> alpha0 = R_NilValue,
    Nullable<NumericVector> beta0 = R_NilValue,
    Nullable<NumericMatrix> alpha0_mat = R_NilValue,
    Nullable<NumericVector> data_scale = R_NilValue
) {
  if (model == "BKP") {

    NumericVector y_R(y);
    NumericVector m_R(m);
    NumericVector alpha0_R(alpha0);
    NumericVector beta0_R(beta0);

    arma::vec y_vec = as<arma::vec>(y_R);
    arma::vec m_vec = as<arma::vec>(m_R);
    arma::vec alpha0_vec = as<arma::vec>(alpha0_R);
    arma::vec beta0_vec = as<arma::vec>(beta0_R);
    arma::vec data_scale_vec;

    if (data_scale.isNotNull()) {
      NumericVector data_scale_R(data_scale);
      data_scale_vec = as<arma::vec>(data_scale_R);
    } else {
      data_scale_vec = arma::ones<arma::vec>(K.n_rows);
    }

    return loss_bkp_arma(
      loss,
      K,
      y_vec,
      m_vec,
      alpha0_vec,
      beta0_vec,
      data_scale_vec
    );

  } else if (model == "DKP") {

    NumericMatrix Y_R(Y);
    NumericMatrix alpha0_R(alpha0_mat);

    arma::mat Y_mat = as<arma::mat>(Y_R);
    arma::mat alpha0_dkp = as<arma::mat>(alpha0_R);

    return loss_dkp_arma(
      loss,
      K,
      Y_mat,
      alpha0_dkp
    );

  } else {
    stop("Unsupported model: " + model);
  }
}
