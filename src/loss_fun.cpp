// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

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
    Nullable<NumericMatrix> alpha0_mat = R_NilValue
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

    arma::vec alpha_n = alpha0_vec + K * y_vec;
    arma::vec beta_n  = beta0_vec + K * (m_vec - y_vec);

    arma::vec pi_hat = alpha_n / (alpha_n + beta_n);
    arma::vec pi_tilde = y_vec / m_vec;

    if (loss == "brier") {

      return arma::as_scalar(
        arma::mean(arma::square(pi_hat - pi_tilde))
      );

    } else if (loss == "log_loss") {

      pi_hat = arma::clamp(pi_hat, 1e-10, 1.0 - 1e-10);

      return -arma::as_scalar(
          arma::mean(
            pi_tilde % arma::log(pi_hat) +
              (1.0 - pi_tilde) % arma::log(1.0 - pi_hat)
          )
      );

    } else {
      stop("Unsupported loss: " + loss);
    }

  } else if (model == "DKP") {

    NumericMatrix Y_R(Y);
    NumericMatrix alpha0_R(alpha0_mat);

    arma::mat Y_mat = as<arma::mat>(Y_R);
    arma::mat alpha0_dkp = as<arma::mat>(alpha0_R);

    arma::mat alpha_n = alpha0_dkp + K * Y_mat;

    arma::vec alpha_row_sums = arma::sum(alpha_n, 1);
    arma::mat pi_hat = alpha_n.each_col() / alpha_row_sums;

    arma::vec Y_row_sums = arma::sum(Y_mat, 1);
    arma::mat pi_tilde = Y_mat.each_col() / Y_row_sums;

    if (loss == "brier") {

      // Keep current C++ definition: average over observations,
      // after summing squared errors over classes.
      return arma::accu(arma::square(pi_hat - pi_tilde)) / Y_mat.n_rows;

    } else if (loss == "log_loss") {

      pi_hat = arma::clamp(pi_hat, 1e-10, 1.0 - 1e-10);

      return -arma::accu(pi_tilde % arma::log(pi_hat)) / Y_mat.n_rows;

    } else {
      stop("Unsupported loss: " + loss);
    }

  } else {
    stop("Unsupported model: " + model);
  }
}
