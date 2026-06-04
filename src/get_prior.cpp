// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List get_prior_rcpp(
    std::string model,
    std::string prior,
    double r0,
    Nullable<NumericVector> p0 = R_NilValue,
    Nullable<NumericVector> y = R_NilValue,
    Nullable<NumericVector> m = R_NilValue,
    Nullable<NumericMatrix> Y = R_NilValue,
    Nullable<NumericMatrix> K = R_NilValue
) {
  if (model == "BKP") {

    int nrowK = 1;
    if (K.isNotNull()) {
      NumericMatrix K_R(K);
      nrowK = K_R.nrow();
    }

    arma::vec alpha0;
    arma::vec beta0;

    if (prior == "noninformative") {

      alpha0 = arma::vec(nrowK, arma::fill::ones);
      beta0  = arma::vec(nrowK, arma::fill::ones);

    } else if (prior == "fixed") {

      NumericVector p0_R(p0);
      double p0_scalar = p0_R[0];

      alpha0 = arma::vec(nrowK, arma::fill::value(r0 * p0_scalar));
      beta0  = arma::vec(nrowK, arma::fill::value(r0 * (1.0 - p0_scalar)));

    } else if (prior == "adaptive") {

      NumericMatrix K_R(K);
      NumericVector y_R(y);
      NumericVector m_R(m);

      arma::mat K_mat = as<arma::mat>(K_R);
      arma::vec y_vec = as<arma::vec>(y_R);
      arma::vec m_vec = as<arma::vec>(m_R);

      arma::vec row_sum_K = arma::sum(K_mat, 1);
      // Equivalent to R's pmax(row_sum_K, 1e-6); no finite upper bound.
      arma::vec denom_K = arma::clamp(row_sum_K, 1e-6, arma::datum::inf);

      // Row-normalized kernel weights: W_ij = K_ij / rowSums(K)_i
      arma::mat W = K_mat.each_col() / denom_K;

      arma::vec p_hat = W * (y_vec / m_vec);
      arma::vec r_hat = r0 * arma::clamp(row_sum_K, 1e-3, arma::datum::inf);

      alpha0 = r_hat % p_hat;
      beta0  = r_hat % (1.0 - p_hat);

      alpha0 = arma::clamp(alpha0, 1e-2, arma::datum::inf);
      beta0  = arma::clamp(beta0,  1e-2, arma::datum::inf);
    }

    return List::create(
      Named("alpha0") = alpha0,
      Named("beta0")  = beta0
    );

  } else {  // model == "DKP"

    int nrowK = 1;
    if (K.isNotNull()) {
      NumericMatrix K_R(K);
      nrowK = K_R.nrow();
    }

    int q;
    if (Y.isNotNull()) {
      NumericMatrix Y_R(Y);
      q = Y_R.ncol();
    } else {
      NumericVector p0_R(p0);
      q = p0_R.size();
    }

    arma::mat alpha0;

    if (prior == "noninformative") {

      alpha0 = arma::mat(nrowK, q, arma::fill::ones);

    } else if (prior == "fixed") {

      NumericVector p0_R(p0);
      arma::vec p0_vec = as<arma::vec>(p0_R);

      alpha0 = arma::mat(nrowK, q, arma::fill::none);
      alpha0.each_row() = (r0 * p0_vec).t();

    } else if (prior == "adaptive") {

      NumericMatrix K_R(K);
      NumericMatrix Y_R(Y);

      arma::mat K_mat = as<arma::mat>(K_R);
      arma::mat Y_mat = as<arma::mat>(Y_R);

      arma::vec row_sum_Y = arma::sum(Y_mat, 1);
      arma::mat Pi = Y_mat.each_col() / row_sum_Y;

      arma::vec row_sum_K = arma::sum(K_mat, 1);
      arma::vec denom_K = arma::clamp(row_sum_K, 1e-6, arma::datum::inf);

      // Row-normalized kernel weights
      arma::mat W = K_mat.each_col() / denom_K;

      arma::mat Pi_hat = W * Pi;
      arma::vec r_hat = r0 * arma::clamp(row_sum_K, 1e-3, arma::datum::inf);

      alpha0 = Pi_hat.each_col() % r_hat;
    }

    alpha0 = arma::clamp(alpha0, 1e-2, arma::datum::inf);

    return List::create(
      Named("alpha0") = alpha0
    );
  }
}
