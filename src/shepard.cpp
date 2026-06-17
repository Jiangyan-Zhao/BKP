// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;

static double squared_distance_row(
    const arma::mat& A,
    const arma::mat& B,
    const arma::uword i,
    const arma::uword j
) {
  double out = 0.0;
  for (arma::uword k = 0; k < A.n_cols; ++k) {
    const double diff = A(i, k) - B(j, k);
    out += diff * diff;
  }
  return out;
}

// [[Rcpp::export]]
arma::vec shepard_m_rcpp(
    const arma::mat& Xquery_norm,
    const arma::mat& Xtrain_norm,
    const arma::vec& m,
    const double power = 2.0
) {
  const arma::uword n_query = Xquery_norm.n_rows;
  const arma::uword n_train = Xtrain_norm.n_rows;

  arma::vec out(n_query);

  for (arma::uword i = 0; i < n_query; ++i) {
    if (i % 1024 == 0) Rcpp::checkUserInterrupt();

    bool has_exact = false;
    double numerator = 0.0;
    double denominator = 0.0;

    for (arma::uword j = 0; j < n_train; ++j) {
      const double dist_sq = squared_distance_row(Xquery_norm, Xtrain_norm, i, j);

      if (dist_sq == 0.0) {
        out[i] = m[j];
        has_exact = true;
        break;
      }

      double weight;
      if (std::abs(power - 2.0) < 1e-9) {
        weight = 1.0 / dist_sq;
      } else {
        weight = std::pow(dist_sq, -0.5 * power);
      }
      numerator += weight * m[j];
      denominator += weight;
    }

    if (!has_exact) {
      out[i] = numerator / denominator;
    }
  }

  return out;
}

// [[Rcpp::export]]
arma::vec shepard_m_loo_rcpp(
    const arma::mat& Xnorm,
    const arma::vec& m,
    const double power = 2.0
) {
  const arma::uword n = Xnorm.n_rows;
  arma::vec out(n);

  for (arma::uword i = 0; i < n; ++i) {
    if (i % 1024 == 0) Rcpp::checkUserInterrupt();

    double numerator = 0.0;
    double denominator = 0.0;

    for (arma::uword j = 0; j < n; ++j) {
      if (i == j) continue;

      const double dist_sq = squared_distance_row(Xnorm, Xnorm, i, j);
      double weight;
      if (std::abs(power - 2.0) < 1e-9) {
        weight = 1.0 / dist_sq;
      } else {
        weight = std::pow(dist_sq, -0.5 * power);
      }

      numerator += weight * m[j];
      denominator += weight;
    }

    out[i] = numerator / denominator;
  }

  return out;
}
