// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat wendland_kernel_rcpp(
    const arma::mat& X1,
    const arma::mat& X2,
    const double theta,
    const int q_wend
) {
  if (theta <= 0.0) {
    stop("'theta' must be strictly positive.");
  }
  if (q_wend < 0) {
    stop("'q_wend' must be nonnegative.");
  }
  if (X1.n_cols != X2.n_cols) {
    stop("'X1' and 'X2' must have the same number of columns.");
  }

  const arma::uword n1 = X1.n_rows;
  const arma::uword n2 = X2.n_rows;
  arma::mat K(n1, n2, arma::fill::zeros);

  for (arma::uword i = 0; i < n1; ++i) {
    for (arma::uword j = 0; j < n2; ++j) {
      const double r = arma::norm(X1.row(i) - X2.row(j), 2);
      const double u = r / theta;
      const double one_minus_u = std::max(0.0, 1.0 - u);
      K(i, j) = (q_wend * u + 1.0) * std::pow(one_minus_u, q_wend);
    }
  }

  return K;
}
