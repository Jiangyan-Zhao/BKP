#include <Rcpp.h>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

inline double lbeta_cpp(double a, double b) {
  return R::lbeta(a, b);
}

inline int qbetabinom_one(double prob, int size, double alpha, double beta) {
  if (!R_finite(prob) || !R_finite(alpha) || !R_finite(beta) ||
      size < 0 || alpha <= 0.0 || beta <= 0.0) {
    return NA_INTEGER;
  }

  if (prob < 0.0 || prob > 1.0) return NA_INTEGER;
  if (prob == 0.0) return 0;
  if (prob == 1.0) return size;
  if (size == 0) return 0;

  // log P(Y = 0)
  const double logp0 = lbeta_cpp(alpha, size + beta) - lbeta_cpp(alpha, beta);
  if (!R_finite(logp0)) return NA_INTEGER;

  // Pass 1: find max log-probability for log-sum-exp stabilization
  double logpk = logp0;
  double max_logp = logpk;

  for (int k = 0; k < size; ++k) {
    logpk += std::log(static_cast<double>(size - k))
    - std::log(static_cast<double>(k + 1))
    + std::log(alpha + k)
    - std::log(beta + size - k - 1.0);

    if (logpk > max_logp) max_logp = logpk;
  }

  // Pass 2: compute stabilized total mass
  double total = 0.0;
  logpk = logp0;

  for (int k = 0; k <= size; ++k) {
    total += std::exp(logpk - max_logp);

    if (k < size) {
      logpk += std::log(static_cast<double>(size - k))
      - std::log(static_cast<double>(k + 1))
      + std::log(alpha + k)
      - std::log(beta + size - k - 1.0);
    }
  }

  const double target = prob * total;

  // Pass 3: find the first k such that F(k) >= prob
  double cdf = 0.0;
  logpk = logp0;

  for (int k = 0; k <= size; ++k) {
    cdf += std::exp(logpk - max_logp);

    if (cdf >= target) return k;

    if (k < size) {
      logpk += std::log(static_cast<double>(size - k))
      - std::log(static_cast<double>(k + 1))
      + std::log(alpha + k)
      - std::log(beta + size - k - 1.0);
    }
  }

  return size;
}

// [[Rcpp::export]]
IntegerVector qbetabinom_rcpp(
    NumericVector prob,
    NumericVector size,
    NumericVector alpha,
    NumericVector beta
) {
  const R_xlen_t n_prob  = prob.size();
  const R_xlen_t n_size  = size.size();
  const R_xlen_t n_alpha = alpha.size();
  const R_xlen_t n_beta  = beta.size();

  const R_xlen_t n = std::max(
    std::max(n_prob, n_size),
    std::max(n_alpha, n_beta)
  );

  if (!((n_prob == 1 || n_prob == n) &&
      (n_size == 1 || n_size == n) &&
      (n_alpha == 1 || n_alpha == n) &&
      (n_beta == 1 || n_beta == n))) {
    stop("'prob', 'size', 'alpha', and 'beta' must have length 1 or the same length.");
  }

  IntegerVector out(n);

  for (R_xlen_t i = 0; i < n; ++i) {
    if (i % 1000 == 0) Rcpp::checkUserInterrupt();

    const double pi = prob[n_prob == 1 ? 0 : i];
    const double si = size[n_size == 1 ? 0 : i];
    const double ai = alpha[n_alpha == 1 ? 0 : i];
    const double bi = beta[n_beta == 1 ? 0 : i];

    if (!R_finite(si) || si < 0.0 || std::floor(si) != si) {
      out[i] = NA_INTEGER;
    } else {
      out[i] = qbetabinom_one(pi, static_cast<int>(si), ai, bi);
    }
  }

  return out;
}
