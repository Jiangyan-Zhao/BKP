// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <string>

using namespace Rcpp;

namespace {

inline double kernel_from_dist_sq(double dist_sq, const std::string& kernel, int d) {
  if (kernel == "gaussian") return std::exp(-dist_sq);
  const double r = std::sqrt(std::max(0.0, dist_sq));
  if (kernel == "matern52") {
    const double sr = std::sqrt(5.0) * r;
    return (1.0 + sr + 5.0 * r * r / 3.0) * std::exp(-sr);
  }
  if (kernel == "matern32") {
    const double sr = std::sqrt(3.0) * r;
    return (1.0 + sr) * std::exp(-sr);
  }
  if (kernel == "wendland") {
    const int q = static_cast<int>(std::floor(static_cast<double>(d) / 2.0)) + 3;
    const double one_minus_r = std::max(0.0, 1.0 - r);
    return (static_cast<double>(q) * r + 1.0) * std::pow(one_minus_r, q);
  }
  Rcpp::stop("Unsupported kernel.");
  return NA_REAL;
}

inline double dist_sq_scaled_global(const NumericMatrix& Xq, const NumericMatrix& Xt,
                                    std::size_t qi, std::size_t tj,
                                    const NumericVector& theta_g, bool isotropic) {
  const std::size_t d = Xq.ncol();
  double out = 0.0;
  if (isotropic) {
    const double th = theta_g[0];
    for (std::size_t k = 0; k < d; k++) {
      const double diff = (Xq(qi, k) - Xt(tj, k)) / th;
      out += diff * diff;
    }
  } else {
    for (std::size_t k = 0; k < d; k++) {
      const double diff = (Xq(qi, k) - Xt(tj, k)) / theta_g[k];
      out += diff * diff;
    }
  }
  return out;
}

inline double dist_sq_scaled_local(const NumericMatrix& Xq, const NumericMatrix& Xt,
                                   std::size_t qi, std::size_t tj, double theta_l) {
  const std::size_t d = Xq.ncol();
  double out = 0.0;
  for (std::size_t k = 0; k < d; k++) {
    const double diff = (Xq(qi, k) - Xt(tj, k)) / theta_l;
    out += diff * diff;
  }
  return out;
}

inline void add_weighted_row(NumericMatrix& accum, NumericVector& sum_w_pi,
                             double& sum_w, const NumericMatrix& Y,
                             std::size_t qi, std::size_t j, double w) {
  const std::size_t q = Y.ncol();
  double y_row_sum = 0.0;
  for (std::size_t c = 0; c < q; c++) {
    const double yjc = Y(j, c);
    accum(qi, c) += w * yjc;
    y_row_sum += yjc;
  }
  sum_w += w;
  if (y_row_sum > 0.0) {
    for (std::size_t c = 0; c < q; c++) sum_w_pi[c] += w * (Y(j, c) / y_row_sum);
  }
}

} // namespace

// [[Rcpp::export]]
Rcpp::List twin_dkp_posterior_rcpp(
    Rcpp::NumericMatrix Xquery_norm,
    Rcpp::NumericMatrix Xtrain_norm,
    Rcpp::NumericMatrix Y,
    Rcpp::IntegerVector g_indices,
    Rcpp::IntegerMatrix local_indices,
    Rcpp::NumericVector theta_g,
    double theta_l,
    std::string global_kernel,
    std::string local_kernel,
    bool isotropic,
    std::string prior,
    double r0,
    Rcpp::NumericVector p0,
    bool store_kernel = false)
{
  const std::size_t t = Xquery_norm.nrow();
  const std::size_t n = Xtrain_norm.nrow();
  const int d = Xquery_norm.ncol();
  const std::size_t q = Y.ncol();
  const std::size_t g = g_indices.size();
  const std::size_t l = local_indices.ncol();

  NumericMatrix alpha0(t, q), alpha_n(t, q), prob(t, q), data_update(t, q);
  NumericMatrix K, K_global, K_local;
  if (store_kernel) {
    K = NumericMatrix(t, n);
    K_global = NumericMatrix(t, n);
    K_local = NumericMatrix(t, n);
  }

  for (std::size_t qi = 0; qi < t; qi++) {
    NumericVector sum_w_pi(q);
    double sum_w = 0.0;

    for (std::size_t kk = 0; kk < g; kk++) {
      const int idx = g_indices[kk];
      if (idx <= 0) continue;
      const std::size_t j = static_cast<std::size_t>(idx - 1);
      const double w = kernel_from_dist_sq(
        dist_sq_scaled_global(Xquery_norm, Xtrain_norm, qi, j, theta_g, isotropic),
        global_kernel, d);
      add_weighted_row(data_update, sum_w_pi, sum_w, Y, qi, j, w);
      if (store_kernel) K_global(qi, j) = w;
    }

    for (std::size_t kk = 0; kk < l; kk++) {
      const int idx = local_indices(qi, kk);
      if (idx <= 0) continue;
      const std::size_t j = static_cast<std::size_t>(idx - 1);
      const double w = kernel_from_dist_sq(
        dist_sq_scaled_local(Xquery_norm, Xtrain_norm, qi, j, theta_l),
        local_kernel, d);
      add_weighted_row(data_update, sum_w_pi, sum_w, Y, qi, j, w);
      if (store_kernel) K_local(qi, j) = w;
    }

    double row_sum_alpha = 0.0;
    for (std::size_t c = 0; c < q; c++) {
      double a0 = 1.0;
      if (prior == "fixed") {
        a0 = r0 * p0[c];
      } else if (prior == "adaptive") {
        const double pi_hat = sum_w_pi[c] / std::max(sum_w, 1e-6);
        const double r_hat = r0 * std::max(sum_w, 1e-3);
        a0 = std::max(1e-2, r_hat * pi_hat);
      } else if (prior != "noninformative") {
        Rcpp::stop("Unsupported prior.");
      }
      if (!R_finite(a0) || a0 <= 0.0) a0 = 1e-2;
      alpha0(qi, c) = a0;
      double an = a0 + data_update(qi, c);
      if (!R_finite(an) || an <= 0.0) Rcpp::stop("Posterior concentration parameters must be positive and finite.");
      alpha_n(qi, c) = an;
      row_sum_alpha += an;
    }
    if (!R_finite(row_sum_alpha) || row_sum_alpha <= 0.0) Rcpp::stop("Posterior row sums must be positive and finite.");
    for (std::size_t c = 0; c < q; c++) prob(qi, c) = alpha_n(qi, c) / row_sum_alpha;

    if (store_kernel) {
      for (std::size_t j = 0; j < n; j++) K(qi, j) = K_global(qi, j) + K_local(qi, j);
    }
  }

  Rcpp::RObject K_out = store_kernel ? Rcpp::wrap(K) : R_NilValue;
  Rcpp::RObject K_global_out = store_kernel ? Rcpp::wrap(K_global) : R_NilValue;
  Rcpp::RObject K_local_out = store_kernel ? Rcpp::wrap(K_local) : R_NilValue;

  return List::create(
    _["alpha0"] = alpha0,
    _["alpha_n"] = alpha_n,
    _["prob"] = prob,
    _["K"] = K_out,
    _["K_global"] = K_global_out,
    _["K_local"] = K_local_out
  );
}
