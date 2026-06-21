// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <string>

using namespace Rcpp;

namespace {

inline double clamp_prob(double x) {
  return std::min(1.0 - 1e-12, std::max(1e-12, x));
}

inline double kernel_from_dist_sq(double dist_sq, const std::string& kernel, int d) {
  if (kernel == "gaussian") {
    return std::exp(-dist_sq);
  }

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

} // namespace

// [[Rcpp::export]]
Rcpp::List twin_bkp_posterior_rcpp(
    Rcpp::NumericMatrix Xquery_norm,
    Rcpp::NumericMatrix Xtrain_norm,
    Rcpp::NumericVector y,
    Rcpp::NumericVector m,
    Rcpp::IntegerVector g_indices,
    Rcpp::IntegerMatrix local_indices,
    Rcpp::NumericVector theta_g,
    double theta_l,
    std::string global_kernel,
    std::string local_kernel,
    bool isotropic,
    std::string prior,
    double r0,
    double p0,
    std::string ess,
    Rcpp::Nullable<Rcpp::NumericVector> m_shepard = R_NilValue,
    bool store_kernel = false)
{
  const std::size_t t = Xquery_norm.nrow();
  const std::size_t n = Xtrain_norm.nrow();
  const int d = Xquery_norm.ncol();
  const std::size_t g = g_indices.size();
  const std::size_t l = local_indices.ncol();

  NumericVector alpha0(t), beta0(t), alpha_n(t), beta_n(t);
  NumericVector scale(t), m_kernel_vec(t), rho_vec(t);
  NumericVector m_target_vec;
  NumericVector m_shepard_vec;
  const bool use_shepard = (ess == "shepard");

  if (use_shepard) {
    m_shepard_vec = Rcpp::as<Rcpp::NumericVector>(m_shepard);
    m_target_vec = NumericVector(t);
  }

  NumericMatrix K, K_global, K_local;
  if (store_kernel) {
    K = NumericMatrix(t, n);
    K_global = NumericMatrix(t, n);
    K_local = NumericMatrix(t, n);
  }

  for (std::size_t qi = 0; qi < t; qi++) {
    double success = 0.0;
    double failure = 0.0;
    double m_kernel = 0.0;
    double rho = 0.0;
    double sum_w = 0.0;
    double sum_w_prop = 0.0;

    for (std::size_t kk = 0; kk < g; kk++) {
      const std::size_t j = static_cast<std::size_t>(g_indices[kk] - 1);
      const double w = kernel_from_dist_sq(
        dist_sq_scaled_global(Xquery_norm, Xtrain_norm, qi, j, theta_g, isotropic),
        global_kernel, d);
      success += w * y[j];
      failure += w * (m[j] - y[j]);
      m_kernel += w * m[j];
      rho = std::max(rho, w);
      sum_w += w;
      sum_w_prop += w * (y[j] / m[j]);
      if (store_kernel) K_global(qi, j) = w;
    }

    for (std::size_t kk = 0; kk < l; kk++) {
      const std::size_t j = static_cast<std::size_t>(local_indices(qi, kk) - 1);
      const double w = kernel_from_dist_sq(
        dist_sq_scaled_local(Xquery_norm, Xtrain_norm, qi, j, theta_l),
        local_kernel, d);
      success += w * y[j];
      failure += w * (m[j] - y[j]);
      m_kernel += w * m[j];
      rho = std::max(rho, w);
      sum_w += w;
      sum_w_prop += w * (y[j] / m[j]);
      if (store_kernel) K_local(qi, j) = w;
    }

    if (prior == "noninformative") {
      alpha0[qi] = 1.0;
      beta0[qi] = 1.0;
    } else if (prior == "fixed") {
      alpha0[qi] = r0 * p0;
      beta0[qi] = r0 * (1.0 - p0);
    } else if (prior == "adaptive") {
      if (sum_w > 0.0) {
        const double p_adapt = clamp_prob(sum_w_prop / sum_w);
        const double r_adapt = r0 * sum_w;
        alpha0[qi] = r_adapt * p_adapt;
        beta0[qi] = r_adapt * (1.0 - p_adapt);
      } else {
        const double p_adapt = clamp_prob(p0);
        const double r_adapt = 0.0;
        alpha0[qi] = r_adapt * p_adapt;
        beta0[qi] = r_adapt * (1.0 - p_adapt);
      }
    }

    double sc = 1.0;
    double target = NA_REAL;
    if (use_shepard) {
      target = rho * m_shepard_vec[qi];
      sc = (m_kernel > 0.0) ? target / m_kernel : 0.0;
      m_target_vec[qi] = target;
    }

    scale[qi] = sc;
    m_kernel_vec[qi] = m_kernel;
    rho_vec[qi] = rho;
    alpha_n[qi] = alpha0[qi] + sc * success;
    beta_n[qi] = beta0[qi] + sc * failure;

    if (store_kernel) {
      for (std::size_t j = 0; j < n; j++) {
        K(qi, j) = K_global(qi, j) + K_local(qi, j);
      }
    }
  }

  Rcpp::RObject m_shepard_out = use_shepard ? Rcpp::wrap(m_shepard_vec) : R_NilValue;
  Rcpp::RObject m_target_out = use_shepard ? Rcpp::wrap(m_target_vec) : R_NilValue;
  Rcpp::RObject K_out = store_kernel ? Rcpp::wrap(K) : R_NilValue;
  Rcpp::RObject K_global_out = store_kernel ? Rcpp::wrap(K_global) : R_NilValue;
  Rcpp::RObject K_local_out = store_kernel ? Rcpp::wrap(K_local) : R_NilValue;

  List ess_info = List::create(
    _["scale"] = scale,
    _["m_kernel"] = m_kernel_vec,
    _["rho"] = rho_vec,
    _["m_shepard"] = m_shepard_out,
    _["m_target"] = m_target_out
  );

  return List::create(
    _["alpha0"] = alpha0,
    _["beta0"] = beta0,
    _["alpha_n"] = alpha_n,
    _["beta_n"] = beta_n,
    _["ess_info"] = ess_info,
    _["K"] = K_out,
    _["K_global"] = K_global_out,
    _["K_local"] = K_local_out
  );
}
