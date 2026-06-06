// -----------------------------------------------------------------------------
// This file implements kernel hyperparameter optimization for BKP and DKP.
//
// This file contains the shared optimization workflow for both Beta Kernel
// Process and Dirichlet Kernel Process models. The BKP and DKP routines are
// kept together because they use the same gamma-scale parameterization,
// candidate generation, and nloptr-based refinement strategy.
// -----------------------------------------------------------------------------
#include <RcppArmadillo.h>
#include <nloptrAPI.h>
#include <limits>
#include <cmath>
#include <algorithm>
#include <vector>

// [[Rcpp::depends(RcppArmadillo, nloptr)]]

using namespace Rcpp;

// -------- BEGIN: declarations from other cpp files --------
arma::mat kernel_matrix_arma(
    const arma::mat& Xm,
    const arma::mat& Xpm,
    const arma::vec& theta,
    const std::string& kernel,
    const bool isotropic,
    const bool symmetric
);

List get_prior_rcpp(
    std::string model,
    std::string prior,
    double r0,
    Nullable<NumericVector> p0,
    Nullable<NumericVector> y,
    Nullable<NumericVector> m,
    Nullable<NumericMatrix> Y,
    Nullable<NumericMatrix> K
);

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
);
// -------- END: declarations from other cpp files --------

static double eval_bkp_loss_from_gamma(
    const arma::vec& gamma,
    const arma::mat& Xnorm,
    const arma::vec& y,
    const arma::vec& m,
    const std::string& prior,
    const double r0,
    const double p0,
    const std::string& loss,
    const std::string& kernel,
    const bool isotropic
) {
  arma::vec theta = arma::exp(std::log(10.0) * gamma);

  arma::mat K = kernel_matrix_arma(Xnorm, Xnorm, theta, kernel, isotropic, true);
  K.diag().zeros();

  NumericVector p0_vec = NumericVector::create(p0);
  NumericVector y_vec = wrap(y);
  NumericVector m_vec = wrap(m);
  NumericMatrix K_mat = wrap(K);

  List prior_par = get_prior_rcpp(
    "BKP", prior, r0, p0_vec, y_vec, m_vec, R_NilValue, K_mat
  );

  arma::vec alpha0 = as<arma::vec>(prior_par["alpha0"]);
  arma::vec beta0  = as<arma::vec>(prior_par["beta0"]);

  double val = loss_fun_rcpp(
    "BKP", loss, K, wrap(y), wrap(m), R_NilValue, wrap(alpha0), wrap(beta0)
  );

  // guard: NaN or Inf -> return large finite value so sort_index won't crash
  if (!std::isfinite(val)) return std::numeric_limits<double>::max();

  return val;
}

static double eval_dkp_loss_from_gamma(
    const arma::vec& gamma,
    const arma::mat& Xnorm,
    const arma::mat& Y,
    const std::string& prior,
    const double r0,
    const arma::vec& p0,
    const std::string& loss,
    const std::string& kernel,
    const bool isotropic
) {
  arma::vec theta = arma::exp(std::log(10.0) * gamma);

  arma::mat K = kernel_matrix_arma(Xnorm, Xnorm, theta, kernel, isotropic, true);
  K.diag().zeros();

  NumericVector p0_vec = wrap(p0);
  NumericMatrix Y_mat = wrap(Y);
  NumericMatrix K_mat = wrap(K);

  List prior_par = get_prior_rcpp(
    "DKP", prior, r0, p0_vec, R_NilValue, R_NilValue, Y_mat, K_mat
  );

  arma::mat alpha0 = as<arma::mat>(prior_par["alpha0"]);

  double val = loss_fun_rcpp(
    "DKP", loss, K, R_NilValue, R_NilValue, wrap(Y),
    R_NilValue, R_NilValue, wrap(alpha0)
  );

  if (!std::isfinite(val)) return std::numeric_limits<double>::max();

  return val;
}

// ---------- NLOPT objective ----------
struct BKPOptData {
  const arma::mat* Xnorm;
  const arma::vec* y;
  const arma::vec* m;
  std::string prior;
  double r0;
  double p0;
  std::string loss;
  std::string kernel;
  bool isotropic;
};

struct DKPOptData {
  const arma::mat* Xnorm;
  const arma::mat* Y;
  std::string prior;
  double r0;
  arma::vec p0;
  std::string loss;
  std::string kernel;
  bool isotropic;
};

static double bkp_nlopt_obj(unsigned n, const double* x, double* grad, void* f_data) {
  if (grad != nullptr) {
    for (unsigned i = 0; i < n; ++i) grad[i] = 0.0; // SBPLX does not use gradient
  }

  BKPOptData* d = reinterpret_cast<BKPOptData*>(f_data);
  arma::vec gamma(n);
  for (unsigned i = 0; i < n; ++i) gamma[i] = x[i];

  return eval_bkp_loss_from_gamma(
    gamma, *(d->Xnorm), *(d->y), *(d->m),
    d->prior, d->r0, d->p0, d->loss, d->kernel, d->isotropic
  );
}

static double dkp_nlopt_obj(unsigned n, const double* x, double* grad, void* f_data) {
  if (grad != nullptr) {
    for (unsigned i = 0; i < n; ++i) grad[i] = 0.0;
  }

  DKPOptData* d = reinterpret_cast<DKPOptData*>(f_data);
  arma::vec gamma(n);
  for (unsigned i = 0; i < n; ++i) gamma[i] = x[i];

  return eval_dkp_loss_from_gamma(
    gamma, *(d->Xnorm), *(d->Y),
    d->prior, d->r0, d->p0, d->loss, d->kernel, d->isotropic
  );
}

static Rcpp::List nloptr_refine(
    arma::vec gamma_init,
    const arma::mat& Xnorm,
    const arma::vec& y,
    const arma::vec& m,
    const std::string& prior,
    const double r0,
    const double p0,
    const std::string& loss,
    const std::string& kernel,
    const bool isotropic,
    const arma::vec& lower,
    const arma::vec& upper,
    const int max_eval
) {
  for (arma::uword i = 0; i < gamma_init.n_elem; ++i) {
    gamma_init[i] = std::max(lower[i], std::min(upper[i], gamma_init[i]));
  }
  std::vector<double> x(gamma_init.begin(), gamma_init.end());
  std::vector<double> lb(lower.begin(), lower.end());
  std::vector<double> ub(upper.begin(), upper.end());

  BKPOptData data{&Xnorm, &y, &m, prior, r0, p0, loss, kernel, isotropic};

  nlopt_opt opt = nlopt_create(NLOPT_LN_SBPLX, static_cast<unsigned>(x.size()));
  nlopt_set_lower_bounds(opt, lb.data());
  nlopt_set_upper_bounds(opt, ub.data());
  nlopt_set_min_objective(opt, bkp_nlopt_obj, &data);
  nlopt_set_maxeval(opt, max_eval);
  nlopt_set_xtol_rel(opt, 1e-6);

  double f_min = std::numeric_limits<double>::infinity();
  nlopt_result rc = nlopt_optimize(opt, x.data(), &f_min);
  nlopt_destroy(opt);

  arma::vec g_opt(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) g_opt[i] = x[i];

  if (!std::isfinite(f_min)) {
    f_min = eval_bkp_loss_from_gamma(g_opt, Xnorm, y, m, prior, r0, p0, loss, kernel, isotropic);
  }

  return List::create(
    Named("gamma") = g_opt,
    Named("value") = f_min,
    Named("status") = static_cast<int>(rc)
  );
}

static Rcpp::List nloptr_refine_dkp(
    arma::vec gamma_init,
    const arma::mat& Xnorm,
    const arma::mat& Y,
    const std::string& prior,
    const double r0,
    const arma::vec& p0,
    const std::string& loss,
    const std::string& kernel,
    const bool isotropic,
    const arma::vec& lower,
    const arma::vec& upper,
    const int max_eval
) {
  for (arma::uword i = 0; i < gamma_init.n_elem; ++i) {
    gamma_init[i] = std::max(lower[i], std::min(upper[i], gamma_init[i]));
  }
  std::vector<double> x(gamma_init.begin(), gamma_init.end());
  std::vector<double> lb(lower.begin(), lower.end());
  std::vector<double> ub(upper.begin(), upper.end());

  DKPOptData data{&Xnorm, &Y, prior, r0, p0, loss, kernel, isotropic};

  nlopt_opt opt = nlopt_create(NLOPT_LN_SBPLX, static_cast<unsigned>(x.size()));
  nlopt_set_lower_bounds(opt, lb.data());
  nlopt_set_upper_bounds(opt, ub.data());
  nlopt_set_min_objective(opt, dkp_nlopt_obj, &data);
  nlopt_set_maxeval(opt, max_eval);
  nlopt_set_xtol_rel(opt, 1e-6);

  double f_min = std::numeric_limits<double>::infinity();
  nlopt_result rc = nlopt_optimize(opt, x.data(), &f_min);
  nlopt_destroy(opt);

  arma::vec g_opt(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) g_opt[i] = x[i];

  if (!std::isfinite(f_min)) {
    f_min = eval_dkp_loss_from_gamma(g_opt, Xnorm, Y, prior, r0, p0, loss, kernel, isotropic);
  }

  return List::create(
    Named("gamma") = g_opt,
    Named("value") = f_min,
    Named("status") = static_cast<int>(rc)
  );
}

// [[Rcpp::export]]
Rcpp::List optimize_bkp_theta_rcpp(
    const arma::mat& Xnorm,
    const arma::vec& y,
    const arma::vec& m,
    const std::string& prior,
    const double r0,
    const double p0,
    const std::string& loss,
    const std::string& kernel,
    const bool isotropic,
    const arma::mat& init_gamma,
    const arma::vec& lower,
    const arma::vec& upper,
    const int max_iter
) {
  const int n_starts = static_cast<int>(init_gamma.n_rows);

  arma::vec best_gamma = init_gamma.row(0).t();
  double best_val = std::numeric_limits<double>::infinity();

  for (int k = 0; k < n_starts; ++k) {
    arma::vec g0 = init_gamma.row(k).t();

    List ref = nloptr_refine(
      g0, Xnorm, y, m, prior, r0, p0, loss, kernel, isotropic,
      lower, upper, max_iter
    );

    const arma::vec gk = as<arma::vec>(ref["gamma"]);
    const double vk = as<double>(ref["value"]);

    if (vk < best_val) {
      best_val = vk;
      best_gamma = gk;
    }
  }

  arma::vec theta_opt = arma::exp(std::log(10.0) * best_gamma);

  return List::create(
    Named("theta_opt") = theta_opt,
    Named("gamma_opt") = best_gamma,
    Named("loss_min") = best_val
  );
}

// [[Rcpp::export]]
Rcpp::List optimize_dkp_theta_rcpp(
    const arma::mat& Xnorm,
    const arma::mat& Y,
    const std::string& prior,
    const double r0,
    const arma::vec& p0,
    const std::string& loss,
    const std::string& kernel,
    const bool isotropic,
    const arma::mat& init_gamma,
    const arma::vec& lower,
    const arma::vec& upper,
    const int max_iter
) {
  const int n_starts = static_cast<int>(init_gamma.n_rows);

  arma::vec best_gamma = init_gamma.row(0).t();
  double best_val = std::numeric_limits<double>::infinity();

  for (int k = 0; k < n_starts; ++k) {
    arma::vec g0 = init_gamma.row(k).t();

    List ref = nloptr_refine_dkp(
      g0, Xnorm, Y, prior, r0, p0, loss, kernel, isotropic,
      lower, upper, max_iter
    );

    const arma::vec gk = as<arma::vec>(ref["gamma"]);
    const double vk = as<double>(ref["value"]);

    if (vk < best_val) {
      best_val = vk;
      best_gamma = gk;
    }
  }

  arma::vec theta_opt = arma::exp(std::log(10.0) * best_gamma);

  return List::create(
    Named("theta_opt") = theta_opt,
    Named("gamma_opt") = best_gamma,
    Named("loss_min") = best_val
  );
}
