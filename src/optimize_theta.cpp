// -----------------------------------------------------------------------------
// This file implements kernel hyperparameter optimization for BKP and DKP.
//
// This file contains the shared optimization workflow for both Beta Kernel
// Process and Dirichlet Kernel Process models. The BKP and DKP routines are
// kept together because they use the same nloptr-based refinement strategy.
//
// The objective functions use pure Armadillo engines only, so the outer
// multi-start loop can be safely parallelized with OpenMP.
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
#include <nloptrAPI.h>
#include <limits>
#include <cmath>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

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

void get_prior_bkp_arma(
    const std::string& prior,
    const double r0,
    const double p0,
    const arma::vec& y,
    const arma::vec& m,
    const arma::mat& K,
    arma::vec& alpha0,
    arma::vec& beta0
);

arma::mat get_prior_dkp_arma(
    const std::string& prior,
    const double r0,
    const arma::vec& p0,
    const arma::mat& Y,
    const arma::mat& K
);

double loss_bkp_arma(
    const std::string& loss,
    const arma::mat& K,
    const arma::vec& y,
    const arma::vec& m,
    const arma::vec& alpha0,
    const arma::vec& beta0,
    const arma::vec& data_scale
);

double loss_dkp_arma(
    const std::string& loss,
    const arma::mat& K,
    const arma::mat& Y,
    const arma::mat& alpha0,
    const arma::vec& data_scale
);

// -------- END: declarations from other cpp files --------


// -----------------------------------------------------------------------------
// Pure C++ objective evaluation from gamma = log10(theta)
// -----------------------------------------------------------------------------

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
    const bool isotropic,
    const std::string& ess,
    const arma::vec& m_shepard_loo
) {
  arma::vec theta = arma::exp(std::log(10.0) * gamma);

  arma::mat K = kernel_matrix_arma(Xnorm, Xnorm, theta, kernel, isotropic, true);
  K.diag().zeros();

  arma::vec alpha0;
  arma::vec beta0;

  get_prior_bkp_arma(prior, r0, p0, y, m, K, alpha0, beta0);

  arma::vec data_scale = arma::ones<arma::vec>(K.n_rows);
  if (ess == "shepard") {
    arma::vec m_kernel = K * m;
    arma::vec rho = arma::max(K, 1);
    arma::vec m_target = rho % m_shepard_loo;
    arma::uvec positive_kernel_mass = arma::find(m_kernel > 0.0);
    data_scale.zeros();
    data_scale.elem(positive_kernel_mass) =
      m_target.elem(positive_kernel_mass) / m_kernel.elem(positive_kernel_mass);
  }

  double val = loss_bkp_arma(loss, K, y, m, alpha0, beta0, data_scale);

  // Guard: NaN or Inf -> return a large finite value.
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
    const bool isotropic,
    const std::string& ess,
    const arma::vec& m_shepard_loo
) {
  arma::vec theta = arma::exp(std::log(10.0) * gamma);

  arma::mat K = kernel_matrix_arma(Xnorm, Xnorm, theta, kernel, isotropic, true);
  K.diag().zeros();

  arma::mat alpha0 = get_prior_dkp_arma(prior, r0, p0, Y, K);

  arma::vec data_scale = arma::ones<arma::vec>(K.n_rows);
  if (ess == "shepard") {
    arma::vec m = arma::sum(Y, 1);
    arma::vec m_kernel = K * m;
    arma::vec rho = arma::max(K, 1);
    arma::vec m_target = rho % m_shepard_loo;
    arma::uvec positive_kernel_mass = arma::find(m_kernel > 0.0);
    data_scale.zeros();
    data_scale.elem(positive_kernel_mass) =
      m_target.elem(positive_kernel_mass) / m_kernel.elem(positive_kernel_mass);
  }

  double val = loss_dkp_arma(loss, K, Y, alpha0, data_scale);

  if (!std::isfinite(val)) return std::numeric_limits<double>::max();

  return val;
}


// -----------------------------------------------------------------------------
// NLOPT objective data
// -----------------------------------------------------------------------------

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
  std::string ess;
  const arma::vec* m_shepard_loo;
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
  std::string ess;
  const arma::vec* m_shepard_loo;
};


struct OptResult {
  arma::vec gamma;
  double value;
  int status;

  OptResult()
    : gamma(), value(std::numeric_limits<double>::infinity()), status(0) {}
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
    d->prior, d->r0, d->p0, d->loss, d->kernel, d->isotropic,
    d->ess, *(d->m_shepard_loo)
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
    d->prior, d->r0, d->p0, d->loss, d->kernel, d->isotropic,
    d->ess, *(d->m_shepard_loo)
  );
}


// -----------------------------------------------------------------------------
// One local nloptr refinement from a single start
// -----------------------------------------------------------------------------

static OptResult nloptr_refine(
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
    const std::string& ess,
    const arma::vec& m_shepard_loo,
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

  BKPOptData data{&Xnorm, &y, &m, prior, r0, p0, loss, kernel, isotropic, ess, &m_shepard_loo};

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
    f_min = eval_bkp_loss_from_gamma(
      g_opt, Xnorm, y, m, prior, r0, p0, loss, kernel, isotropic,
      ess, m_shepard_loo
    );
  }

  OptResult out;
  out.gamma = g_opt;
  out.value = f_min;
  out.status = static_cast<int>(rc);

  return out;
}


static OptResult nloptr_refine_dkp(
    arma::vec gamma_init,
    const arma::mat& Xnorm,
    const arma::mat& Y,
    const std::string& prior,
    const double r0,
    const arma::vec& p0,
    const std::string& loss,
    const std::string& kernel,
    const bool isotropic,
    const std::string& ess,
    const arma::vec& m_shepard_loo,
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

  DKPOptData data{&Xnorm, &Y, prior, r0, p0, loss, kernel, isotropic, ess, &m_shepard_loo};

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
    f_min = eval_dkp_loss_from_gamma(
      g_opt, Xnorm, Y, prior, r0, p0, loss, kernel, isotropic, ess, m_shepard_loo
    );
  }

  OptResult out;
  out.gamma = g_opt;
  out.value = f_min;
  out.status = static_cast<int>(rc);

  return out;
}


// -----------------------------------------------------------------------------
// Exported multi-start optimization wrappers
// -----------------------------------------------------------------------------

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
    const int max_iter,
    const int n_threads = 1,
    const std::string& ess = "none",
    Nullable<NumericVector> m_shepard_loo = R_NilValue
) {
  const int n_starts = static_cast<int>(init_gamma.n_rows);

  arma::vec m_shepard_loo_vec;
  if (ess == "shepard") {
    if (m_shepard_loo.isNull()) {
      stop("'m_shepard_loo' must be provided when ess = 'shepard'.");
    }
    NumericVector m_shepard_loo_R(m_shepard_loo);
    m_shepard_loo_vec = as<arma::vec>(m_shepard_loo_R);
    if (m_shepard_loo_vec.n_elem != m.n_elem) {
      stop("'m_shepard_loo' must have the same length as 'm'.");
    }
  } else {
    m_shepard_loo_vec = arma::ones<arma::vec>(m.n_elem);
  }

  int n_threads_used = 1;

#ifdef _OPENMP
  n_threads_used = std::max(1, std::min(n_threads, n_starts));
#else
  n_threads_used = 1;
#endif

  std::vector<OptResult> results(n_starts);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads_used)
#endif
  for (int k = 0; k < n_starts; ++k) {
    try {
      arma::vec g0 = init_gamma.row(k).t();

      results[k] = nloptr_refine(
        g0, Xnorm, y, m,
        prior, r0, p0, loss, kernel, isotropic, ess, m_shepard_loo_vec,
        lower, upper, max_iter
      );

    } catch (...) {
      results[k].gamma = init_gamma.row(k).t();
      results[k].value = std::numeric_limits<double>::infinity();
      results[k].status = -999;
    }
  }

  arma::vec best_gamma = init_gamma.row(0).t();
  double best_val = std::numeric_limits<double>::infinity();

  arma::vec all_loss(n_starts);
  arma::ivec all_status(n_starts);

  for (int k = 0; k < n_starts; ++k) {
    all_loss[k] = results[k].value;
    all_status[k] = results[k].status;

    if (results[k].value < best_val) {
      best_val = results[k].value;
      best_gamma = results[k].gamma;
    }
  }

  arma::vec theta_opt = arma::exp(std::log(10.0) * best_gamma);

  return List::create(
    Named("theta_opt") = theta_opt,
    Named("gamma_opt") = best_gamma,
    Named("loss_min") = best_val,
    Named("all_loss") = all_loss,
    Named("all_status") = all_status,
    Named("n_threads") = n_threads_used
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
    const int max_iter,
    const int n_threads = 1,
    const std::string& ess = "none",
    Nullable<NumericVector> m_shepard_loo = R_NilValue
) {
  const int n_starts = static_cast<int>(init_gamma.n_rows);

  arma::vec m_shepard_loo_vec;
  if (ess == "shepard") {
    if (m_shepard_loo.isNull()) {
      stop("'m_shepard_loo' must be provided when ess = 'shepard'.");
    }
    NumericVector m_shepard_loo_R(m_shepard_loo);
    m_shepard_loo_vec = as<arma::vec>(m_shepard_loo_R);
    if (m_shepard_loo_vec.n_elem != Y.n_rows) {
      stop("'m_shepard_loo' must have the same length as the number of rows in 'Y'.");
    }
  } else {
    m_shepard_loo_vec = arma::ones<arma::vec>(Y.n_rows);
  }

  int n_threads_used = 1;

#ifdef _OPENMP
  n_threads_used = std::max(1, std::min(n_threads, n_starts));
#else
  n_threads_used = 1;
#endif

  std::vector<OptResult> results(n_starts);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads_used)
#endif
  for (int k = 0; k < n_starts; ++k) {
    try {
      arma::vec g0 = init_gamma.row(k).t();

      results[k] = nloptr_refine_dkp(
        g0, Xnorm, Y,
        prior, r0, p0, loss, kernel, isotropic, ess, m_shepard_loo_vec,
        lower, upper, max_iter
      );

    } catch (...) {
      results[k].gamma = init_gamma.row(k).t();
      results[k].value = std::numeric_limits<double>::infinity();
      results[k].status = -999;
    }
  }

  arma::vec best_gamma = init_gamma.row(0).t();
  double best_val = std::numeric_limits<double>::infinity();

  arma::vec all_loss(n_starts);
  arma::ivec all_status(n_starts);

  for (int k = 0; k < n_starts; ++k) {
    all_loss[k] = results[k].value;
    all_status[k] = results[k].status;

    if (results[k].value < best_val) {
      best_val = results[k].value;
      best_gamma = results[k].gamma;
    }
  }

  arma::vec theta_opt = arma::exp(std::log(10.0) * best_gamma);

  return List::create(
    Named("theta_opt") = theta_opt,
    Named("gamma_opt") = best_gamma,
    Named("loss_min") = best_val,
    Named("all_loss") = all_loss,
    Named("all_status") = all_status,
    Named("n_threads") = n_threads_used
  );
}
