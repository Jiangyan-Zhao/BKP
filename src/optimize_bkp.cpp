#include <RcppArmadillo.h>
#include <nloptrAPI.h>
#include <limits>
#include <cmath>
#include <algorithm>
#include <vector>

// [[Rcpp::depends(RcppArmadillo, nloptr)]]

using namespace Rcpp;

// ---- declarations from other cpp files ----
arma::mat kernel_matrix_rcpp(SEXP X, SEXP Xprime, NumericVector theta, std::string kernel, bool isotropic);

List get_prior_bkp_noninformative_rcpp(int n);
List get_prior_bkp_fixed_rcpp(int n, double r0, double p0);
List get_prior_bkp_adaptive_rcpp(const arma::mat& K, const arma::vec& y, const arma::vec& m, double r0);

double loss_fun_brier_bkp_rcpp(
  const arma::mat& K, const arma::vec& y, const arma::vec& m,
  const arma::vec& alpha0, const arma::vec& beta0
);

double loss_fun_logloss_bkp_rcpp(
  const arma::mat& K, const arma::vec& y, const arma::vec& m,
  const arma::vec& alpha0, const arma::vec& beta0
);

arma::mat get_prior_dkp_noninformative_rcpp(int n, int q);
arma::mat get_prior_dkp_fixed_rcpp(int n, double r0, const arma::vec& p0);
arma::mat get_prior_dkp_adaptive_rcpp(const arma::mat& K, const arma::mat& Y, double r0);

double loss_fun_brier_dkp_rcpp(
  const arma::mat& K, const arma::mat& Y, const arma::mat& alpha0
);

double loss_fun_logloss_dkp_rcpp(
  const arma::mat& K, const arma::mat& Y, const arma::mat& alpha0
);

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

  NumericVector theta_nv(theta.begin(), theta.end());
  arma::mat K = kernel_matrix_rcpp(wrap(Xnorm), R_NilValue, theta_nv, kernel, isotropic);
  K.diag().zeros();

  List prior_par;
  if (prior == "noninformative") {
    prior_par = get_prior_bkp_noninformative_rcpp(Xnorm.n_rows);
  } else if (prior == "fixed") {
    prior_par = get_prior_bkp_fixed_rcpp(Xnorm.n_rows, r0, p0);
  } else if (prior == "adaptive") {
    prior_par = get_prior_bkp_adaptive_rcpp(K, y, m, r0);
  } else {
    stop("Unsupported prior: " + prior);
  }

  arma::vec alpha0 = as<arma::vec>(prior_par["alpha0"]);
  arma::vec beta0  = as<arma::vec>(prior_par["beta0"]);

  double val = std::numeric_limits<double>::infinity();

  if (loss == "brier") {
    val = loss_fun_brier_bkp_rcpp(K, y, m, alpha0, beta0);
  } else if (loss == "log_loss") {
    val = loss_fun_logloss_bkp_rcpp(K, y, m, alpha0, beta0);
  } else {
    stop("Unsupported loss: " + loss);
  }

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

  NumericVector theta_nv(theta.begin(), theta.end());
  arma::mat K = kernel_matrix_rcpp(wrap(Xnorm), R_NilValue, theta_nv, kernel, isotropic);
  K.diag().zeros();

  arma::mat alpha0;
  if (prior == "noninformative") {
    alpha0 = get_prior_dkp_noninformative_rcpp(Xnorm.n_rows, Y.n_cols);
  } else if (prior == "fixed") {
    alpha0 = get_prior_dkp_fixed_rcpp(Xnorm.n_rows, r0, p0);
  } else if (prior == "adaptive") {
    alpha0 = get_prior_dkp_adaptive_rcpp(K, Y, r0);
  } else {
    stop("Unsupported prior: " + prior);
  }

  double val = std::numeric_limits<double>::infinity();

  if (loss == "brier") {
    val = loss_fun_brier_dkp_rcpp(K, Y, alpha0);
  } else if (loss == "log_loss") {
    val = loss_fun_logloss_dkp_rcpp(K, Y, alpha0);
  } else {
    stop("Unsupported loss: " + loss);
  }

  if (!std::isfinite(val)) return std::numeric_limits<double>::max();

  return val;
}

static List generate_anisotropic_candidates(
    const int p,
    const int n_candidates,
    const double lower,
    const double upper
) {
  arma::mat cand(n_candidates, p, arma::fill::zeros);

  // first candidate = zero vector (theta = 1 for all dims)
  cand.row(0).zeros();

  if (n_candidates > 1) {
    // Latin Hypercube Sampling for remaining n_candidates-1 points
    const int n_lhs = n_candidates - 1;
    arma::mat lhs(n_lhs, p);

    for (int j = 0; j < p; ++j) {
      arma::uvec perm = arma::randperm(n_lhs);
      for (int i = 0; i < n_lhs; ++i) {
        double u = (static_cast<double>(perm[i]) + R::runif(0.0, 1.0)) / static_cast<double>(n_lhs);
        lhs(i, j) = lower + u * (upper - lower);
      }
    }

    cand.rows(1, n_candidates - 1) = lhs;
  }

  return List::create(Named("cand") = cand);
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
    const int n_grid,
    const int n_starts,
    const int max_iter,
    const double g_lower,
    const double g_upper
) {
  if (n_grid < 5) stop("n_grid must be >= 5.");
  if (n_starts < 1) stop("n_starts must be >= 1.");
  if (max_iter < 1) stop("max_iter must be >= 1.");
  if (Xnorm.n_rows == 0 || Xnorm.n_cols == 0) stop("Xnorm must be non-empty.");
  if (static_cast<int>(y.n_elem) != static_cast<int>(Xnorm.n_rows)) stop("length(y) must equal nrow(Xnorm).");
  if (static_cast<int>(m.n_elem) != static_cast<int>(Xnorm.n_rows)) stop("length(m) must equal nrow(Xnorm).");

  arma::vec gamma_opt;
  double loss_min = std::numeric_limits<double>::infinity();

  if (isotropic) {
    // ---- coarse grid in 1D ----
    arma::vec grid = arma::linspace(g_lower, g_upper, n_grid);
    arma::vec vals(n_grid, arma::fill::value(std::numeric_limits<double>::max()));

    for (int i = 0; i < n_grid; ++i) {
      arma::vec g(1);
      g[0] = grid[i];
      vals[i] = eval_bkp_loss_from_gamma(g, Xnorm, y, m, prior, r0, p0, loss, kernel, true);
    }

    arma::uvec ord = arma::sort_index(vals, "ascend");
    const int k_starts = std::min(n_starts, n_grid);

    arma::vec best_gamma(1);
    best_gamma[0] = grid[ord[0]];
    double best_val = std::numeric_limits<double>::infinity();

    for (int k = 0; k < k_starts; ++k) {
      const int idx = static_cast<int>(ord[k]);
      const int i_left  = std::max(0, idx - 1);
      const int i_right = std::min(n_grid - 1, idx + 1);

      arma::vec g0(1);
      g0[0] = grid[idx];

      arma::vec lower(1), upper(1);
      lower[0] = grid[i_left];
      upper[0] = grid[i_right];

      List ref = nloptr_refine(
        g0, Xnorm, y, m, prior, r0, p0, loss, kernel, true,
        lower, upper, max_iter
      );

      const arma::vec gk = as<arma::vec>(ref["gamma"]);
      const double vk = as<double>(ref["value"]);

      if (vk < best_val) {
        best_val = vk;
        best_gamma = gk;
      }
    }

    gamma_opt = best_gamma;
    loss_min = best_val;

  } else {
    // ---- anisotropic: random coarse candidates + multi-start local refine ----
    RNGScope scope;

    const int p = static_cast<int>(Xnorm.n_cols);
    const int n_candidates = std::max(10, n_grid);

    List cand_obj = generate_anisotropic_candidates(p, n_candidates, g_lower, g_upper);
    arma::mat cand = as<arma::mat>(cand_obj["cand"]);
    arma::vec vals(n_candidates, arma::fill::value(std::numeric_limits<double>::max()));

    for (int i = 0; i < n_candidates; ++i) {
      arma::vec g = cand.row(i).t();
      vals[i] = eval_bkp_loss_from_gamma(g, Xnorm, y, m, prior, r0, p0, loss, kernel, false);
    }

    arma::uvec ord = arma::sort_index(vals, "ascend");
    const int k_starts = std::min(n_starts, n_candidates);

    arma::vec lower(p), upper(p);
    lower.fill(g_lower);
    upper.fill(g_upper);

    arma::vec best_gamma = cand.row(ord[0]).t();
    double best_val = std::numeric_limits<double>::infinity();

    for (int k = 0; k < k_starts; ++k) {
      const int idx = static_cast<int>(ord[k]);
      arma::vec g0 = cand.row(idx).t();

      List ref = nloptr_refine(
        g0, Xnorm, y, m, prior, r0, p0, loss, kernel, false,
        lower, upper, max_iter
      );

      const arma::vec gk = as<arma::vec>(ref["gamma"]);
      const double vk = as<double>(ref["value"]);

      if (vk < best_val) {
        best_val = vk;
        best_gamma = gk;
      }
    }

    gamma_opt = best_gamma;
    loss_min = best_val;
  }

  arma::vec theta_opt = arma::exp(std::log(10.0) * gamma_opt);

  return List::create(
    Named("theta_opt") = theta_opt,
    Named("gamma_opt") = gamma_opt,
    Named("loss_min") = loss_min
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
    const int n_grid,
    const int n_starts,
    const int max_iter,
    const double g_lower,
    const double g_upper
) {
  if (n_grid < 5) stop("n_grid must be >= 5.");
  if (n_starts < 1) stop("n_starts must be >= 1.");
  if (max_iter < 1) stop("max_iter must be >= 1.");
  if (Xnorm.n_rows == 0 || Xnorm.n_cols == 0) stop("Xnorm must be non-empty.");
  if (Y.n_rows != Xnorm.n_rows) stop("nrow(Y) must equal nrow(Xnorm).");
  if (prior == "fixed" && static_cast<int>(p0.n_elem) != static_cast<int>(Y.n_cols)) {
    stop("For fixed prior, length(p0) must equal ncol(Y).");
  }

  arma::vec gamma_opt;
  double loss_min = std::numeric_limits<double>::infinity();

  if (isotropic) {
    arma::vec grid = arma::linspace(g_lower, g_upper, n_grid);
    arma::vec vals(n_grid, arma::fill::value(std::numeric_limits<double>::max()));

    for (int i = 0; i < n_grid; ++i) {
      arma::vec g(1);
      g[0] = grid[i];
      vals[i] = eval_dkp_loss_from_gamma(g, Xnorm, Y, prior, r0, p0, loss, kernel, true);
    }

    arma::uvec ord = arma::sort_index(vals, "ascend");
    const int k_starts = std::min(n_starts, n_grid);

    arma::vec best_gamma(1);
    best_gamma[0] = grid[ord[0]];
    double best_val = std::numeric_limits<double>::infinity();

    for (int k = 0; k < k_starts; ++k) {
      const int idx = static_cast<int>(ord[k]);
      const int i_left = std::max(0, idx - 1);
      const int i_right = std::min(n_grid - 1, idx + 1);

      arma::vec g0(1);
      g0[0] = grid[idx];

      arma::vec lower(1), upper(1);
      lower[0] = grid[i_left];
      upper[0] = grid[i_right];

      List ref = nloptr_refine_dkp(
        g0, Xnorm, Y, prior, r0, p0, loss, kernel, true,
        lower, upper, max_iter
      );

      const arma::vec gk = as<arma::vec>(ref["gamma"]);
      const double vk = as<double>(ref["value"]);

      if (vk < best_val) {
        best_val = vk;
        best_gamma = gk;
      }
    }

    gamma_opt = best_gamma;
    loss_min = best_val;
  } else {
    RNGScope scope;

    const int p = static_cast<int>(Xnorm.n_cols);
    const int n_candidates = std::max(10, n_grid);

    List cand_obj = generate_anisotropic_candidates(p, n_candidates, g_lower, g_upper);
    arma::mat cand = as<arma::mat>(cand_obj["cand"]);
    arma::vec vals(n_candidates, arma::fill::value(std::numeric_limits<double>::max()));

    for (int i = 0; i < n_candidates; ++i) {
      arma::vec g = cand.row(i).t();
      vals[i] = eval_dkp_loss_from_gamma(g, Xnorm, Y, prior, r0, p0, loss, kernel, false);
    }

    arma::uvec ord = arma::sort_index(vals, "ascend");
    const int k_starts = std::min(n_starts, n_candidates);

    arma::vec lower(p), upper(p);
    lower.fill(g_lower);
    upper.fill(g_upper);

    arma::vec best_gamma = cand.row(ord[0]).t();
    double best_val = std::numeric_limits<double>::infinity();

    for (int k = 0; k < k_starts; ++k) {
      const int idx = static_cast<int>(ord[k]);
      arma::vec g0 = cand.row(idx).t();

      List ref = nloptr_refine_dkp(
        g0, Xnorm, Y, prior, r0, p0, loss, kernel, false,
        lower, upper, max_iter
      );

      const arma::vec gk = as<arma::vec>(ref["gamma"]);
      const double vk = as<double>(ref["value"]);

      if (vk < best_val) {
        best_val = vk;
        best_gamma = gk;
      }
    }

    gamma_opt = best_gamma;
    loss_min = best_val;
  }

  arma::vec theta_opt = arma::exp(std::log(10.0) * gamma_opt);

  return List::create(
    Named("theta_opt") = theta_opt,
    Named("gamma_opt") = gamma_opt,
    Named("loss_min") = loss_min
  );
}