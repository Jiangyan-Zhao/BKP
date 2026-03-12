#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>
#include <cmath>
#include <random>

using namespace Rcpp;

// squared Euclidean distance
inline double dist2_row(const NumericMatrix& X, int i, int j) {
  int d = X.ncol();
  double s = 0.0;
  for (int k = 0; k < d; ++k) {
    double diff = X(i, k) - X(j, k);
    s += diff * diff;
  }
  return s;
}

inline double dist2_point_row(const NumericMatrix& X, const NumericVector& x, int j) {
  int d = X.ncol();
  double s = 0.0;
  for (int k = 0; k < d; ++k) {
    double diff = x[k] - X(j, k);
    s += diff * diff;
  }
  return s;
}

// k nearest indices from candidate set to row "anchor"
std::vector<int> k_nearest_from_set(
    const NumericMatrix& X,
    int anchor,
    const std::vector<int>& candidates,
    int k
) {
  std::vector< std::pair<double,int> > dd;
  dd.reserve(candidates.size());
  for (int idx : candidates) {
    dd.push_back({dist2_row(X, anchor, idx), idx});
  }
  std::sort(dd.begin(), dd.end(),
            [](const std::pair<double,int>& a, const std::pair<double,int>& b){ return a.first < b.first; });

  k = std::min(k, (int)dd.size());
  std::vector<int> out;
  out.reserve(k);
  for (int i = 0; i < k; ++i) out.push_back(dd[i].second);
  return out;
}

// nearest remaining point to row "anchor"
int nearest_from_set(const NumericMatrix& X, int anchor, const std::vector<int>& candidates) {
  double best = std::numeric_limits<double>::infinity();
  int best_idx = candidates[0];
  for (int idx : candidates) {
    double d = dist2_row(X, anchor, idx);
    if (d < best) {
      best = d;
      best_idx = idx;
    }
  }
  return best_idx;
}

double min_pairwise_dist2(const NumericMatrix& X, const std::vector<int>& idx) {
  if (idx.size() < 2) return 0.0;
  double best = std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < idx.size(); ++i) {
    for (size_t j = i + 1; j < idx.size(); ++j) {
      double d = dist2_row(X, idx[i], idx[j]);
      if (d < best) best = d;
    }
  }
  return best;
}

std::vector<int> setdiff_idx(int n, const std::vector<int>& removed) {
  std::vector<char> mark(n, 0);
  for (int x : removed) mark[x] = 1;
  std::vector<int> out;
  out.reserve(n - removed.size());
  for (int i = 0; i < n; ++i) if (!mark[i]) out.push_back(i);
  return out;
}

// [[Rcpp::export]]
List get_twin_indices_rcpp(
    NumericMatrix data,
    int g,
    Nullable<int> v = R_NilValue,
    int runs = 10,
    int seed = 123
) {
  RNGScope scope;

  std::mt19937 rng(seed);

  int n = data.nrow();
  if (g <= 0) stop("'g' must be positive.");
  g = std::min(g, n);

  int v_eff = v.isNotNull() ? as<int>(v) : 2 * g;
  v_eff = std::max(1, std::min(v_eff, std::max(1, n - g)));

  int r = std::max(1, (int)std::ceil((double)n / (double)g));

  std::vector< std::vector<int> > all_g(runs);
  std::vector<double> run_score(runs, -1.0);

  Function set_seed("set.seed");
  set_seed(seed);

  for (int run = 0; run < runs; ++run) {
    std::vector<int> remaining(n);
    std::iota(remaining.begin(), remaining.end(), 0);

    int start = (int)std::floor(R::runif(0.0, (double)n));
    int pos = std::min(std::max(start, 0), n - 1);

    std::vector<int> chosen;
    chosen.reserve(g);

    while (!remaining.empty() && (int)chosen.size() < g) {
      int k_now = std::min(r, (int)remaining.size());
      std::vector<int> nn = k_nearest_from_set(data, pos, remaining, k_now);

      // add nearest as representative
      chosen.push_back(nn[0]);

      // remove cluster nn from remaining
      std::vector<char> rm(n, 0);
      for (int x : nn) rm[x] = 1;
      std::vector<int> new_remaining;
      new_remaining.reserve(remaining.size());
      for (int x : remaining) if (!rm[x]) new_remaining.push_back(x);
      remaining.swap(new_remaining);

      if (remaining.empty() || (int)chosen.size() >= g) break;

      // move to nearest remaining point from the farthest in nn
      int far_idx = nn.back();
      pos = nearest_from_set(data, far_idx, remaining);
    }

    // if not enough points, random fill
    if ((int)chosen.size() < g) {
      std::vector<int> pool = setdiff_idx(n, chosen);
      std::shuffle(pool.begin(), pool.end(), rng);
      int need = g - (int)chosen.size();
      for (int i = 0; i < need && i < (int)pool.size(); ++i) chosen.push_back(pool[i]);
    }

    // keep exactly g unique
    std::sort(chosen.begin(), chosen.end());
    chosen.erase(std::unique(chosen.begin(), chosen.end()), chosen.end());
    if ((int)chosen.size() > g) chosen.resize(g);

    // pad if unique shrink
    if ((int)chosen.size() < g) {
      std::vector<int> pool = setdiff_idx(n, chosen);
      int need = g - (int)chosen.size();
      for (int i = 0; i < need && i < (int)pool.size(); ++i) chosen.push_back(pool[i]);
    }

    run_score[run] = min_pairwise_dist2(data, chosen);
    all_g[run] = chosen;
  }

  int best_run = std::distance(run_score.begin(), std::max_element(run_score.begin(), run_score.end()));
  std::vector<int> g_idx0 = all_g[best_run]; // 0-based

  // theta_l: nearest global distance for each point, take 99th percentile
  std::vector<double> mind(n, std::numeric_limits<double>::infinity());
  for (int i = 0; i < n; ++i) {
    for (int j : g_idx0) {
      double d = dist2_row(data, i, j);
      if (d < mind[i]) mind[i] = d;
    }
  }
  std::sort(mind.begin(), mind.end());
  int q01 = (int)std::floor((n - (int)g_idx0.size()) * 0.01);
  int pos_q = std::max(0, n - q01 - 1);
  double theta_l = std::sqrt(mind[pos_q]);

  // vIndices: from non-global, random sample (simple/robust)
  std::vector<int> non_g = setdiff_idx(n, g_idx0);
  std::shuffle(non_g.begin(), non_g.end(), rng);
  if ((int)non_g.size() > v_eff) non_g.resize(v_eff);

  // convert to 1-based for R
  IntegerVector gIndices(g_idx0.size());
  for (size_t i = 0; i < g_idx0.size(); ++i) gIndices[i] = g_idx0[i] + 1;

  IntegerVector vIndices(non_g.size());
  for (size_t i = 0; i < non_g.size(); ++i) vIndices[i] = non_g[i] + 1;

  return List::create(
    _["gIndices"] = gIndices,
    _["theta_l"]  = theta_l,
    _["vIndices"] = vIndices
  );
}