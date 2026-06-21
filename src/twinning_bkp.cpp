/*
 src/twinning_bkp.cpp

 Adapted from twingp/src/twinning.cpp.

 Original copyright:
 Copyright 2023 Akhil Vakayil (akhilv@gatech.edu).
 License: Apache-2.0.

 Modifications for BKP:
 - Separate twin_data from Xnorm.
 - Use twin_data for Twinning selection.
 - Use Xnorm for input-space scoring and theta_l.
 - Remove vIndices, since validation indices are TwinGP-specific.
 - Return R-facing indices as 1-based integers.

 Design principle:
 R validates inputs. C++ performs computation.
 */

// [[Rcpp::plugins("cpp11")]]
#include <memory>
#include <vector>
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <nanoflann.hpp>

using namespace Rcpp;

class DF2
{
private:
  std::shared_ptr<Rcpp::NumericMatrix> df_;
  bool subset_ = false;
  std::vector<std::size_t>* indices_ = nullptr;

public:
  void import_data(Rcpp::NumericMatrix& df)
  {
    df_ = std::make_shared<Rcpp::NumericMatrix>(Rcpp::transpose(df));
  }

  std::size_t kdtree_get_point_count() const
  {
    return subset_ ? indices_->size() : df_->cols();
  }

  double kdtree_get_pt(const std::size_t idx, const std::size_t dim) const
  {
    return subset_ ? (*df_)(dim, indices_->at(idx)) : (*df_)(dim, idx);
  }

  const double* get_row(const std::size_t idx) const
  {
    return &(*df_)(0, idx);
  }

  void subset_on(std::vector<std::size_t>* indices)
  {
    subset_ = true;
    indices_ = indices;
  }

  void subset_off()
  {
    subset_ = false;
    indices_ = nullptr;
  }

  template <class BBOX>
  bool kdtree_get_bbox(BBOX&) const
  {
    return false;
  }
};

typedef nanoflann::KDTreeSingleIndexDynamicAdaptor<
  nanoflann::L2_Adaptor<double, DF2>,
  DF2,
  -1,
  std::size_t
> kdTree;


class BKPTwinKDTree
{
private:
  std::size_t dim_twin_;
  std::size_t dim_x_;
  std::size_t N_;
  std::size_t r_;
  std::size_t runs_;
  std::vector<std::size_t> u1_;
  std::size_t leaf_size_;

  DF2 twin_data_;
  DF2 x_data_;
  Rcpp::List returns_;

public:
  BKPTwinKDTree(
    Rcpp::NumericMatrix& twin_data,
    Rcpp::NumericMatrix& Xnorm,
    std::size_t r,
    std::size_t runs,
    std::vector<std::size_t>& u1,
    std::size_t leaf_size)
    :
    dim_twin_(twin_data.cols()),
    dim_x_(Xnorm.cols()),
    N_(twin_data.rows()),
    r_(r),
    runs_(runs),
    u1_(u1),
    leaf_size_(leaf_size)
  {
    twin_data_.import_data(twin_data);
    x_data_.import_data(Xnorm);
  }

  Rcpp::List twin()
  {
    std::vector<std::vector<std::size_t> > all_indices;
    std::vector<double> min_dist;

    all_indices.resize(runs_);
    min_dist.resize(runs_);

#pragma omp parallel for
    for (std::size_t run = 0; run < runs_; run++)
    {
      kdTree tree(
          dim_twin_,
          twin_data_,
          nanoflann::KDTreeSingleIndexAdaptorParams(leaf_size_)
      );

      nanoflann::KNNResultSet<double> resultSet(r_);
      std::vector<std::size_t> index(r_);
      std::vector<double> distance(r_);

      nanoflann::KNNResultSet<double> resultSet_next_u(1);
      std::size_t index_next_u;
      double distance_next_u;

      std::vector<std::size_t> indices;
      indices.reserve(N_ / r_ + 1);

      std::size_t position = u1_[run];

      while (true)
      {
        resultSet.init(&index[0], &distance[0]);
        tree.findNeighbors(resultSet, twin_data_.get_row(position));

        indices.push_back(index[0]);

        for (std::size_t i = 0; i < r_; i++) {
          tree.removePoint(index[i]);
        }

        resultSet_next_u.init(&index_next_u, &distance_next_u);
        tree.findNeighbors(
          resultSet_next_u,
          twin_data_.get_row(index[r_ - 1])
        );

        position = index_next_u;

        // Same stopping rule as twingp. Do not use nanoflann tree-size APIs.
        if (N_ - indices.size() * r_ <= r_)
        {
          indices.push_back(position);
          break;
        }
      }

      all_indices[run] = indices;
      min_dist[run] = input_space_min_dist(indices);
    }

    int max_index =
      std::max_element(min_dist.begin(), min_dist.end()) - min_dist.begin();

    std::vector<std::size_t> g_indices0 = all_indices[max_index];

    Rcpp::IntegerVector g_indices(g_indices0.size());
    for (std::size_t i = 0; i < g_indices0.size(); i++) {
      g_indices[i] = static_cast<int>(g_indices0[i]) + 1;
    }

    returns_["g_indices"] = g_indices;
    returns_["theta_l"] = theta_l(g_indices0);
    returns_["min_dist"] = min_dist[max_index];
    returns_["best_run"] = max_index + 1;

    return returns_;
  }

  double input_space_min_dist(std::vector<std::size_t> indices)
  {
    if (indices.size() < 2) {
      return 0.0;
    }

    double dist;
    double min_dist_value = std::numeric_limits<double>::max();

    for (std::size_t i = 0; i < indices.size(); i++)
    {
      for (std::size_t j = i + 1; j < indices.size(); j++)
      {
        const double* u = x_data_.get_row(indices[i]);
        const double* v = x_data_.get_row(indices[j]);

        dist = 0.0;

        for (std::size_t k = 0; k < dim_x_; k++) {
          double diff = *(u + k) - *(v + k);
          dist += diff * diff;
        }

        if (dist < min_dist_value) {
          min_dist_value = dist;
        }
      }
    }

    return std::sqrt(min_dist_value);
  }

  double theta_l(std::vector<std::size_t> g_indices0)
  {
    // Build a KD-tree on the selected global points only.
    // Keep subset_on during all neighbour queries. Calling subset_off before
    // querying would make kdtree_get_pt use the first |G| rows instead of G.
    x_data_.subset_on(&g_indices0);

    kdTree tree(
        dim_x_,
        x_data_,
        nanoflann::KDTreeSingleIndexAdaptorParams(leaf_size_)
    );

    std::vector<double> distances;
    distances.resize(N_);

#pragma omp parallel for
    for (std::size_t i = 0; i < N_; i++)
    {
      nanoflann::KNNResultSet<double> resultSet(1);
      std::size_t index;
      double distance;

      resultSet.init(&index, &distance);
      tree.findNeighbors(resultSet, x_data_.get_row(i));

      // nanoflann L2_Adaptor returns squared Euclidean distance.
      distances[i] = std::sqrt(std::max(0.0, distance));
    }

    x_data_.subset_off();

    return *std::max_element(distances.begin(), distances.end());
  }
};


// [[Rcpp::export]]
Rcpp::List twin_select_global_rcpp(
    Rcpp::NumericMatrix twin_data,
    Rcpp::NumericMatrix Xnorm,
    std::size_t r,
    std::size_t runs,
    Rcpp::IntegerVector u1,
    std::size_t leaf_size = 8)
{
  std::vector<std::size_t> u1_cpp(runs);

  for (std::size_t i = 0; i < runs; i++) {
    u1_cpp[i] = static_cast<std::size_t>(u1[i] - 1);
  }

  BKPTwinKDTree tree(
      twin_data,
      Xnorm,
      r,
      runs,
      u1_cpp,
      leaf_size
  );

  return tree.twin();
}
