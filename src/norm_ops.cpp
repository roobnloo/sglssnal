#include "norm_ops.h"

using namespace arma;

arma::vec proximal_l1(const arma::vec &w, double lam) {
  if (lam < 0) {
    Rcpp::stop("lam must be non-negative");
  }
  return clamp(w - lam, 0.0, datum::inf) - clamp(-w - lam, 0.0, datum::inf);
}

arma::vec proximal_l2(const arma::vec &z, double lam, const GroupStruct &gs) {
  if (lam < 0) {
    Rcpp::stop("lam must be non-negative");
  }

  arma::vec result = z;

  for (int i = 0; i < gs.num_group; ++i) {
    int kstart = gs.ind(0, i);
    int kend = gs.ind(1, i);
    int w = lam * gs.ind(2, i);

    double l2norm = norm(z.subvec(kstart, kend), 2);
    if (l2norm > w) {
      result.subvec(kstart, kend) = (1 - w / l2norm) * z.subvec(kstart, kend);
    } else {
      result.subvec(kstart, kend).zeros();
    }
  }

  return result;
}

arma::vec projection_l2(const arma::vec &z, double lam, const GroupStruct &gs,
                        arma::vec &grp_norms) {
  if (lam < 0) {
    Rcpp::stop("lam must be non-negative");
  }
  if (grp_norms.n_elem != gs.num_group) {
    Rcpp::stop(
        "grp_norms must have the same number of elements as gs.num_group");
  }

  arma::vec Pz(size(z));

  for (int i = 0; i < gs.num_group; ++i) {
    int kstart = gs.ind(0, i);
    int kend = gs.ind(1, i);
    double w = lam * gs.ind(2, i);

    arma::vec group_z = z.subvec(kstart, kend);
    double l2norm = norm(group_z, 2);

    if (l2norm > w) {
      Pz.subvec(kstart, kend) = (w / l2norm) * group_z;
      grp_norms(i) = l2norm;
    } else {
      Pz.subvec(kstart, kend) = group_z;
      grp_norms(i) = 0;
    }
  }

  return Pz;
}