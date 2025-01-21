#include "group_struct.h"
#include "norm_ops.h"
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

void mat_ssn(const arma::vec &u, const arma::sp_mat &A, double lam1,
             double lam2, const GroupStruct &gs, double sig, arma::mat &V,
             arma::vec &proj2pv, arma::vec &prox1u) {
  prox1u = proximal_l1(u, lam1);
  vec pv = gs.pma * prox1u;
  int n = A.n_rows;
  vec grp_norms(gs.num_group);
  proj2pv = projection_l2(pv, lam2, gs, grp_norms);
  if (!any(prox1u) || !any(grp_norms)) {
    V = eye(n, n);
  }

  V = zeros(n, n);

  for (uint k = 0; k < gs.num_group; k++) {
    if (grp_norms(k) <= datum::eps) {
      continue;
    }
    vec vk = gs.get_group_subview(prox1u, k);
    double cw = lam2 * gs.ind(2, k);

    double par1 = sig * cw / grp_norms(k);
    double par2 = par1 / (grp_norms(k) * grp_norms(k));

    sp_mat Al = gs.get_group_subview(A, k);
    uvec indvk = find(vk);
    Al = Al.cols(indvk);

    sp_mat M1 = Al * Al.t();
    V += (sig - par1) * M1;

    if (any(abs(vk) > datum::eps)) {
      vec pv = vk(indvk);
      vec Asl = Al * pv;
      mat M2 = Asl * Asl.t();
      V += par2 * M2;
    }
  }
  V.diag() += 1;
}

// [[Rcpp::export]]
void mat_ssn_interface(const arma::vec &u, const arma::sp_mat &A, double lam1,
                       double lam2, double sig, const arma::sp_mat &pma,
                       const arma::uvec g, const arma::mat &ind, uint num_group,
                       arma::mat &V, arma::vec &proj2pv, arma::vec &prox_vec) {
  GroupStruct gs = {pma, g, ind, num_group};
  mat_ssn(u, A, lam1, lam2, gs, sig, V, proj2pv, prox_vec);
}