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

bool mat2_ssn(const arma::vec &u, const arma::sp_mat &A, double lam1,
              double lam2, const GroupStruct &gs, double sig, arma::mat &V2,
              arma::sp_mat &D, uint &sp_dim) {
  vec prox1u = proximal_l1(u, lam1);
  vec pv = gs.pma * prox1u;
  int n = A.n_rows;
  vec grp_norms(gs.num_group);
  vec proj2pv = projection_l2(pv, lam2, gs, grp_norms);
  if (!any(prox1u) || !any(grp_norms)) {
    return false;
  }

  uvec nz_prox1u_ids = find(prox1u);
  uvec nz_grp_ids = find(grp_norms);
  sp_mat B(n, nz_prox1u_ids.n_elem + nz_grp_ids.n_elem);
  sp_mat C(n, nz_grp_ids.n_elem);
  uint s_start = 0;
  uint i = 0;
  for (int k : nz_grp_ids) {
    vec vk = gs.get_group_subview(prox1u, k);
    uvec indvk = find(vk);

    double cw = lam2 * gs.ind(2, k);
    double par1 = sig * cw / grp_norms(k);
    double par2 = par1 / (grp_norms(k) * grp_norms(k));

    sp_mat Al = gs.get_group_subview(A, k);
    Al = Al.cols(indvk);

    sp_mat Bl = sqrt(sig - par1) * Al;
    uint lenind1 = Bl.n_cols;
    uint s_end = s_start + lenind1 - 1;
    B.cols(s_start, s_end) = Bl;
    s_start = s_end + 1;

    if (any(abs(vk) > datum::eps)) {
      vec pv = vk(indvk);
      vec Asl = Al * pv;
      C.col(i) = sqrt(par2) * Asl;
      i++;
    }
  }

  sp_dim = s_start + nz_grp_ids.n_elem;
  B.cols(s_start, sp_dim - 1) = C;
  D = B.cols(0, sp_dim - 1);
  V2 = D.t() * D;
  V2.diag() += 1;

  return true;
}

// [[Rcpp::export]]
List conjgrad_linsolver_impl(const arma::sp_mat &A, const arma::vec &rhs,
                             const arma::vec &u, double lam1, double lam2,
                             const arma::sp_mat &pma, const arma::uvec g,
                             const arma::mat &ind, uint num_group, int density,
                             double sig) {
  GroupStruct gs = {pma, g, ind, num_group};
  int n = rhs.n_elem;
  int solver = 1; // 1: direct, 2: direct woodbury-formula, 3: pcg
  int dn = 10000;

  if (n <= 1000) {
    solver = 1;
  }

  if (density <= n && density <= dn) {
    solver = 2;
  }

  arma::vec dy;
  double resnrm = 0;
  int solve_ok = 0;

  if (solver == 1) {
    arma::mat V;
    arma::vec proj2pv, prox1u;
    mat_ssn(u, A, lam1, lam2, gs, sig, V, proj2pv, prox1u);
    dy = solve(V, rhs);
    solve_ok = 1;
  }

  if (solver == 2) {
    arma::mat V2;
    arma::sp_mat D;
    uint sp_dim;
    if (mat2_ssn(u, A, lam1, lam2, gs, sig, V2, D, sp_dim)) {
      arma::vec rhstmp = D.t() * rhs;
      dy = solve(V2, rhstmp);
      dy = D * dy;
      dy = rhs - dy;
    } else {
      dy = rhs;
    }
    solve_ok = 1;
  }

  return List::create(Named("dy") = dy, Named("resnrm") = resnrm,
                      Named("solve_ok") = solve_ok);
}