#include "conjgrad_linsolver.h"
#include "group_struct.h"
#include "norm_ops.h"
#include "psqmr.h"
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

template <typename MatType>
void mat_ssn(const arma::vec &u, const MatType &A, double lam1, double lam2,
             const GroupStruct &gs, double sig, arma::mat &V,
             arma::vec &proj2pv, arma::vec &prox1u) {
  prox1u = proximal_l1(u, lam1);
  vec pv = gs.pma * prox1u;
  int n = A.n_rows;
  int nvars = A.n_cols;
  vec grp_norms(gs.num_group);
  proj2pv = projection_l2(pv, lam2, gs, grp_norms);
  if (!any(prox1u) || !any(grp_norms)) {
    V = eye(n, n);
    return;
  }

  V = zeros(n, n);
  vec DD = zeros(nvars);

  for (uword k = 0; k < gs.num_group; k++) {
    if (grp_norms(k) <= datum::eps) {
      continue;
    }
    vec vk = gs.get_group_subview(prox1u, k);
    double cw = lam2 * gs.ind(2, k);

    uvec subview_col_ids = gs.elem_ids.subvec(gs.ind(0, k), gs.ind(1, k));
    uvec indvk = find(vk);
    subview_col_ids = subview_col_ids(indvk);
    double par1 = sig * cw / grp_norms(k);
    DD(subview_col_ids) += sig - par1;

    if (any(abs(vk) > datum::eps)) {
      vec Asl = A.cols(subview_col_ids) * vk(indvk);
      mat M2 = Asl * Asl.t();
      V += par1 / (grp_norms(k) * grp_norms(k)) * M2;
    }
  }
  V.diag() += 1;
  uvec indDD = find(DD);
  MatType A_dd = A.cols(indDD);
  V += A_dd * diagmat(DD(indDD)) * A_dd.t();
}

// Builds the sparse Woodbury factor D (n x sp_dim) shared by solver 2
// (direct Woodbury) and solver 3 (PCG on the reduced system): D^T D + I is
// the reduced system either gets solved directly (solver 2) or iteratively
// via matrix-free matvecs (solver 3). Returns false (matching mat2_ssn's/
// MATLAB's inverted-"id_yes" convention) when there's no active correction
// to apply, i.e. the Newton system is just the identity.
template <typename MatType>
bool build_woodbury_D(const arma::vec &u, const MatType &A, double lam1,
                       double lam2, const GroupStruct &gs, double sig,
                       sp_mat &D, uword &sp_dim) {
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
  uword s_start = 0;
  uword i = 0;
  for (int k : nz_grp_ids) {
    vec vk = gs.get_group_subview(prox1u, k);
    uvec indvk = find(vk);

    double cw = lam2 * gs.ind(2, k);
    double par1 = sig * cw / grp_norms(k);
    double par2 = par1 / (grp_norms(k) * grp_norms(k));

    MatType Al = gs.get_group_subview(A, k);
    Al = Al.cols(indvk);

    MatType Bl = sqrt(sig - par1) * Al;
    uword lenind1 = Bl.n_cols;
    uword s_end = s_start + lenind1 - 1;
    B.cols(s_start, s_end) = Bl;
    s_start = s_end + 1;

    if (any(abs(vk) > datum::eps)) {
      vec vk_active = vk(indvk);
      vec Asl = Al * vk_active;
      C.col(i) = sqrt(par2) * Asl;
      i++;
    }
  }

  sp_dim = s_start + nz_grp_ids.n_elem;
  B.cols(s_start, sp_dim - 1) = C;
  D = B.cols(0, sp_dim - 1);

  return true;
}

template <typename MatType>
bool mat2_ssn(const arma::vec &u, const MatType &A, double lam1, double lam2,
              const GroupStruct &gs, double sig, arma::mat &V2, sp_mat &D,
              uword &sp_dim) {
  if (!build_woodbury_D(u, A, lam1, lam2, gs, sig, D, sp_dim)) {
    return false;
  }
  V2 = D.t() * D;
  V2.diag() += 1;
  return true;
}

template <typename MatType>
List conjgrad_linsolver(const MatType &A, const arma::vec &rhs,
                        const arma::vec &u, double lam1, double lam2,
                        const GroupStruct &gs, int density, double sig,
                        List &par) {
  int n = rhs.n_elem;
  // 1: direct, 2: direct woodbury-formula, 3: pcg. Default matches the
  // MATLAB reference's Linsolver_CG.m exactly (it defaults to 3, not 1) --
  // without this, the largest/densest problems that neither of the two
  // `if`s below reassigns would otherwise silently fall through to the
  // most expensive dense path, backwards from intent.
  int solver = 3;
  int dn = 10000;

  if (n <= 1000) {
    solver = 1;
  }

  if (density <= n && density <= dn) {
    solver = 2;
  }

  // Ported verbatim from Linsolver_CG.m's third dispatch condition.
  if ((n > 5000 && density >= 1000) || (n > 2000 && density > 5000) ||
      (n > 100 && density > 8000)) {
    solver = 3;
  }

  arma::vec dy;
  arma::vec resnrm = {0.0};
  int solve_ok = 0;

  if (solver == 1) {
    arma::mat V;
    arma::vec proj2pv, prox1u;
    mat_ssn(u, A, lam1, lam2, gs, sig, V, proj2pv, prox1u);
    dy = solve(V, rhs, solve_opts::likely_sympd);
    solve_ok = 1;
  }

  if (solver == 2) {
    arma::mat V2;
    arma::sp_mat D;
    uword sp_dim;
    if (mat2_ssn(u, A, lam1, lam2, gs, sig, V2, D, sp_dim)) {
      arma::vec rhstmp = D.t() * rhs;
      dy = solve(V2, rhstmp, solve_opts::likely_sympd);
      dy = D * dy;
      dy = rhs - dy;
    } else {
      dy = rhs;
    }
    solve_ok = 1;
  }

  if (solver == 3) {
    arma::sp_mat D;
    uword sp_dim;
    if (!build_woodbury_D(u, A, lam1, lam2, gs, sig, D, sp_dim)) {
      dy = rhs;
      solve_ok = 1;
    } else {
      arma::vec rhstmp = D.t() * rhs;
      // Never materializes anything larger than the sparse D built above --
      // the whole point of this path. Reduced system is D^T D + I, applied
      // matrix-free.
      MatVecFun Afun = [&D](const arma::vec &x) -> arma::vec {
        return x + D.t() * (D * x);
      };
      List psqmr_result = psqmrGL(Afun, rhstmp, par);
      arma::vec dy_reduced = as<arma::vec>(psqmr_result["x"]);
      resnrm = as<arma::vec>(psqmr_result["resnrm"]);
      solve_ok = as<int>(psqmr_result["solve_ok"]);
      dy = rhs - D * dy_reduced;
    }
  }

  return List::create(Named("dy") = dy, Named("resnrm") = resnrm,
                      Named("solve_ok") = solve_ok);
}

// Explicit instantiations
template List conjgrad_linsolver<sp_mat>(const sp_mat &A, const vec &rhs,
                                         const vec &u, double lam1, double lam2,
                                         const GroupStruct &gs, int density,
                                         double sig, List &par);

template List conjgrad_linsolver<mat>(const mat &A, const vec &rhs,
                                      const vec &u, double lam1, double lam2,
                                      const GroupStruct &gs, int density,
                                      double sig, List &par);

// The two interfaces below exist solely to give testthat a seam onto this
// otherwise-internal C++ routine (no @export roxygen tag, so NAMESPACE
// never lists them -- reachable only via sglssnal:::, mirroring how
// sglssnal_main_interface/_dense are used elsewhere in this package).
// [[Rcpp::export]]
List conjgrad_linsolver_interface(const arma::sp_mat &A, const arma::vec &rhs,
                                  const arma::vec &u, double lam1,
                                  double lam2, const List &gs_list,
                                  int density, double sig, List par) {
  uvec G = as<uvec>(gs_list["G"]);
  mat ind = as<mat>(gs_list["ind"]);
  uword num_group = ind.n_cols;
  GroupStruct gs = {as<sp_mat>(gs_list["pma"]), G, ind, num_group};
  return conjgrad_linsolver(A, rhs, u, lam1, lam2, gs, density, sig, par);
}

// [[Rcpp::export]]
List conjgrad_linsolver_interface_dense(const arma::mat &A,
                                        const arma::vec &rhs,
                                        const arma::vec &u, double lam1,
                                        double lam2, const List &gs_list,
                                        int density, double sig, List par) {
  uvec G = as<uvec>(gs_list["G"]);
  mat ind = as<mat>(gs_list["ind"]);
  uword num_group = ind.n_cols;
  GroupStruct gs = {as<sp_mat>(gs_list["pma"]), G, ind, num_group};
  return conjgrad_linsolver(A, rhs, u, lam1, lam2, gs, density, sig, par);
}