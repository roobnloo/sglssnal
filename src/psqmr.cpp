#include "psqmr.h"

using namespace Rcpp;
using namespace arma;

namespace {

// Mirrors precondfun in psqmrGL.m, restricted to modes 0 (identity) and 1
// (scalar diagonal) -- see psqmr.h for why mode 1 is scalar-only here.
// Modes 2 (low-rank) and 4 (Cholesky-based) are not ported: nothing in this
// codebase produces the state (par$V/par$Vt/par$d, or a factor object) they
// need, and mode 3 (arbitrary function handle) has no analog here since A
// is always materialized, unlike upstream's matrix-free Amap/ATmap path.
vec precondfun(const List &par, const vec &r) {
  int precond = 0;
  if (par.containsElementNamed("precond")) {
    precond = as<int>(par["precond"]);
  }
  if (precond == 1) {
    double invdiagM = as<double>(par["invdiagM"]);
    return invdiagM * r;
  }
  return r;
}

} // namespace

List psqmrGL(const MatVecFun &matvec, const vec &b, List &par,
             const vec &x0) {
  int N = b.n_elem;
  int maxit = std::max(5000, (int)std::ceil(std::sqrt((double)N)));
  double tol = 1e-6 * norm(b, 2);
  int stagnate_check = 20;
  int miniter = 0;

  if (par.containsElementNamed("maxit")) {
    maxit = as<int>(par["maxit"]);
  }
  if (par.containsElementNamed("tol")) {
    tol = as<double>(par["tol"]);
  }
  if (par.containsElementNamed("stagnate_check_psqmr")) {
    stagnate_check = as<int>(par["stagnate_check_psqmr"]);
  }
  if (par.containsElementNamed("minitpsqmr")) {
    miniter = as<int>(par["minitpsqmr"]);
  }

  int solve_ok = 1;
  const double tiny = 1e-30;

  vec x = (x0.n_elem == (uword)N) ? x0 : zeros<vec>(N);

  vec Aq = (norm(x, 2) > 0) ? matvec(x) : zeros<vec>(N);
  vec r = b - Aq;
  double err = norm(r, 2);
  double minres = err;

  // resnrm_vec(0) is the initial residual norm; resnrm_vec(k) is the
  // residual norm after iteration k (1-based, matching psqmrGL.m's own
  // resnrm(k+1) 1-based storage exactly, just shifted to 0-based indexing).
  vec resnrm_vec = zeros<vec>(maxit + 1);
  resnrm_vec(0) = err;

  vec q = precondfun(par, r);
  double tau_old = norm(q, 2);
  double rho_old = dot(r, q);
  double theta_old = 0;
  vec d = zeros<vec>(N);
  vec res = r;
  vec Ad = zeros<vec>(N);

  int iter = 0;
  int iter_final = 0;
  for (iter = 1; iter <= maxit; iter++) {
    iter_final = iter;
    Aq = matvec(q);
    double sigma = dot(q, Aq);
    if (std::abs(sigma) < tiny) {
      solve_ok = 2;
      break;
    }
    double alpha = rho_old / sigma;
    r = r - alpha * Aq;

    vec u = precondfun(par, r);
    double theta = norm(u, 2) / tau_old;
    double c = 1 / std::sqrt(1 + theta * theta);
    double tau = tau_old * theta * c;
    double gam = c * c * theta_old * theta_old;
    double eta = c * c * alpha;
    d = gam * d + eta * q;
    x = x + d;

    Ad = gam * Ad + eta * Aq;
    res = res - Ad;
    err = norm(res, 2);
    resnrm_vec(iter) = err;
    if (err < minres) {
      minres = err;
    }

    if (err < tol && iter > miniter && dot(b, x) > 0) {
      break;
    }

    if (iter > stagnate_check && iter > 10) {
      vec numer = resnrm_vec.subvec(iter - 10, iter);
      vec denom = resnrm_vec.subvec(iter - 11, iter - 1);
      vec ratio = numer / denom;
      if (ratio.min() > 0.997 && ratio.max() < 1.003) {
        solve_ok = -1;
        break;
      }
    }

    if (std::abs(rho_old) < tiny) {
      solve_ok = 2;
      break;
    }
    double rho = dot(r, u);
    double beta = rho / rho_old;
    q = u + beta * q;

    rho_old = rho;
    tau_old = tau;
    theta_old = theta;
  }

  // Matches psqmrGL.m's own post-loop check exactly: this overwrites
  // whatever solve_ok was set to inside the loop if the cap was reached on
  // the same iteration a break condition also fired -- a faithful port,
  // not a fix, of upstream's behavior.
  if (iter_final == maxit) {
    solve_ok = -2;
  }

  vec Ax = b - res;
  vec resnrm = resnrm_vec.head(iter_final + 1);

  return List::create(Named("x") = x, Named("Ax") = Ax,
                       Named("resnrm") = resnrm, Named("solve_ok") = solve_ok);
}
