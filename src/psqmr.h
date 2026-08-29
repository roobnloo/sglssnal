#ifndef PSQMR_H
#define PSQMR_H

#include <RcppArmadillo.h>
#include <functional>

using MatVecFun = std::function<arma::vec(const arma::vec &)>;

// Preconditioned symmetric QMR, ported from solvers/psqmrGL.m in the
// reference MATLAB implementation (Zhang et al. 2020), itself adapted from
// the SDPNAL package (Xinyuan Zhao, Defeng Sun, Kim-Chuan Toh, 2008). Solves
// matvec(x) = b for a symmetric positive semidefinite operator without ever
// materializing it -- `matvec` is the only way the operator is accessed.
//
// `par` may contain (all optional, matching psqmrGL.m's own defaults):
//   "tol"                  stopping tolerance on the residual norm
//                          (default 1e-6*norm(b))
//   "maxit"                iteration cap (default max(5000, ceil(sqrt(N))))
//   "stagnate_check_psqmr" window size for the stagnation check (default 20)
//   "minitpsqmr"           minimum iterations before accepting convergence
//                          (default 0)
//   "precond"              0 = identity (default), 1 = scalar diagonal
//   "invdiagM"             scalar preconditioner value, required if
//                          precond == 1 -- this codebase only ever produces
//                          a scalar (not a per-coordinate diagonal), so
//                          that's all this port supports
//
// `solve_ok` mirrors psqmrGL.m exactly: 1 converged normally, -1 stagnated,
// -2 hit the iteration cap, 2 breakdown (near-zero pivot).
Rcpp::List psqmrGL(const MatVecFun &matvec, const arma::vec &b,
                    Rcpp::List &par, const arma::vec &x0 = arma::vec());

#endif
