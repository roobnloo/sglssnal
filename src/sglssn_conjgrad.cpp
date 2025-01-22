#include "group_struct.h"
#include "norm_ops.h"
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

List findstep_impl(const arma::vec &b, double sig, double psi_y0,
                   const arma::vec &u0, const arma::vec &Prox_u0,
                   const arma::vec &sigProx_u0, const arma::vec &z0,
                   const arma::vec &y0, const arma::vec &Aty0,
                   const arma::vec &dy, const arma::vec &Atdy,
                   const arma::vec &lambda, GroupStruct &gs, double tol,
                   int stepop) {
  int printlevel = 1;
  int maxit = std::ceil(std::log(1 / (tol + datum::eps)) / std::log(2.0));
  double c1 = 1e-4;
  double c2 = 0.9;

  double tmp1 = dot(dy, y0 + b);
  double tmp2 = sum(square(dy));

  double cross = dot(Atdy, Prox_u0);
  double g0 = sig * cross - tmp1;
  if (g0 <= 0) {
    double alp = 0;
    if (printlevel) {
      Rcpp::Rcout << "\n Need an ascent direction, " << g0 << "  " << std::endl;
    }
    return List::create(
        Named("psi_y") = psi_y0, Named("u") = u0, Named("Prox_u") = Prox_u0,
        Named("sigProx_u") = sigProx_u0, Named("z") = z0, Named("y") = y0,
        Named("Aty") = Aty0, Named("alp") = alp, Named("iter") = 0);
  }

  double alp = 1, alpconst = 0.5, LB = 0, UB = 1, psi_y = psi_y0;
  double gLB, gUB;
  vec u = u0, Prox_u = Prox_u0, sigProx_u = sigProx_u0, z = z0, y = y0,
      Aty = Aty0;

  int iter = 0;
  while (iter < maxit) {
    if (iter > 0) {
      alp = alpconst * (LB + UB);
    }
    y = y0 + alp * dy;
    u = u0 - alp * Atdy;
    Prox_u = proximal_combo(u, lambda[0], lambda[1], gs);
    z = u - Prox_u;
    double galp = tmp1 + alp * tmp2;
    sigProx_u = sig * Prox_u;
    galp = dot(Atdy, sigProx_u) - galp;
    psi_y =
        -(dot(b, y) + 0.5 * sum(square(y)) + 0.5 * sig * sum(square(Prox_u)));

    if (printlevel > 1) {
      Rcpp::Rcout << "\n --------------------------------------- " << std::endl;
      Rcpp::Rcout << " alp = " << alp << ", psi_y = " << psi_y
                  << ", psi_y0 = " << psi_y0 << std::endl;
      Rcpp::Rcout << " --------------------------------------- " << std::endl;
    }

    if (iter == 0) {
      gLB = g0;
      gUB = galp;
      if (sign(gLB) * sign(gUB) > 0) {
        if (printlevel)
          Rcpp::Rcout << "|";
        Aty = Aty0 + alp * Atdy;
        return List::create(
            Named("psi_y") = psi_y, Named("u") = u, Named("Prox_u") = Prox_u,
            Named("sigProx_u") = sigProx_u, Named("z") = z, Named("y") = y,
            Named("Aty") = Aty, Named("alp") = alp, Named("iter") = iter);
      }
    }

    if (std::abs(galp) < c2 * std::abs(g0) &&
        (psi_y - psi_y0 - c1 * alp * g0 >
         1e-8 / std::max(1.0, std::abs(psi_y0)))) {
      if (stepop == 1 || (stepop == 2 && std::abs(galp) < tol)) {
        if (printlevel)
          Rcpp::Rcout << ":";
        Aty = Aty0 + alp * Atdy;
        return List::create(
            Named("psi_y") = psi_y, Named("u") = u, Named("Prox_u") = Prox_u,
            Named("sigProx_u") = sigProx_u, Named("z") = z, Named("y") = y,
            Named("Aty") = Aty, Named("alp") = alp, Named("iter") = iter);
      }
    }

    if (sign(galp) * sign(gUB) < 0) {
      LB = alp;
      gLB = galp;
    } else if (sign(galp) * sign(gLB) < 0) {
      UB = alp;
      gUB = galp;
    }
    iter++;
  }

  if (iter == maxit) {
    Aty = Aty0 + alp * Atdy;
  }
  if (printlevel)
    Rcpp::Rcout << "m";
  return List::create(Named("psi_y") = psi_y, Named("u") = u,
                      Named("Prox_u") = Prox_u, Named("sigProx_u") = sigProx_u,
                      Named("z") = z, Named("y") = y, Named("Aty") = Aty,
                      Named("alp") = alp, Named("iter") = iter);
}

// [[Rcpp::export]]
List findstep_interface(const arma::vec &b, double sig, double psi_y0,
                        const arma::vec &u0, const arma::vec &Prox_u0,
                        const arma::vec &sigProx_u0, const arma::vec &z0,
                        const arma::vec &y0, const arma::vec &Aty0,
                        const arma::vec &dy, const arma::vec &Atdy,
                        const arma::vec &lambda, const arma::sp_mat &pma,
                        const arma::uvec &g, const arma::mat &ind,
                        uint num_group, double tol, int stepop) {

  GroupStruct gs = {pma, g, ind, num_group};
  return findstep_impl(b, sig, psi_y0, u0, Prox_u0, sigProx_u0, z0, y0, Aty0,
                       dy, Atdy, lambda, gs, tol, stepop);
}