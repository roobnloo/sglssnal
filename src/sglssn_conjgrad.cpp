#include "sglssn_conjgrad.h"
#include "conjgrad_linsolver.h"
#include "norm_ops.h"
#include <RcppArmadillo.h>
#include <ios>
using namespace Rcpp;
using namespace arma;

List findstep_impl(const arma::vec &b, double sig, double psi_y0,
                   const arma::vec &u0, const arma::vec &Prox_u0,
                   const arma::vec &sigProx_u0, const arma::vec &z0,
                   const arma::vec &y0, const arma::vec &Aty0,
                   const arma::vec &dy, const arma::vec &Atdy, double lam1,
                   double lam2, const GroupStruct &gs, double tol, int stepop,
                   bool printsub) {

  if (printsub) {
    Rcout << std::scientific << std::setprecision(2) << std::setw(3);
  }
  int maxit = std::ceil(std::log(1 / (tol + datum::eps)) / std::log(2.0));
  double c1 = 1e-4;
  double c2 = 0.9;

  double tmp1 = dot(dy, y0 + b);
  double tmp2 = sum(square(dy));

  double cross = dot(Atdy, Prox_u0);
  double g0 = sig * cross - tmp1;
  if (g0 <= 0) {
    double alp = 0;
    if (printsub) {
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
    Prox_u = proximal_combo(u, lam1, lam2, gs);
    z = u - Prox_u;
    double galp = tmp1 + alp * tmp2;
    sigProx_u = sig * Prox_u;
    galp = dot(Atdy, sigProx_u) - galp;
    psi_y =
        -(dot(b, y) + 0.5 * sum(square(y)) + 0.5 * sig * sum(square(Prox_u)));

    // if (printsub) {
    //   Rcpp::Rcout << "\n --------------------------------------- " <<
    //   std::endl; Rcpp::Rcout << " alp = " << alp << ", psi_y = " << psi_y
    //               << ", psi_y0 = " << psi_y0 << std::endl;
    //   Rcpp::Rcout << " --------------------------------------- " <<
    //   std::endl;
    // }

    if (iter == 0) {
      gLB = g0;
      gUB = galp;
      if (sign(gLB) * sign(gUB) > 0) {
        if (printsub)
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
        if (printsub)
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
  if (printsub)
    Rcpp::Rcout << "m";
  return List::create(Named("psi_y") = psi_y, Named("u") = u,
                      Named("Prox_u") = Prox_u, Named("sigProx_u") = sigProx_u,
                      Named("z") = z, Named("y") = y, Named("Aty") = Aty,
                      Named("alp") = alp, Named("iter") = iter);
}

List sglssn_conjgrad(const arma::vec &y0, const arma::vec &Aty0,
                     const arma::vec &x0, const arma::vec &Ax0,
                     const arma::sp_mat &A, const arma::vec &b, double lam1,
                     double lam2, const GroupStruct &gs, List &par,
                     bool printsub, int maxitersub, double tol) {
  int breakyes = 0;
  double tiny = 1e-10;
  int maxitpsqmr = 500;
  int precond = 0;
  double const3 = 0.7;

  double sig = as<double>(par["sigma"]);
  double tolconst = as<double>(par["tolconst"]);
  double normb = 1 + norm(b, 2);

  arma::vec y = y0;
  arma::vec Aty = Aty0;
  arma::vec u = x0 / sig - Aty;
  arma::vec Prox_u = proximal_combo(u, lam1, lam2, gs);
  arma::vec x = sig * Prox_u;
  arma::vec Ax;
  arma::vec z = u - Prox_u;
  double psi_y =
      -(dot(b, y) + 0.5 * sum(square(y)) + 0.5 * sig * sum(square(Prox_u)));

  NumericVector priminf_vec(maxitersub), dualinf_vec(maxitersub),
      psiy_vec(maxitersub);
  IntegerVector solve_ok_vec(maxitersub), findstep_vec(maxitersub),
      psqmr_vec(maxitersub);
  double av_findstep;

  for (int itersub = 1; itersub <= maxitersub; ++itersub) {
    Ax = A * x;
    arma::vec Grad = -y - b + Ax;
    double normGrad = norm(Grad, 2);
    double normGrad_b = normGrad / normb;
    double &priminf_sub = normGrad_b;
    arma::vec Rd = Aty + z;
    double normRd = norm(Rd, 2);
    double dualinf_sub = normRd / (1 + norm(z, 2));

    double tolsubconst = (priminf_sub < tol && dualinf_sub < tol) ? 0.9 : 1e-2;
    double tolsub =
        std::max(std::min(1.0, tolconst * dualinf_sub), tolsubconst * tol);

    priminf_vec[itersub - 1] = priminf_sub;
    dualinf_vec[itersub - 1] = dualinf_sub;
    psiy_vec[itersub - 1] = psi_y;

    if (printsub) {
      Rcpp::Rcout << "\n      " << itersub << "  " << psi_y << "  "
                  << priminf_sub << "   " << dualinf_sub << " " << tolconst;
    }

    if (normGrad_b < tolsub && itersub > 1) {
      if (printsub) {
        Rcpp::Rcout << "\n       good termination in subproblem: "
                    << " dualinfes = " << dualinf_sub
                    << ", normGrad = " << normGrad_b << ", tolsub = " << tolsub
                    << std::endl;
      }
      breakyes = -1;
      break;
    }

    par["epsilon"] = std::min(1e-3, 0.1 * normGrad_b);
    par["precond"] = precond;
    if (precond == 1) {
      par["invdiagM"] = 1 / (1 + sig);
    }

    if (dualinf_sub > 1e-3 || itersub <= 5) {
      maxitpsqmr = std::max(maxitpsqmr, 200);
    } else if (dualinf_sub > 1e-4) {
      maxitpsqmr = std::max(maxitpsqmr, 300);
    } else if (dualinf_sub > 1e-5) {
      maxitpsqmr = std::max(maxitpsqmr, 400);
    } else if (dualinf_sub > 5e-6) {
      maxitpsqmr = std::max(maxitpsqmr, 500);
    }

    double prim_ratio =
        (itersub > 1) ? priminf_sub / priminf_vec[itersub - 2] : 0;
    double dual_ratio =
        (itersub > 1) ? dualinf_sub / dualinf_vec[itersub - 2] : 0;

    arma::vec &rhs = Grad;
    double tolpsqmr = std::min(5e-3, 0.001 * normGrad);
    double const2 = 1;
    if (itersub > 1 &&
        (prim_ratio > 0.5 || priminf_sub > 0.1 * priminf_vec[0])) {
      const2 *= 0.5;
    }
    if (dual_ratio > 1.1)
      const2 *= 0.5;
    tolpsqmr *= const2;
    par["tol"] = tolpsqmr;
    par["maxit"] = maxitpsqmr;
    int nnz = sum(abs(x) > tol);
    par["nnz"] = nnz;

    List result = conjgrad_linsolver(A, rhs, u, lam1, lam2, gs, nnz, sig);
    arma::vec dy = as<arma::vec>(result["dy"]);
    arma::vec resnrm = as<arma::vec>(result["resnrm"]);
    int solve_ok = as<int>(result["solve_ok"]);

    arma::vec Atdy = A.t() * dy;
    int iterpsqmr = resnrm.n_elem - 1;
    if (printsub) {
      Rcpp::Rcout << " | " << tolpsqmr << " " << resnrm[resnrm.n_elem - 1]
                  << " " << iterpsqmr << "  " << const2 << " " << nnz;
    }

    par["iter"] = itersub;

    int stepop = (itersub <= 3 && dualinf_sub > 1e-4) || (itersub < 3) ? 1 : 2;
    double steptol = 1e-5;

    List step_result =
        findstep_impl(b, sig, psi_y, u, Prox_u, x, z, y, Aty, dy, Atdy, lam1,
                      lam2, gs, steptol, stepop, printsub);
    psi_y = as<double>(step_result["psi_y"]);
    u = as<arma::vec>(step_result["u"]);
    Prox_u = as<arma::vec>(step_result["Prox_u"]);
    x = as<arma::vec>(step_result["sigProx_u"]);
    z = as<arma::vec>(step_result["z"]);
    y = as<arma::vec>(step_result["y"]);
    Aty = as<arma::vec>(step_result["Aty"]);
    double alp = as<double>(step_result["alp"]);
    int iterstep = as<int>(step_result["iter"]);

    solve_ok_vec[itersub - 1] = solve_ok;
    psqmr_vec[itersub - 1] = iterpsqmr;
    findstep_vec[itersub - 1] = iterstep;
    av_findstep = mean(findstep_vec);

    if (alp < tiny)
      breakyes = 11;
    double psiy_ratio =
        (itersub > 1)
            ? (psi_y - psiy_vec[itersub - 2]) /
                  (std::abs(psi_y) + std::numeric_limits<double>::epsilon())
            : 1;

    if (printsub) {
      Rcpp::Rcout << " " << alp << " " << iterstep;
      if (psiy_ratio < 0)
        Rcpp::Rcout << "-";
    }

    if (itersub > 4) {
      IntegerVector idx = seq(std::max(0, itersub - 4), itersub - 1);
      NumericVector tmp = priminf_vec[idx];
      double ratio = min(tmp) / max(tmp);

      IntegerVector solve_ok_subvec = solve_ok_vec[idx];
      IntegerVector psqmr_subvec = psqmr_vec[idx];
      bool all_stagnant =
          std::all_of(solve_ok_subvec.begin(), solve_ok_subvec.end(),
                      [](int x) { return x <= -1; });
      if (all_stagnant && ratio > 0.9 &&
          min(psqmr_subvec) == max(psqmr_subvec) && max(tmp) < 5 * tol) {
        Rcpp::Rcout << "#";
        breakyes = 1;
      }

      double priminf_1half =
          min(priminf_vec[seq(0, ceil(itersub * const3) - 1)]);
      double priminf_2half =
          min(priminf_vec[seq(std::ceil(itersub * const3), itersub - 1)]);
      double priminf_best = min(priminf_vec[seq(0, itersub - 2)]);
      double priminf_ratio =
          priminf_vec[itersub - 1] / priminf_vec[itersub - 2];
      solve_ok_subvec = solve_ok_vec[seq(0, itersub - 1)];
      int stagnate_count =
          std::count_if(solve_ok_subvec.begin(), solve_ok_subvec.end(),
                        [](int x) { return x <= -1; });
      IntegerVector idx2 = seq(std::max(0, itersub - 8), itersub - 1);
      solve_ok_subvec = solve_ok_vec[idx2];
      all_stagnant = std::all_of(solve_ok_subvec.begin(), solve_ok_subvec.end(),
                                 [](int x) { return x == -1; });

      if (itersub >= 10 && all_stagnant && priminf_best < 1e-2 &&
          dualinf_sub < 1e-3) {
        tmp = priminf_vec[idx2];
        ratio = min(tmp) / max(tmp);
        if (ratio > 0.5) {
          if (printsub)
            Rcpp::Rcout << "##";
          breakyes = 2;
        }
      }

      if (itersub >= 15 && priminf_1half < std::min(2e-3, priminf_2half) &&
          dualinf_sub < 0.8 * dualinf_vec[0] && dualinf_sub < 1e-3 &&
          stagnate_count >= 3) {
        if (printsub)
          Rcpp::Rcout << "###";
        breakyes = 3;
      }

      if (itersub >= 15 && priminf_ratio < 0.1 &&
          priminf_sub < 0.8 * priminf_1half &&
          dualinf_sub < std::min(1e-3, 2 * priminf_sub) &&
          (priminf_sub < 2e-3 || (dualinf_sub < 1e-5 && priminf_sub < 5e-3)) &&
          stagnate_count >= 3) {
        if (printsub)
          Rcpp::Rcout << " $$";
        breakyes = 4;
      }

      if (itersub >= 10 && dualinf_sub > 5 * min(dualinf_vec) &&
          priminf_sub > 2 * min(priminf_vec)) {
        if (printsub)
          Rcpp::Rcout << "$$$";
        breakyes = 5;
      }

      if (itersub >= 20) {
        vec dualinf_ratioall =
            dualinf_vec[seq(1, itersub - 1)] / dualinf_vec[seq(0, itersub - 2)];
        arma::uvec idx = find(dualinf_ratioall > 1);
        if (idx.n_elem >= 3) {
          double dualinf_increment = mean(dualinf_ratioall.elem(idx));
          if (dualinf_increment > 1.25) {
            if (printsub)
              Rcpp::Rcout << "^^";
            breakyes = 6;
          }
        }
      }

      if (breakyes > 0) {
        x = sig * Prox_u;
        Ax = A * x;
        break;
      }
    }
  }

  List info;
  info["maxCG"] = (int)max(psqmr_vec);
  info["avgCG"] = (double)sum(psqmr_vec) / maxitersub;
  info["breakyes"] = breakyes;
  info["itersub"] = maxitersub;
  info["tolconst"] = tolconst;

  List runhist = List::create(
      Named("priminf") = priminf_vec, Named("dualinf") = dualinf_vec,
      Named("psiy") = psiy_vec, Named("findstep") = findstep_vec,
      Named("psqmr") = psqmr_vec, Named("solve_ok") = solve_ok_vec,
      Named("av_findstep") = av_findstep);

  return List::create(Named("y") = y, Named("z") = z, Named("Aty") = Aty,
                      Named("Prox_u") = Prox_u, Named("x") = x,
                      Named("Ax") = Ax, Named("par") = par,
                      Named("runhist") = runhist, Named("info") = info);
}
