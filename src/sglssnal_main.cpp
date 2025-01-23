#include "group_struct.h"
#include "norm_ops.h"
#include "sglssn_conjgrad.h"
#include <RcppArmadillo.h>
#include <chrono>
#include <iomanip>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

using TimePoint = std::chrono::time_point<std::chrono::system_clock>;
TimePoint time_now() { return std::chrono::system_clock::now(); }
double time_diff(TimePoint end, TimePoint start) {
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  return duration.count() / 1000.0;
}

int sigma_fun(int iter) {
  if (iter < 10) {
    return 2;
  } else if (iter < 20) {
    return 3;
  } else if (iter < 200) {
    return 3;
  } else {
    return 10;
  }
}

List sglssnal_main(const arma::sp_mat &A, const arma::vec &b, double lam1,
                   double lam2, const GroupStruct &gs, const List &parmain,
                   const arma::vec &y0, const arma::vec &z0,
                   const arma::vec &x0) {
  double Lip = as<double>(parmain["Lip"]);
  int maxit = as<int>(parmain["maxit"]);
  bool printyes = as<bool>(parmain["printyes"]);
  bool printsub = as<bool>(parmain["printsub"]);
  double stoptol = as<double>(parmain["stoptol"]);
  double stoptol_gap = as<double>(parmain["stoptol_gap"]);
  int stopopt = as<int>(parmain["stopopt"]);
  int p = as<int>(parmain["p"]);
  int n = as<int>(parmain["n"]);
  bool stop = false;

  TimePoint tstart = time_now();
  double ttime;

  double sigmaLip = 1 / Lip;

  vec y = y0, z = z0, x = x0;

  arma::vec Aty = A.t() * y;
  arma::vec Ax = A * x;
  arma::vec Px = gs.pma * x;

  arma::vec obj(2);
  obj[0] = 0.5 * sum(square(Ax - b)) + lam2 * group_l2_norm(Px, gs) +
           lam1 * norm(x, 1);
  obj[1] = 0;
  double normb = 1 + norm(b, 2);
  double sig = std::max(1 / sqrt(Lip),
                        std::min({1.0, sigmaLip, lam1, 1 / lam1, 1 / lam2}));
  arma::vec Rp = Ax - b - y;
  arma::vec Rd = Aty + z;
  double primfeas = norm(Rp, 2) / normb;
  double dualfeas = norm(Rd, 2) / (1 + norm(z, 2));
  double maxfeas = std::max(primfeas, dualfeas);
  double relgap =
      std::abs(obj[0] - obj[1]) / (1 + std::abs(obj[0]) + std::abs(obj[1]));

  if (printyes) {
    Rcout << "\n n=" << p << ", m=" << n << ", tol=" << stoptol
          << ", parameters:c1=" << lam1 << ", c2=" << lam2 << "\n";
    Rcout << " ---------------------------------------------------\n";
    Rcout << " iter|  [pinfeas   dinfeas]    relgap |    pobj      dobj    "
          << " | time(s)  |   sigma  |\n";
    Rcout << "*****************************************************************"
             "***********************************\n";
    Rcout << std::scientific << std::setprecision(2) << std::setw(3) << "    "
          << 0 << "|  [" << primfeas << " " << dualfeas << "]  " << relgap
          << " | " << std::setprecision(4) << obj[0] << " " << obj[1] << " | "
          << std::setprecision(2) << std::setw(3)
          << time_diff(time_now(), tstart) << " | " << sig << " |";
  }

  List par_ncg;
  par_ncg["tolconst"] = 0.5;
  par_ncg["p"] = p;
  int maxitersub = 10;
  int prim_win = 0;
  int dual_win = 0;

  List ssncgop;
  ssncgop["tol"] = stoptol;
  ssncgop["precond"] = 0;
  ssncgop["printsub"] = printsub;

  double sigmamax = 1e6;
  double sigmamin = 1e-4;

  List runhist;
  NumericVector primfeas_vec, dualfeas_vec, sigma_vec, primobj_vec, dualobj_vec,
      gap_vec, relgap_vec, ttime_vec, ratio_seq_vec;
  IntegerVector itersub_vec;

  double eta_1 = 0, eta_tmp = 0;
  int final_iter = -1;
  for (int iter = 1; iter <= maxit; ++iter) {
    par_ncg["sigma"] = sig;

    if (dualfeas < 1e-5) {
      maxitersub = std::max(maxitersub, 30);
    } else if (dualfeas < 1e-3) {
      maxitersub = std::max(maxitersub, 30);
    } else if (dualfeas < 1e-1) {
      maxitersub = std::max(maxitersub, 20);
    }
    ssncgop["maxitersub"] = maxitersub;

    List result = sglssn_conjgrad(y, Aty, x, Ax, A, b, lam1, lam2, gs, par_ncg,
                                  ssncgop["printsub"], ssncgop["maxitersub"],
                                  ssncgop["tol"]);
    y = as<arma::vec>(result["y"]);
    z = as<arma::vec>(result["z"]);
    x = as<arma::vec>(result["x"]);
    Ax = as<arma::vec>(result["Ax"]);
    Aty = as<arma::vec>(result["Aty"]);

    par_ncg = result["par"];
    List runhist_ncg = result["runhist"];
    List info_ncg = result["info"];
    int breakyes = as<int>(info_ncg["breakyes"]);

    if (breakyes < 0) {
      par_ncg["tolconst"] =
          std::max(as<double>(par_ncg["tolconst"]) / 1.06, 1e-3);
    }

    Rd = Aty + z;
    Px = gs.pma * x;
    double normRd = norm(Rd, 2);
    dualfeas = normRd / (1 + norm(z, 2));
    arma::vec Axb = Ax - b;
    Rp = Axb - y;
    double normRp = norm(Rp, 2);
    primfeas = normRp / normb;

    double lasso = lam2 * group_l2_norm(Px, gs) + lam1 * sum(abs(x));
    double dualobj = -sum(square(y)) / 2 - dot(b, y);
    double primobj = sum(square(Axb)) / 2 + lasso;

    double gap = primobj - dualobj;
    relgap = std::abs(gap) / (1 + std::abs(primobj) + std::abs(dualobj));

    if (stopopt == 1) {
      stop = std::max(dualfeas, relgap) < stoptol;
    } else if (stopopt == 2) {
      eta_tmp = std::max(dualfeas, primfeas);
      if (eta_tmp < stoptol) {
        eta_1 = norm(x - proximal_combo(x + z, lam1, lam2, gs), 2) /
                (1 + norm(x, 2));
        stop = eta_1 < stoptol;
      }
    } else if (stopopt == 3) {
      stop = (dualfeas < stoptol) && (relgap < stoptol_gap);
    } else if (stopopt == 4) {
      stop = (dualfeas < stoptol) && (gap < stoptol * sum(square(b)));
    }

    ttime = time_diff(time_now(), tstart);
    primfeas_vec.push_back(primfeas);
    dualfeas_vec.push_back(dualfeas);
    sigma_vec.push_back(sig);
    primobj_vec.push_back(primobj);
    dualobj_vec.push_back(dualobj);
    gap_vec.push_back(gap);
    relgap_vec.push_back(relgap);
    ttime_vec.push_back(ttime);
    itersub_vec.push_back(as<int>(info_ncg["itersub"]));

    if (printyes) {
      Rcout << "\n"
            << std::right << std::setw(5) << iter << std::left
            << std::scientific << std::setprecision(2) << std::setw(3) << "|  ["
            << primfeas << " " << dualfeas << "]  " << relgap << " | "
            << std::setprecision(4) << primobj << " " << dualobj << " | "
            << std::setprecision(2) << std::setw(3) << ttime << " | " << sig
            << " |";
    }

    if (stop || (iter == maxit)) {
      std::string termination = "converged";
      if (iter == maxit)
        termination = "maxiter reached";
      runhist["termination"] = termination;
      runhist["iter"] = iter;
      runhist["nnz"] = cardcal(x);
      obj[0] = primobj;
      obj[1] = dualobj;
      final_iter = iter;
      break;
    }

    double ratio = primfeas / (dualfeas + datum::eps);
    ratio_seq_vec.push_back(ratio);

    if (ratio < 1) {
      prim_win++;
    } else {
      dual_win++;
    }

    int sigma_update_iter = sigma_fun(iter);
    double sigmascale = 1.0;
    double av_findstep = as<double>(runhist_ncg["av_findstep"]);
    if (primfeas > 100 * stoptol) {
      if (av_findstep > 5) {
        sigmascale = sqrt(3);
      } else {
        sigmascale = 3;
      }
    } else {
      if (av_findstep > 5) {
        sigmascale = sqrt(5);
      } else {
        sigmascale = 5;
      }
    }

    if ((iter % sigma_update_iter == 0) && breakyes < 0) {
      if (prim_win > std::max(1.0, 1.2 * dual_win)) {
        prim_win = 0;
        sig = std::min(sigmamax, sig * sigmascale);
      } else if (dual_win > std::max(1.0, 1.2 * prim_win)) {
        dual_win = 0;
        sig = std::max(sigmamin, sig / sigmascale);
      }
      if (breakyes >= 0 && iter >= 10) {
        sig = std::max(sigmamin, 2 * sig / sigmascale);
      }
    }
  }
  double eta =
      norm(x - proximal_combo(x + z, lam1, lam2, gs), 2) / (1 + norm(x, 2));
  runhist["primfeas"] = primfeas_vec;
  runhist["dualfeas"] = dualfeas_vec;
  runhist["sigma"] = sigma_vec;
  runhist["primobj"] = primobj_vec;
  runhist["dualobj"] = dualobj_vec;
  runhist["gap"] = gap_vec;
  runhist["relgap"] = relgap_vec;
  runhist["ttime"] = ttime_vec;
  runhist["ratio_seq"] = ratio_seq_vec;
  runhist["itersub"] = itersub_vec;

  List info;

  if (stopopt == 2) {
    double kktres = std::max(eta_tmp, eta_1);
    runhist["kktres"] = kktres;
    info["kktres"] = kktres;
  }

  info["maxfeas"] = maxfeas;
  info["iter"] = final_iter;
  info["ttime"] = ttime;
  info["termination"] = runhist["termination"];
  info["relgap"] = relgap;
  info["msg"] = runhist["termination"];
  info["eta"] = eta;

  return List::create(Named("y") = y, Named("z") = z, Named("x") = x,
                      Named("obj") = obj, Named("info") = info,
                      Named("runhist") = runhist);
}

// [[Rcpp::export]]
List sglssnal_main_interface(const arma::sp_mat &A, const arma::vec &b,
                             double lam1, double lam2, const List &gs_list,
                             const List &parmain, const arma::vec &y0,
                             const arma::vec &z0, const arma::vec &x0) {
  uvec G = as<uvec>(gs_list["G"]);
  mat ind = as<mat>(gs_list["ind"]);
  uint num_group = ind.n_cols;
  GroupStruct gs = {as<sp_mat>(gs_list["pma"]), G, ind, num_group};
  return sglssnal_main(A, b, lam1, lam2, gs, parmain, y0, z0, x0);
}