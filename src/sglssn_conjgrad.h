#ifndef SGLSSN_CONJGRAD_H
#define SGLSSN_CONJGRAD_H
#include "group_struct.h"
#include <RcppArmadillo.h>

template <typename MatType>
Rcpp::List sglssn_conjgrad(const arma::vec &y0, const arma::vec &Aty0,
                           const arma::vec &x0, const arma::vec &Ax0,
                           const MatType &A, const arma::vec &b, double lam1,
                           double lam2, const GroupStruct &gs, Rcpp::List &par,
                           bool printsub, int maxitersub, double tol);

// Explicit instantiations declaration
extern template Rcpp::List sglssn_conjgrad<arma::sp_mat>(
    const arma::vec &y0, const arma::vec &Aty0, const arma::vec &x0,
    const arma::vec &Ax0, const arma::sp_mat &A, const arma::vec &b,
    double lam1, double lam2, const GroupStruct &gs, Rcpp::List &par,
    bool printsub, int maxitersub, double tol);

extern template Rcpp::List
sglssn_conjgrad<arma::mat>(const arma::vec &y0, const arma::vec &Aty0,
                           const arma::vec &x0, const arma::vec &Ax0,
                           const arma::mat &A, const arma::vec &b, double lam1,
                           double lam2, const GroupStruct &gs, Rcpp::List &par,
                           bool printsub, int maxitersub, double tol);

#endif