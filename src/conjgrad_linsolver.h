#ifndef CONJGRAD_LINSOLVER_H
#define CONJGRAD_LINSOLVER_H

#include "group_struct.h"
#include <RcppArmadillo.h>

template <typename MatType>
Rcpp::List conjgrad_linsolver(const MatType &A, const arma::vec &rhs,
                              const arma::vec &u, double lam1, double lam2,
                              const GroupStruct &gs, int density, double sig);

// Explicit instantiations declaration
extern template Rcpp::List conjgrad_linsolver<arma::sp_mat>(
    const arma::sp_mat &A, const arma::vec &rhs, const arma::vec &u,
    double lam1, double lam2, const GroupStruct &gs, int density, double sig);

extern template Rcpp::List
conjgrad_linsolver<arma::mat>(const arma::mat &A, const arma::vec &rhs,
                              const arma::vec &u, double lam1, double lam2,
                              const GroupStruct &gs, int density, double sig);

#endif