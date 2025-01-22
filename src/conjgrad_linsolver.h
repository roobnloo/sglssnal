#ifndef CONJGRAD_LINSOLVER_H
#define CONJGRAD_LINSOLVER_H

#include "group_struct.h"
#include <RcppArmadillo.h>

Rcpp::List conjgrad_linsolver_impl(const arma::sp_mat &A, const arma::vec &rhs,
                                   const arma::vec &u, double lam1, double lam2,
                                   const GroupStruct &gs, int density,
                                   double sig);

#endif