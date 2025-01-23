#ifndef NORM_OPS_H
#define NORM_OPS_H
#include "group_struct.h"
#include <RcppArmadillo.h>

using namespace arma;

arma::vec proximal_l1(const arma::vec &w, double lam);
arma::vec proximal_l2(const arma::vec &z, double lam, const GroupStruct &gs);
arma::vec proximal_combo(const arma::vec &v, double lam1, double lam2,
                         const GroupStruct &gs);
arma::vec projection_l2(const arma::vec &z, double lam, const GroupStruct &gs,
                        arma::vec &grp_norms);
double group_l2_norm(const arma::vec &z, const GroupStruct &gs);
int cardcal(const arma::vec &x, double r = 0.999);

#endif