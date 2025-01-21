#ifndef GROUP_STRUCT_H
#define GROUP_STRUCT_H
#include <RcppArmadillo.h>

using namespace arma;

struct GroupStruct {
  sp_mat pma;
  uvec elem_ids;
  mat ind;
  uword num_group;

  // returns the subvector of x corresponding to the group with id group_id
  vec get_group_subview(const vec &x, int group_id) const {
    int kstart = ind(0, group_id);
    int kend = ind(1, group_id);
    return x(elem_ids.subvec(kstart, kend));
  }

  sp_mat get_group_subview(const sp_mat &x, int group_id) const {
    int kstart = ind(0, group_id);
    int kend = ind(1, group_id);
    return x.cols(elem_ids.subvec(kstart, kend));
  }
};

#endif