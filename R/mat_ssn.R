mat_ssn <- function(u, A, c, P, sig) {
  v <- numeric(length(u))
  pma <- P$matrix
  proj2pv <- numeric(ncol(pma))
  n <- nrow(A)
  V <- matrix(0, n, n)

  PP <- P
  PP$G <- P$G - 1L
  PP$ind[1, ] <- P$ind[1, ] - 1L
  PP$ind[2, ] <- P$ind[2, ] - 1L
  mat_ssn_interface(
    as.numeric(u), A, c[1], c[2], sig, PP$matrix, PP$G, PP$ind,
    PP$num_group, V, proj2pv, v
  )

  return(list(V = V, Proj2_Pv = proj2pv, v = v))
}
