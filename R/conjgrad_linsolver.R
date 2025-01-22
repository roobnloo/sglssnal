conjgrad_linsolver <- function(A, rhs, u, c, P, par) {
  PP <- P
  PP$G <- P$G - 1L
  PP$ind[1, ] <- P$ind[1, ] - 1L
  PP$ind[2, ] <- P$ind[2, ] - 1L

  result <- conjgrad_linsolver_impl(
    A, rhs, as.numeric(u), c[1], c[2], PP$matrix, PP$G, PP$ind,
    PP$num_group, par$nnz, par$sigma
  )

  return(result)
}
