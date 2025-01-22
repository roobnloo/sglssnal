mat2_ssn <- function(u, A, c, P, sig) {
  PP <- P
  PP$G <- P$G - 1L
  PP$ind[1, ] <- P$ind[1, ] - 1L
  PP$ind[2, ] <- P$ind[2, ] - 1L
  mat2_result <- mat2_ssn_interface(
    as.numeric(u), A, c[1], c[2], sig, PP$matrix, PP$G, PP$ind,
    PP$num_group
  )

  return(mat2_result)
}
