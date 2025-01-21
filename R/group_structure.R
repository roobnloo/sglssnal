group_structure_matrix <- function(i, G, ind, n, len_group) {
  tmp <- len_group[i]
  I <- 1:tmp
  J <- G[ind[1, i]:ind[2, i]]
  V <- rep(1, tmp)
  Pi <- Matrix::sparseMatrix(i = I, j = J, x = V, dims = c(tmp, n))
  return(Pi)
}

group_structure <- function(n, G, ind) {
  num_group <- ncol(ind)
  len_group <- (ind[2, ] - ind[1, ]) + 1

  P <- list()
  P$num_group <- num_group
  P$len_group <- len_group
  P$ntotal <- sum(len_group)
  P$ind <- ind
  P$G <- G

  I <- 1:P$ntotal
  J <- G
  V <- rep(1, P$ntotal)
  pma <- Matrix::sparseMatrix(i = I, j = J, x = V)
  P$matrix <- pma

  P$Pi <- \(i) group_structure_matrix(i, G, ind, n, len_group)
  P$times <- \(x) pma %*% x
  P$trans <- \(x) crossprod(pma, x)

  P$ProjL2 <- \(z, c1) projection_l2(z, c1, ind, num_group)
  P$ProxL2 <- \(z, c1) proximal_l2(z, c1, ind, num_group)
  P$Lasso_fx <- \(x) group_l2_norm(pma %*% x, ind, num_group)
  P$Lasso_fz <- \(z) group_l2_norm(z, ind, num_group)

  return(P)
}
