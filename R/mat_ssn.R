mat_ssn <- function(u, A, c, P, sig) {
  c1 <- c[1]
  c2 <- c[2]
  v <- proximal_l1(u, c1)
  Pv <- P$times(v)
  n <- nrow(A)
  result <- P$ProjL2(Pv, c2)
  Proj2_Pv <- result$Pz
  grp_nrm <- result$indicate
  supp_v2 <- v != 0

  if (all(!supp_v2)) {
    V <- diag(n)
    return(list(V = V, Proj2_Pv = Proj2_Pv, v = v))
  }

  if (all(grp_nrm == 0)) {
    V <- diag(n)
    return(list(V = V, Proj2_Pv = Proj2_Pv, v = v))
  }

  V <- matrix(0, n, n)

  for (k in which(grp_nrm != 0)) {
    G_k <- P$G[P$ind[1, k]:P$ind[2, k]]
    v_k <- v[G_k]
    indvk <- v_k != 0

    cw <- c2 * P$ind[3, k]
    par1 <- cw / grp_nrm[k]
    par1 <- sig * par1
    par2 <- par1 / grp_nrm[k]^2

    Al <- A[, G_k]
    Al <- Al[, indvk]

    M1 <- tcrossprod(Al)
    V <- V + (sig - par1) * M1

    if (any(abs(v_k) > .Machine$double.eps)) {
      pv <- v_k[indvk]
      Asl <- Al %*% pv
      M2 <- tcrossprod(Asl)
      V <- V + par2 * M2
    }
  }

  for (i in 1:n) {
    V[i, i] <- V[i, i] + 1
  }

  return(list(V = V, Proj2_Pv = Proj2_Pv, v = v))
}
