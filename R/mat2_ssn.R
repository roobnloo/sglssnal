mat2_ssn <- function(u, A, c, P, sig) {
  c2 <- c[2]
  v <- proximal_l1(u, c[1])
  Pv <- P$times(v)
  result <- P$ProjL2(Pv, c2)
  grp_nrm <- result$indicate
  I_grp <- which(grp_nrm != 0)
  supp_v <- v != 0
  r <- sum(supp_v)
  r2 <- length(I_grp)
  V2 <- NULL
  D <- NULL
  sp_dim <- 0
  id_yes <- 0

  if (!(r && r2)) {
    id_yes <- 1
    return(list(V2 = V2, D = D, sp_dim = sp_dim, id_yes = id_yes))
  }

  B <- Matrix::Matrix(0, nrow(A), r + r2, sparse = TRUE)
  C <- Matrix::Matrix(0, nrow(A), r2)
  s_start <- 1
  i <- 1

  for (k in I_grp) {
    G_k <- P$G[P$ind[1, k]:P$ind[2, k]]
    v_k <- v[G_k]
    indvk <- v_k != 0

    cw <- c2 * P$ind[3, k]
    par1 <- cw / grp_nrm[k]
    par1 <- sig * par1
    par2 <- par1 / grp_nrm[k]^2

    Al <- A[, G_k]
    Al <- Al[, indvk]

    Bl <- sqrt(sig - par1) * Al
    lenind1 <- ncol(Bl)
    s_end <- s_start + lenind1 - 1
    B[, s_start:s_end] <- as.numeric(Bl)
    s_start <- s_end + 1

    if (any(abs(v_k) > .Machine$double.eps)) {
      pv <- v_k[indvk]
      Asl <- Al %*% pv
      cl <- sqrt(par2) * Asl
      C[, i] <- cl
      i <- i + 1
    }
  }

  sp_dim <- s_start + r2 - 1
  B[, s_start:sp_dim] <- C
  D <- B[, 1:sp_dim]
  V2 <- crossprod(D)
  V2 <- as.matrix(V2)
  for (i in seq_len(sp_dim)) {
    V2[i, i] <- V2[i, i] + 1
  }

  return(list(V2 = V2, D = D, sp_dim = sp_dim, id_yes = id_yes))
}
