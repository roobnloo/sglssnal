proximal_l2 <- function(z, c1, ind, num_group) {
  m <- length(z)
  Pz <- numeric(m)
  cw <- c1 * ind[3, 1:num_group]
  kstart <- ind[1, 1:num_group]
  kend <- ind[2, 1:num_group]

  for (j in 1:num_group) {
    nrm <- l2norm(z, kstart[j], kend[j])
    if (nrm > cw[j]) {
      Pz[kstart[j]:kend[j]] <- z[kstart[j]:kend[j]] * (1 - cw[j] / nrm)
    } else {
      Pz[kstart[j]:kend[j]] <- 0
    }
  }

  return(Pz)
}
