projection_l2 <- function(z, c1, ind, num_group) {
  m <- length(z)
  Pz <- numeric(m)
  indicate <- numeric(num_group)

  kstart <- ind[1, 1:num_group]
  kend <- ind[2, 1:num_group]
  cw <- c1 * ind[3, 1:num_group]

  # Calculate norms for each group
  nrms <- sapply(1:num_group, \(j) l2norm(z, kstart[j], kend[j]))

  # Update indicate vector
  indicate <- ifelse(nrms > cw, nrms, 0)

  # Apply the thresholding operation
  for (j in 1:num_group) {
    idx <- kstart[j]:kend[j]
    if (nrms[j] > cw[j]) {
      Pz[idx] <- z[idx] * (cw[j] / nrms[j])
    } else {
      Pz[idx] <- z[idx]
    }
  }

  return(list(Pz = Pz, indicate = indicate))
}
