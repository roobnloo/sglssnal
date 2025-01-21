l2norm <- function(z, kstart, kend) {
  nrm2 <- sum(z[kstart:kend]^2)
  nrm <- sqrt(nrm2)
  return(nrm)
}


group_l2_norm <- function(z, ind, num_group) {
  norms <- apply(ind, 2, \(col) {
    kstart <- col[1]
    kend <- col[2]
    weight <- col[3]
    weight * l2norm(z, kstart, kend)
  })

  return(sum(norms))
}
