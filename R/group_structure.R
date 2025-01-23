group_structure <- function(n, G, ind, pfgroup) {
  len_group <- (ind[2, ] - ind[1, ]) + 1

  gs <- list()
  gs$len_group <- len_group
  gs$ntotal <- sum(len_group)


  I <- 1:gs$ntotal
  J <- G
  V <- rep(1, gs$ntotal)
  gs$pma <- Matrix::sparseMatrix(i = I, j = J, x = V)

  # This code adjusts for zero-based indexing in C++
  G <- G - 1L
  ind[1, ] <- ind[1, ] - 1L
  ind[2, ] <- ind[2, ] - 1L

  indw <- rbind(ind, pfgroup)

  gs$ind <- indw
  gs$G <- G

  return(gs)
}
