#' Build the internal group structure from a membership vector
#' @description Converts a user-facing length-`n` group-membership vector
#'   into the group-major representation the C++ solvers expect (a sparse
#'   permutation matrix `pma`, a 0-based concatenated index vector `G`, and a
#'   3-row `ind` matrix of 0-based start/end positions plus per-group
#'   weight). `pfgroup` is indexed by `sort(unique(group))` order.
#' @param n Number of variables (`length(group)`, i.e. `ncol(A)`).
#' @param group Length-`n` group-membership vector.
#' @param pfgroup Penalty factor for each group, ordered by
#'   `sort(unique(group))`.
#' @noRd
group_structure <- function(n, group, pfgroup) {
  stopifnot("length(group) must equal n" = length(group) == n)
  stopifnot("group must not contain NA" = !anyNA(group))
  ug <- sort(unique(group))
  g <- length(ug)
  stopifnot(
    "length(pfgroup) must equal the number of unique groups" =
      length(pfgroup) == g
  )

  # Integer group code per coordinate (1..g). Sort by code, breaking ties by
  # original position, so each group's member order is preserved regardless
  # of which order() method gets selected internally.
  codes <- match(group, ug)
  ord <- order(codes, seq_along(codes))
  sizes <- tabulate(codes, nbins = g)

  gs <- list()
  gs$len_group <- sizes
  gs$ntotal <- sum(sizes)

  I <- seq_len(gs$ntotal)
  J <- ord
  V <- rep(1, gs$ntotal)
  gs$pma <- Matrix::sparseMatrix(i = I, j = J, x = V)

  end <- cumsum(sizes)
  start <- c(1L, end[-g] + 1L)

  # This code adjusts for zero-based indexing in C++
  gs$G <- ord - 1L
  gs$ind <- rbind(start - 1L, end - 1L, pfgroup)

  return(gs)
}
