#' Largest eigenvalue of A A^T, used as the Lipschitz constant of the
#' least-squares gradient.
#' @description A single, unconditional call regardless of whether `A` is
#'   sparse or dense. `which = "LA"` and `which = "LM"` are equivalent here
#'   -- `tcrossprod(A)` is always positive semidefinite, so its largest
#'   eigenvalue by magnitude and its largest eigenvalue algebraically are
#'   the same value by construction -- and `"LA"` is not a valid `which`
#'   for `RSpectra::eigs()`'s plain-`matrix` (dense) method, so `"LM"`
#'   (the default) is the only choice that works for both.
#'
#'   No `n` argument is passed to `eigs()` either, deliberately: it's only
#'   read by `eigs.function` (the matrix-free-operator method, for when
#'   `A` is a bare function rather than a concrete matrix -- never the
#'   case here). Checking `args(RSpectra:::eigs.matrix)` /
#'   `args(RSpectra:::eigs.dgCMatrix)` shows neither method even declares
#'   an `n` parameter; anything passed as `n` for a concrete-matrix `A`
#'   silently falls into `...` and is discarded, verified empirically
#'   (passing `n` correct, wildly wrong, or omitted all give the
#'   bit-identical eigenvalue for both dense and sparse input).
#' @param A n x p design matrix (dense or sparse).
#' @noRd
compute_lip <- function(A) {
  eigsopt <- list(retvec = FALSE)
  RSpectra::eigs(tcrossprod(A), k = 1, opts = eigsopt)$values
}
