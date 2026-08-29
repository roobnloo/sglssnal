#' Largest eigenvalue of A A^T, used as the Lipschitz constant of the
#' least-squares gradient.
#' @description A single, unconditional call regardless of whether `A` is
#'   sparse or dense. `which = "LA"` and `which = "LM"` are equivalent here
#'   -- `tcrossprod(A)` is always positive semidefinite, so its largest
#'   eigenvalue by magnitude and its largest eigenvalue algebraically are
#'   the same value by construction -- and `"LA"` is not a valid `which`
#'   for `RSpectra::eigs()`'s plain-`matrix` (dense) method, so `"LM"`
#'   (the default) is the only choice that works for both. The R-version
#'   branch this replaced existed because `crossprod`/`tcrossprod` only
#'   became proper S3/S4 generics in R 4.4.0; DESCRIPTION now requires
#'   that as a minimum, so the pre-4.4.0 densify-everything fallback is
#'   dead code.
#' @param A n x p design matrix (dense or sparse).
#' @param n Number of rows of `A`.
#' @noRd
compute_lip <- function(A, n) {
  eigsopt <- list(retvec = FALSE)
  RSpectra::eigs(tcrossprod(A), k = 1, opts = eigsopt, n = n)$values
}
