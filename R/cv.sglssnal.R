#' Cross-validation for sglssnal
#' @inheritParams sglssnal
#' @param lambdas Vector of lambda parameters to tune.
#' @param alphas Vector of alpha parameters to tune.
#' @param nfolds Number of folds for cross-validation.
#'   Ignored if `foldid` is provided.
#' @param foldid Vector of integers specifying the fold for each observation.
#'   If `NULL`, a default assignment is used.
#' @param stoptolcv Tolerance for convergence in cross-validation folds.
#'   May be set to a smaller value than `stoptol` for faster convergence.
#'   Default is `1e-4`.
#' @param ... Additional arguments passed to [sglssnal()].
#' @export
cv.sglssnal <- function(
    A, b, grp_vec, grp_idx, lambdas, alphas,
    nfolds = 5, foldid = NULL, printyes = TRUE,
    stoptol = 1e-6, stoptolcv = 1e-4, ...) {
  nlambda <- length(lambdas)
  nalpha <- length(alphas)
  n <- length(b)
  p <- ncol(A)
  stopifnot("nrow(A) must be equal to length(b)" = nrow(A) == length(b))
  if (!is.null(foldid)) {
    stopifnot(
      "length(foldid) must be equal to nrow(A)" = length(foldid) == nrow(A)
    )
    stopifnot(
      "foldid must be a contiguous vector starting from 1" =
        all(foldid >= 1) && all(diff(sort(unique(foldid))) == 1)
    )
  }
  A <- Matrix::Matrix(A, sparse = TRUE)

  cv_err <- matrix(0, nlambda, nalpha)

  if (is.null(foldid)) {
    fold_sizes <- rep(floor(n / nfolds), nfolds)
    fold_sizes[nfolds] <- fold_sizes[nfolds] + (n %% nfolds)
    foldid <- rep(1:nfolds, fold_sizes)
  } else {
    stopifnot("foldid must be a vector of integers" = is.integer(foldid))
    stopifnot("length(foldid) must equal nrow(A)" = length(foldid) == nrow(A))
    stopifnot()
    nfolds <- max(foldid)
  }

  message("Cross-validation with ", nfolds, " folds")
  tstart <- Sys.time()
  for (t in 1:nfolds) {
    tstart_fold <- Sys.time()
    foldidx1 <- foldid != t
    Atrain <- A[foldidx1, , drop = FALSE]
    btrain <- b[foldidx1]
    foldidx2 <- foldid == t
    Atest <- A[foldidx2, , drop = FALSE]
    btest <- b[foldidx2]

    eigsopt <- list(issym = TRUE)
    Lip <- RSpectra::eigs(
      tcrossprod(Atrain),
      k = 1, which = "LA", opts = eigsopt, n = n
    )$values

    for (k in 1:nalpha) {
      y0 <- rep(0, length(btrain))
      z0 <- rep(0, p)
      x0 <- rep(0, p)

      for (i in 1:nlambda) {
        lambda <- lambdas[i]
        alpha <- alphas[k]
        result <- sglssnal(
          Atrain, btrain, grp_vec, grp_idx, lambda, alpha,
          Lip = Lip, y0 = y0, z0 = z0, x0 = x0, printyes = FALSE,
          stoptol = stoptolcv,
          ...
        )
        y0 <- result$y
        z0 <- result$z
        x0 <- result$x
        res <- Atest %*% x0 - btest
        cv_err[i, k] <- cv_err[i, k] + sum(res^2)
      }
    }
    message(
      sprintf(
        "Fold %d - %3.2fs", t,
        as.numeric(difftime(Sys.time(), tstart_fold, units = "secs"))
      )
    )
  }
  message(
    sprintf(
      "Total time: %3.2fs",
      as.numeric(difftime(Sys.time(), tstart, units = "secs"))
    )
  )

  min_error <- which(cv_err == min(cv_err), arr.ind = TRUE)
  min_lambda_id <- min_error[1, 1]
  min_alpha_id <- min_error[1, 2]
  message(sprintf(
    "Minimum error at (lambda, alpha) index (%d, %d)",
    min_lambda_id, min_alpha_id
  ))

  result <- sglssnal(
    A, b, grp_vec, grp_idx, lambdas[min_lambda_id], alphas[min_alpha_id],
    printyes = printyes, stoptol = stoptol, ...
  )

  cv_info <- list(
    cvm = cv_err,
    cv_idx = c("cv_lambda_id" = min_lambda_id, "cv_alpha_id" = min_alpha_id)
  )
  result$cv_info <- cv_info
  return(result)
}
