#' Cross-validation for sglssnal
#' @description Perform cross-validation for the sglssnal algorithm over
#'   a grid of lambda and alpha values.
#' @inheritParams sglssnal
#' @param alphas Vector of alpha parameters to tune. Default is `0.75`.
#' @param lambdas Vector of lambda parameters to tune. If `NULL`, a
#'   path is automatically generated based on `norm(t(A) %*% b, "I")` and
#'   `lambda_min_ratio`. If supplied, `nlambda` and `lambda_min_ratio` are ignored_
#'   The values are sorted in decreasing order.
#' @param nlambda Number of lambda values to use.
#' @param lambda_min_ratio Minimum ratio of the smallest to largest lambda.
#' @param nfolds Number of folds for cross-validation.
#'   Ignored if `foldid` is provided.
#' @param foldid Vector of integers specifying the fold for each observation.
#'   If `NULL`, a default assignment is used.
#' @param stoptolcv Tolerance for convergence in cross-validation folds.
#'   May be set to a smaller value than `stoptol` for faster convergence.
#'   Default is `1e-4`.
#' @param quietall Suppress all messages.
#' @param ... Additional arguments passed to [sglssnal()].
#' @return List of class `cv.sglssanl, sglssnal`, containing the following components:
#' * `obj`: Vector containing primal and dual objective values at the solution.
#' * `x`: The primal variable of interest.
#' * `y`: The y dual variable.
#' * `z`: The z dual variable.
#' * `info`: List containing information about the optimization process such as
#'   relative gap, iteration number, and run time.
#' * `cv_info`: List containing the tuning grid, cross-validation errors, and
#'   the inidices of the optimal lambda and alpha.
#' @export
cv.sglssnal <- function(
    A, b, grp_vec, grp_idx, alphas = 0.75, lambdas = NULL,
    nlambda = 100, lambda_min_ratio = 1e-3,
    nfolds = 5, foldid = NULL, printyes = TRUE,
    stoptol = 1e-6, stoptolcv = 1e-4, quietall = FALSE, ...) {
  nalpha <- length(alphas)
  n <- length(b)
  p <- ncol(A)
  if (quietall) {
    printyes <- FALSE
  }
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
  if (is.null(lambdas)) {
    if (nlambda < 0) {
      stop("nlambda must be a positive integer")
    }
    alpha_min <- min(alphas)
    if (alpha_min == 0) {
      alpha_min <- 1e-2
    }
    # This choice of lam_max ensures a sparse lasso solution exists on the path across all alphas.
    lam_max <- norm(crossprod(A, b), "I") / alpha_min
    lambdas <- lam_max *
      exp(seq(log(1), log(lambda_min_ratio), length.out = nlambda))
  } else {
    lambdas <- sort(lambdas, decreasing = TRUE)
    nlambda <- length(lambdas)
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
    stopifnot(
      "foldid must be a contiguous vector starting from 1" =
        all(foldid >= 1) && all(diff(sort(unique(foldid))) == 1)
    )
    nfolds <- max(foldid)
  }

  if (!quietall) {
    message(sprintf(
      "Cross-validating sglssnal over a (%dx%d) grid with %d folds...",
      nlambda, nalpha, nfolds
    ))
  }
  tstart <- Sys.time()
  for (t in 1:nfolds) {
    tstart_fold <- Sys.time()
    foldidx1 <- foldid != t
    Atrain <- A[foldidx1, , drop = FALSE]
    btrain <- b[foldidx1]
    foldidx2 <- foldid == t
    Atest <- A[foldidx2, , drop = FALSE]
    btest <- b[foldidx2]

    eigsopt <- list(retvec = FALSE)
    Lip <- NULL
    if (getRversion() < as.numeric_version("4.4.0")) {
      Lip <- RSpectra::eigs(
        tcrossprod(as.matrix(Atrain)),
        k = 1, opts = eigsopt, n = n
      )$values
    } else {
      Lip <- RSpectra::eigs(
        tcrossprod(Atrain),
        k = 1, which = "LA", opts = eigsopt, n = n
      )$values
    }

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
    if (!quietall) {
      message(
        sprintf(
          "Fold %d - %3.2fs", t,
          as.numeric(difftime(Sys.time(), tstart_fold, units = "secs"))
        )
      )
    }
  }
  if (!quietall) {
    message(
      sprintf(
        "Total time: %3.2fs",
        as.numeric(difftime(Sys.time(), tstart, units = "secs"))
      )
    )
  }

  min_error <- which(cv_err == min(cv_err), arr.ind = TRUE, useNames = FALSE)
  min_lambda_id <- min_error[1, 1]
  min_alpha_id <- min_error[1, 2]
  if (!quietall) {
    message(sprintf(
      "Minimum error at (lambda, alpha) index (%d, %d)",
      min_lambda_id, min_alpha_id
    ))
  }

  result <- sglssnal(
    A, b, grp_vec, grp_idx, lambdas[min_lambda_id], alphas[min_alpha_id],
    printyes = printyes, stoptol = stoptol, ...
  )

  cv_info <- list(
    alphas = alphas,
    lambdas = lambdas,
    cvm = cv_err,
    cv_idx = c("cv_lambda_id" = min_lambda_id, "cv_alpha_id" = min_alpha_id)
  )
  result$cv_info <- cv_info
  class(result) <- c("cv.sglssnal", "sglssnal")
  return(result)
}
