#' Run Sparse-Group Lasso via Semismooth Newton Augmented Lagrangian
#' @description Fits a sparse-group lasso model using
#'   second-order information to solve the dual problem.
#'   The penalty function is given by
#'   \deqn{\Phi(x) = \lambda_1 ||x||_1 + \lambda_2 \sum_{i=1}^g w_i ||x_{G_i}||_2}
#'   where \eqn{G_i} is the \eqn{i}-th group and \eqn{w_i} its penalty factor.
#'   The primal problem is given by
#'   \deqn{\min_{x \in \mathbb{R}^p}\; \frac{1}{2} ||Ax - b||_2^2 + \Phi(x)}
#'   while the dual problem has the form
#'   \deqn{\max_{y \in \mathbb{R}^n, z \in \mathbb{R}^p}\; -\langle b, y\rangle - \frac{1}{2}||y||_2^2 - \Phi^\ast(z) \text{s.t.} A^\top y + z = 0}
#'   where \eqn{\Phi^\ast(z)} denotes the Fenchel conjugate of \eqn{\Phi}.
#'   The algorithm is based on the work of Zhang et al. (2020).
#' @param A \eqn{n \times p} design matrix.
#' @param b \eqn{n} response vector.
#' @param grp_vec Vector of indicies of variables in each group.
#'   If there are \eqn{g} groups and `G_i` contains the indices of the
#'   `i`-th group, then `grp_vec` should be the concatenated vector
#'   `c(G_1, G_2, ..., G_g)`.
#' @param grp_idx \eqn{2 \times g} matrix indexing the groups in `grp_vec`.
#'   `grp_vec[grp_idx[1, i]:grp_idx[2, i]]` are the indices of the `i`-th
#'   group.
#' @param lambda Vector of penalty parameters or a single value. If `NULL`, a
#'   path is automatically generated based on `norm(t(A) %*% b, "I")` and
#'   `lambda_min_ratio`. The path is fit using warm starts from larger to smaller
#'   lambda values.
#' @param alpha Determines the relative weight of the \eqn{\ell_1} and
#'   group \eqn{\ell_2} penalty. Must be in \eqn{[0, 1]}.
#' @param nlambda Number of lambda values to use when `lambda` is `NULL`. Default is 100.
#' @param lambda_min_ratio Minimum ratio of the smallest to largest lambda when `lambda` is `NULL`.
#'   Default is 1e-4.
#' @param pfgroup Penalty factor for each group in the group lasso.
#'   Default is a vector of ones, indicating no weighting for any group.
#' @param standardize Whether to standardize the columns of A to have unit norm.
#'   Default is `TRUE`.
#' @param intercept Centers mean function at 0 through \eqn{y - \bar{y}}.
#'   Default is `TRUE`.
#' @param stoptol Tolerance for stopping criteria. Default is `1e-6`.
#' @param stopopt Stopping criteria. 1: relative duality gap and feasibility,
#'   2: KKT conditions, 3: dual feasibility and relative duality gap,
#'   4: dual feasibility and absolute duality gap. Default is `2L`.
#' @param printyes Print progress in main loop. Default is `TRUE`.
#' @param printsub Print progress in subproblem. Default is `FALSE`.
#' @param maxit Maximum number of iterations. Default is `5000L`.
#' @param Lip Lipschitz constant for the step size.
#'   Automatically computed as the maximum eigenvalue of AA' if `NULL`.
#' @param y0 optional initialization vector for dual variable `y`
#' @param z0 optional initialization vector for dual variable `z`
#' @param x0 optional initialization vector for primal variable `x`
#' @return List of class `sglssnal` containing the following components:
#' * `obj`: Matrix containing primal and dual objective values for each lambda.
#' * `x`: Matrix of primal variables, with each column corresponding to a lambda.
#' * `x0`: Vector of intercept values for each lambda.
#' * `y`: Matrix of y dual variables, with each column corresponding to a lambda.
#' * `z`: Matrix of z dual variables, with each column corresponding to a lambda.
#' * `info`: List containing information about the optimization process for each lambda.
#' * `lambda`: Vector of lambda values used.
#' @importFrom Rcpp sourceCpp
#' @useDynLib sglssnal
#' @references Zhang, Y., Zhang, N., Sun, D., & Toh, K. C. (2020).
#'   \emph{An efficient Hessian based algorithm for solving large-scale
#'   sparse group Lasso problems.} Mathematical Programming, 179, 223-263.
#'   \doi{https://doi.org/10.1007/s10107-018-1329-6}.
#' @export
sglssnal <- function(
    A, b, grp_vec, grp_idx, lambda = NULL,
    nlambda = 100, lambda_min_ratio = 1e-4, alpha = 0.05,
    pfgroup = rep(1, ncol(grp_idx)), intercept = TRUE,
    standardize = TRUE, stoptol = 1e-6, stopopt = 2L,
    printyes = TRUE, printsub = FALSE, maxit = 5000L, Lip = NULL,
    y0 = NULL, z0 = NULL, x0 = NULL) {
  stopifnot("length(alpha) must be 1" = length(alpha) == 1)
  # Generate lambda sequence if lambda is NULL
  if (is.null(lambda)) {
    if (nlambda <= 0) {
      stop("nlambda must be a positive integer")
    }
    alpha_min <- alpha
    if (alpha_min == 0) {
      alpha_min <- 1e-2
    }
    # This choice of lam_max ensures a sparse lasso solution exists on the path
    lam_max <- norm(crossprod(A, b), "I") / alpha_min
    lambda <- lam_max *
      exp(seq(log(1), log(lambda_min_ratio), length.out = nlambda))
  }

  # Convert lambda to a vector if it's a single value
  if (length(lambda) == 1) {
    lambda <- c(lambda)
  }

  # Sort lambda in descending order for the regularization path
  lambda <- sort(lambda, decreasing = TRUE)
  nlambda <- length(lambda)

  stopifnot("lambda values must be positive" = all(lambda > 0))
  stopifnot("alpha must be in [0, 1]" = alpha >= 0 & alpha <= 1)
  stopifnot("nrow(A) must be equal to length(b)" = nrow(A) == length(b))
  stopifnot("length(pfgroup) must be equal to ncol(grp_idx)" = length(pfgroup) == ncol(grp_idx))
  stopifnot("stopopt must be one of 1, 2, 3, or 4" = stopopt %in% c(1L, 2L, 3L, 4L))
  stopifnot("maxit must be a positive integer" = maxit > 0)
  stopifnot("stoptol must be a positive number" = stoptol > 0)
  stopifnot("grp_idx must be within the range of grp_vec" = max(grp_idx) <= max(grp_vec))

  n <- length(b)
  p <- ncol(A)

  # Prepare matrices for storing results
  x_mat <- matrix(0, nrow = p, ncol = nlambda)
  x0_vec <- numeric(nlambda)
  y_mat <- matrix(0, nrow = n, ncol = nlambda)
  z_mat <- matrix(0, nrow = p, ncol = nlambda)

  # Preprocessing steps
  b_mean <- mean(b)
  A_sd <- sqrt(Matrix::colSums(A^2))
  A_sd[A_sd < sqrt(.Machine$double.eps)] <- 1

  if (intercept) {
    b <- b - b_mean
  }

  if (standardize) {
    if (inherits(A, "sparseMatrix")) {
      A <- A %*% Matrix::Diagonal(x = 1 / A_sd)
    } else {
      A <- sweep(A, 2, A_sd, "/")
    }
  }

  if (is.null(Lip)) {
    eigsopt <- list(retvec = FALSE)
    tstartLip <- Sys.time()
    if (getRversion() < as.numeric_version("4.4.0")) {
      Lip <- RSpectra::eigs(
        tcrossprod(as.matrix(A)),
        k = 1, opts = eigsopt, n = n
      )$values
    } else if (inherits(A, "sparseMatrix")) {
      Lip <- RSpectra::eigs(
        tcrossprod(A),
        k = 1, which = "LA", opts = eigsopt, n = n
      )$values
    } else {
      Lip <- RSpectra::eigs(
        tcrossprod(A),
        k = 1, which = "LM", opts = eigsopt, n = n
      )$values
    }

    if (printyes) {
      message(sprintf(
        "\n Lip = %3.2e, time = %3.2f",
        Lip, as.numeric(difftime(Sys.time(), tstartLip, units = "secs"))
      ))
    }
  }

  # Initialize warm start variables
  y <- rep(0, n)
  z <- rep(0, p)
  x <- rep(0, p)

  if (!is.null(y0) && !is.null(z0) && !is.null(x0)) {
    if (length(y0) == n && length(z0) == p && length(x0) == p) {
      y <- y0
      z <- z0
      x <- x0
    } else {
      stop(sprintf(
        "y0, z0, and x0 must have dimensions %d, %d, and %d respectively.",
        n, p, p
      ))
    }
  }

  parmain <- list(
    stoptol = stoptol,
    Lip = Lip,
    maxit = maxit,
    stopopt = stopopt,
    printyes = printyes,
    printsub = printsub,
    p = p,
    n = n
  )

  gs <- group_structure(p, grp_vec, grp_idx, pfgroup)

  # Store info for all lambda values
  all_info <- list()
  all_obj <- matrix(0,
    nrow = 2, ncol = nlambda,
    dimnames = list(c("primal objective", "dual objective"), NULL)
  )

  # Loop through lambda values
  for (i in 1:nlambda) {
    if (printyes) {
      message(sprintf("\nFitting model for lambda = %g (%d/%d)", lambda[i], i, nlambda))
    }

    lambda1 <- alpha * lambda[i]
    lambda2 <- lambda[i] * (1 - alpha)

    result <- NULL
    if (inherits(A, "sparseMatrix")) {
      result <- sglssnal_main_interface(
        A, b, lambda1, lambda2, gs, parmain, y, z, x, intercept
      )
    } else {
      result <- sglssnal_main_interface_dense(
        A, b, lambda1, lambda2, gs, parmain, y, z, x, intercept
      )
    }

    # Get results for current lambda and use as warm start for next lambda
    y <- result$y
    z <- result$z
    x <- result$x

    # Apply hard threshold for numerical stability
    x[abs(x) <= .Machine$double.eps] <- 0

    # Store results for current lambda
    if (standardize) {
      x_mat[, i] <- x / A_sd
    } else {
      x_mat[, i] <- x
    }

    y_mat[, i] <- y
    z_mat[, i] <- z

    # Calculate intercept
    if (intercept) {
      x0_vec[i] <- b_mean + result$intercept
    }

    # Store optimization info
    info_main <- result$info
    runhist_main <- result$runhist

    iter <- info_main$iter
    msg <- info_main$msg
    if (iter == maxit) msg <- "maximum iteration reached"
    dualfeas <- runhist_main$dualfeas[iter]
    primfeas <- runhist_main$primfeas[iter]
    maxfeas <- max(primfeas, dualfeas)

    dualobj <- runhist_main$dualobj[iter]
    primobj <- runhist_main$primobj[iter]
    all_obj[1, i] <- primobj
    all_obj[2, i] <- dualobj

    all_info[[i]] <- list(
      relgap = info_main$relgap,
      iter = iter,
      time_seconds = info_main$ttime,
      eta = info_main$eta,
      maxfeas = maxfeas,
      nnz = runhist_main$nnz,
      msg = msg,
      mse = mean((as.numeric(b - A %*% x))^2)
    )

    if (printyes) {
      message("\n****************************************")
      message(sprintf(" SSNAL       : %s", msg))
      message(sprintf(" iteration   : %d", iter))
      message(sprintf(" time(s)     : %3.2f", info_main$ttime))
      message(sprintf(" prim_obj    : %4.8e", primobj))
      message(sprintf(" dual_obj    : %4.8e", dualobj))
      message(sprintf(" relgap      : %4.5e", info_main$relgap))
      message(sprintf(" primfeas    : %3.2e", primfeas))
      message(sprintf(" dualfeas    : %3.2e", dualfeas))
      message(sprintf(" eta         : %3.2e", info_main$eta))
      message(sprintf(" nnz         : %d", runhist_main$nnz))
    }
  }

  stepnames <- paste("s", seq_len(nlambda) - 1, sep = "")
  vnames <- colnames(A)
  if (is.null(vnames)) {
    vnames <- paste("V", seq_len(p), sep = "")
  }
  dimnames(x_mat) <- list(vnames, stepnames)
  x <- Matrix::drop0(x_mat)

  fit <- list(
    lambda = lambda,
    alpha = alpha,
    obj = all_obj,
    x0 = x0_vec,
    x = x,
    y = y_mat,
    z = z_mat,
    intercept = intercept,
    info = all_info
  )
  class(fit) <- "sglssnal"
  return(fit)
}
