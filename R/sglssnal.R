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
#' @param lambda Penalty parameter together with `alpha`.
#'   `alpha * lambda` is the \eqn{\ell_1} penalty and `(1 - alpha) * lambda`
#'   is the group \eqn{\ell_2} penalty.
#' @param alpha Determines the relative weight of the \eqn{\ell_1} and
#'   group \eqn{\ell_2} penalty. Must be in \eqn{[0, 1]}.
#' @param pfgroup Penalty factor for each group in the group lasso.
#'   Default is a vector of ones, indicating no weighting for any group.
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
#' * `obj`: Vector containing primal and dual objective values at the solution.
#' * `x`: The primal variable of interest.
#' * `y`: The y dual variable.
#' * `z`: The z dual variable.
#' * `info`: List containing information about the optimization process such as
#'   relative gap, iteration number, and run time.
#' @importFrom Rcpp sourceCpp
#' @useDynLib sglssnal
#' @references Zhang, Y., Zhang, N., Sun, D., & Toh, K. C. (2020).
#'   \emph{An efficient Hessian based algorithm for solving large-scale
#'   sparse group Lasso problems.} Mathematical Programming, 179, 223-263.
#'   \doi{https://doi.org/10.1007/s10107-018-1329-6}.
#' @export
sglssnal <- function(
    A, b, grp_vec, grp_idx, lambda, alpha,
    pfgroup = rep(1, ncol(grp_idx)), stoptol = 1e-6, stopopt = 2L,
    printyes = TRUE, printsub = FALSE, maxit = 5000L, Lip = NULL,
    y0 = NULL, z0 = NULL, x0 = NULL) {
  stopifnot("lambda must be positive" = lambda > 0)
  stopifnot("alpha must be in [0, 1]" = alpha >= 0 & alpha <= 1)
  stopifnot("nrow(A) must be equal to length(b)" = nrow(A) == length(b))
  stopifnot("length(pfgroup) must be equal to ncol(grp_idx)" = length(pfgroup) == ncol(grp_idx))
  stopifnot("stopopt must be one of 1, 2, 3, or 4" = stopopt %in% c(1L, 2L, 3L, 4L))
  stopifnot("maxit must be a positive integer" = maxit > 0)
  stopifnot("stoptol must be a positive number" = stoptol > 0)

  lambda1 <- alpha * lambda
  lambda2 <- lambda * (1 - alpha)

  A <- Matrix::Matrix(A, sparse = TRUE)
  n <- length(b)
  p <- ncol(A)

  if (is.null(Lip)) {
    eigsopt <- list(retvec = FALSE)
    tstartLip <- Sys.time()
    Lip <- NULL
    if (getRversion() < as.numeric_version("4.4.0")) {
      Lip <- RSpectra::eigs(
        tcrossprod(as.matrix(A)),
        k = 1, opts = eigsopt, n = n
      )$values
    } else {
      Lip <- RSpectra::eigs(
        tcrossprod(A),
        k = 1, which = "LA", opts = eigsopt, n = n
      )$values
    }
    if (printyes) {
      message(sprintf(
        "\n Lip = %3.2e, time = %3.2f",
        Lip, as.numeric(difftime(Sys.time(), tstartLip, units = "secs"))
      ))
    }
  }

  y <- rep(0, n)
  z <- rep(0, p)
  x <- z
  if (!is.null(y0) && !is.null(z0) && !is.null(x0)) {
    if (length(y0) == n && length(z0) == p && length(x0) == p) {
      y <- y0
      z <- z0
      x <- x0
    } else {
      stop(sprintf(
        "y0, z0, and x0 must have dimensions &d, &d, and &d respectively.",
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
  result <- sglssnal_main_interface(
    A, b, lambda1, lambda2, gs, parmain, y, z, x
  )
  y <- result$y
  z <- result$z
  x <- result$x
  x[abs(x) <= .Machine$double.eps] <- 0 # hard threshold for numerical stability
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
  obj <- c("primal objective" = primobj, "dual objective" = dualobj)

  info <- list(
    relgap = info_main$relgap,
    iter = iter,
    time_seconds = info_main$ttime,
    eta = info_main$eta,
    maxfeas = maxfeas,
    nnz = runhist_main$nnz,
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

  fit <- list(
    obj = obj, x = x, y = y, z = z, info = info
  )
  class(fit) <- "sglssnal"
  return(fit)
}
