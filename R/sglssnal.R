#' Run Sparse-Group Lasso via Semismooth Newton Augmented Lagrangian
#' @description Fits a sparse-group lasso model using
#'   second-order information to solve the dual problem.
#' @param Ainput \eqn{n \times p} design matrix.
#' @param b \eqn{n} response vector.
#' @param lambda1 The \eqn{\ell_1} penalty.
#' @param lambda2 The group-wise \eqn{\ell_2} penalty.
#' @param grp_vec Vector of indicies of variables in each group.
#'   If there are \eqn{g} groups and `G_i` contains the indices of the
#'   `i`-th group, then `grp_vec` should be the concatenated vector
#'   `c(G_1, G_2, ..., G_g)`.
#' @param grp_idx \eqn{2 \times g} matrix indexing the groups in `grp_vec`.
#'   `grp_vec[grp_idx[1, i]:grp_idx[2, i]]` are the indices of the `i`-th
#'   group.
#' @param pfgroup Penalty factor for each group in the group lasso.
#'   Default is a vector of ones, indicating no weighting.
#' @param stoptol Tolerance for stopping criteria. Default is `1e-6`.
#' @param stopopt Stopping criteria. 1: relative duality gap and feasibility,
#'   2: KKT conditions, 3: dual feasibility and relative duality gap,
#'   4: dual feasibility and absolute duality gap. Default is `2L`.
#' @param printyes Print progress in main loop. Default is `TRUE`.
#' @param printsub Print progress in subproblem. Default is `FALSE`.
#' @param maxit Maximum number of iterations. Default is `5000L`.
#' @param y0 optional initialization vector
#' @param z0 optional initialization vector
#' @param x0 optional initialization vector
#' @return List containing the following components:
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
    Ainput, b, lambda1, lambda2, grp_vec, grp_idx,
    pfgroup = rep(1, ncol(grp_idx)), stoptol = 1e-6, stopopt = 2L,
    printyes = TRUE, printsub = FALSE, maxit = 5000L,
    y0 = NULL, z0 = NULL, x0 = NULL) {
  stopifnot("lambda1 and lambda2 must be nonnegative" = lambda1 >= 0 && lambda2 >= 0)
  stopifnot("nrow(A) must be equal to length(b)" = nrow(Ainput) == length(b))
  stopifnot("length(pfgroup) must be equal to ncol(grp_idx)" = length(pfgroup) == ncol(grp_idx))
  stopifnot("stopopt must be one of 1, 2, 3, or 4" = stopopt %in% c(1L, 2L, 3L, 4L))
  stopifnot("maxit must be a positive integer" = maxit > 0)
  stopifnot("stoptol must be a positive number" = stoptol > 0)

  A <- Matrix::Matrix(Ainput, sparse = TRUE)
  n <- length(b)
  p <- ncol(A)

  eigsopt <- list(issym = TRUE)
  tstartLip <- Sys.time()
  Lip <- RSpectra::eigs(
    tcrossprod(A),
    k = 1, which = "LA", opts = eigsopt, n = n
  )$values
  message(sprintf(
    "\n Lip = %3.2e, time = %3.2f",
    Lip, as.numeric(difftime(Sys.time(), tstartLip, units = "secs"))
  ))

  y <- rep(0, n)
  z <- rep(0, p)
  x <- z
  if (!is.null(y0) && !is.null(z0) && !is.null(x0)) {
    y <- y0
    z <- z0
    x <- x0
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
    maxfeas = maxfeas
  )
  message("\n****************************************")
  message(sprintf(" SSNAL       : %s", msg))
  message(sprintf(" iteration   : %d", iter))
  message(sprintf(" time(s)     : %3.2f", info_main$ttime))
  if (printyes) {
    message(sprintf(" prim_obj    : %4.8e", primobj))
    message(sprintf(" dual_obj    : %4.8e", dualobj))
    message(sprintf(" relgap      : %4.5e", info_main$relgap))
    message(sprintf(" primfeas    : %3.2e", primfeas))
    message(sprintf(" dualfeas    : %3.2e", dualfeas))
    message(sprintf(" eta         : %3.2e", info_main$eta))
    message(sprintf(" nnz         : %d", runhist_main$nnz))
  }

  return(list(
    obj = obj, x = x, y = y, z = z, info = info
  ))
}
