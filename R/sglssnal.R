#' Run Sparse-Group Lasso via Semismooth Newton Augmented Lagrangian
#' @description Fits a sparse-group lasso model using second-order information.
#' @param Ainput \eqn{n \times p} design matrix
#' @param b \eqn{n} response vector
#' @param lambda1 the \eqn{\ell_1} penalty
#' @param lambda2 the group-wise \eqn{\ell_2} penalty
#' @param G vector of group indices
#' @param ind \eqn{3 \times g} matrix describing the groups in `G`;
#'   `G[ind[1, i]:ind[2, i]]` are the indices of the `i`-th group and
#'   `ind[3, i]` is the weight of the group.
#' @param options control options
#' @param y0 optional initialization vector
#' @param z0 optional initialization vector
#' @param x0 optional initialization vector
#' @return list
#' @importFrom Rcpp sourceCpp
#' @useDynLib sglssnal
#' @export
sglssnal <- function(
    Ainput, b, lambda1, lambda2, G, ind, options,
    y0 = NULL, z0 = NULL, x0 = NULL) {
  stopifnot("lambda1 and lambda2 must be nonnegative" = lambda1 >= 0 && lambda2 >= 0)
  stopifnot("nrow(A) must be equal to length(b)" = nrow(Ainput) == length(b))

  stoptol <- 1e-6
  stoptol_gap <- stoptol
  printyes <- TRUE
  printsub <- TRUE
  maxit <- 5000
  stopopt <- 2 # 1:relgap+feas 2:kkt
  if (!is.null(options$stoptol)) stoptol <- options$stoptol
  if (!is.null(options$stopopt)) stopopt <- options$stopopt
  if (!is.null(options$printyes)) printyes <- options$printyes
  if (!is.null(options$printsub)) printsub <- options$printsub
  if (!printyes) printsub <- FALSE
  if (!is.null(options$maxit)) maxit <- options$maxit
  if (!is.null(options$stoptol_gap)) stoptol_gap <- options$stoptol_gap

  A <- Ainput
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
    stoptol_gap = stoptol_gap,
    Lip = Lip,
    maxit = maxit,
    stopopt = stopopt,
    printyes = printyes,
    printsub = printsub,
    p = p,
    n = n
  )
  if (!is.null(options$sigma)) parmain$sigma <- options$sigma

  gs <- group_structure(p, G, ind)
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
  obj <- c(primobj, dualobj)

  runhist <- list(
    totaltime = info_main$ttime,
    primobj = primobj,
    dualobj = dualobj,
    maxfeas = max(runhist_main$primfeas[iter], runhist_main$dualfeas[iter]),
    eta = info_main$eta
  )
  info <- list(
    relgap = info_main$relgap,
    iter = iter,
    time = info_main$ttime,
    eta = info_main$eta,
    obj = obj,
    maxfeas = maxfeas
  )
  if (printyes) {
    message("\n****************************************")
    message(sprintf(" SSNAL       : %s", msg))
    message(sprintf(" iteration   : %d", iter))
    message(sprintf(" time        : %3.2f", info_main$ttime))
    message(sprintf(" prim_obj    : %4.8e", primobj))
    message(sprintf(" dual_obj    : %4.8e", dualobj))
    message(sprintf(" relgap      : %4.5e", info_main$relgap))
    message(sprintf(" primfeas    : %3.2e", primfeas))
    message(sprintf(" dualfeas    : %3.2e", dualfeas))
    message(sprintf(" eta         : %3.2e", info_main$eta))
    message(sprintf(" nnz         : %d", cardcal(x)$k))
  }

  return(list(
    obj = obj, y = y, z = z, x = x, info = info, runhist = runhist
  ))
}
