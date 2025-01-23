#' Run Sparse-Group Lasso via Semismooth Newton Augmented Lagrangian
#' @description Fits a sparse-group lasso model using second-order information.
#' @param Ainput \eqn{n \times p} design matrix
#' @param b \eqn{n} response vector
#' @param lambda vector of length 2, `lambda[1]` is the
#'   \eqn{\ell_1} penalty and `lambda[2]` is the \eqn{\ell_2} penalty
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
    Ainput, b, lambda, G, ind, options,
    y0 = NULL, z0 = NULL, x0 = NULL) {
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

  n <- length(b)

  if (length(lambda) != 2) stop("lambda must be a vector of length 2")
  if (lambda[1] < 0 || lambda[2] < 0) stop("lambda must be non-negative")
  if (n != nrow(Ainput)) stop("nrow(A) must be equal to length(b)")

  normb <- 1 + norm(b, "2")
  A <- Ainput
  p <- ncol(A)

  tstart <- Sys.time()
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

  # Group map
  P <- group_structure(p, G, ind)
  PP <- P
  PP$G <- P$G - 1L
  PP$ind[1, ] <- P$ind[1, ] - 1L
  PP$ind[2, ] <- P$ind[2, ] - 1L

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

  result <- sglssnal_main_interface(
    A, b, lambda[1], lambda[2], PP$matrix, PP$G, PP$ind, PP$num_group,
    parmain, y, z, x
  )
  y <- result$y
  z <- result$z
  x <- result$x
  info_main <- result$info

  iter <- info_main$iter
  ttime <- as.numeric(difftime(Sys.time(), tstart, units = "secs"))
  msg <- info_main$msg
  if (iter == maxit) msg <- "maximum iteration reached"
  Aty <- crossprod(A, y)
  Ax <- A %*% x
  Rd <- Aty + z
  Px <- P$matrix %*% x
  normRd <- norm(Rd, "2")
  dualfeas <- normRd / (1 + norm(z, "2"))
  Axb <- Ax - b
  Rp <- Axb - y
  normRp <- norm(Rp, "2")
  primfeas <- normRp / normb
  maxfeas <- max(primfeas, dualfeas)

  eta <- sqrt(sum((x - proximal_combo(x + z, lambda, P))^2)) / (1 + sqrt(sum(x^2)))

  lasso <- lambda[2] * P$Lasso_fz(Px) + lambda[1] * sum(abs(x))
  dualobj <- -sum(y^2) / 2 - sum(b * y)
  primobj <- sum(Axb^2) / 2 + lasso
  obj <- c(primobj, dualobj)
  gap <- primobj - dualobj
  relgap <- abs(gap) / (1 + abs(primobj) + abs(dualobj))

  runhist <- list(
    totaltime = ttime,
    primobj = primobj,
    dualobj = dualobj,
    maxfeas = maxfeas,
    eta = eta
  )
  info <- list(
    relgap = relgap,
    iter = iter,
    time = ttime,
    eta = eta,
    obj = obj,
    maxfeas = maxfeas
  )
  if (printyes) {
    message("\n****************************************")
    message(sprintf(" SSNAL       : %s", msg))
    message(sprintf(" iteration   : %d", iter))
    message(sprintf(" time        : %3.2f", ttime))
    message(sprintf(" prim_obj    : %4.8e", primobj))
    message(sprintf(" dual_obj    : %4.8e", dualobj))
    message(sprintf(" relgap      : %4.5e", relgap))
    message(sprintf(" primfeas    : %3.2e", primfeas))
    message(sprintf(" dualfeas    : %3.2e", dualfeas))
    message(sprintf(" eta         : %3.2e", eta))
    message(sprintf(" nnz         : %d", cardcal(x)$k))
  }

  return(list(obj = obj, y = y, z = z, x = x, info = info, runhist = runhist))
}
