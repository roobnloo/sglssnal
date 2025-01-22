sglssnal_main <- function(A, b, lambda, P, parmain, y, z, x) {
  tstart <- parmain$tstart
  Lip <- parmain$Lip
  maxit <- parmain$maxit
  printyes <- parmain$printyes
  stoptol <- parmain$stoptol
  stoptol_gap <- parmain$stoptol_gap
  stopopt <- parmain$stopopt
  p <- parmain$p
  n <- parmain$n
  c1 <- lambda[1]
  c2 <- lambda[2]
  stop <- 0

  sigmaLip <- 1 / Lip

  Aty <- crossprod(A, y)
  Ax <- A %*% x
  Px <- P$matrix %*% x

  obj <- numeric(2)
  obj[1] <- 0.5 * sum((Ax - b)^2) + c2 * P$Lasso_fz(Px) + c1 * sum(abs(x))
  obj[2] <- 0
  normb <- 1 + norm(b, type = "2")
  sig <- max(1 / sqrt(Lip), min(c(1, sigmaLip, c1, 1 / c1, 1 / c2)))
  Rp <- Ax - b - y
  Rd <- Aty + z
  primfeas <- norm(Rp, type = "2") / normb
  dualfeas <- norm(Rd, type = "2") / (1 + norm(z, type = "2"))
  maxfeas <- max(primfeas, dualfeas)
  relgap <- abs(obj[1] - obj[2]) / (1 + abs(obj[1]) + abs(obj[2]))

  if (printyes) {
    message(sprintf(
      "\n n=%d, m=%d, tol=%1.1e, parameters:c1=%4.3f, c2=%4.3f\n",
      p, n, stoptol, c1, c2
    ))
    message(" ---------------------------------------------------")
    message(" iter|  [pinfeas  dinfeas]    relgap |    pobj          dobj       | time | sigma |")
    message("****************************************************************************************************")
    message(sprintf(
      " #%3.1d| [%3.2e %3.2e]  %- 3.2e | %- 5.4e %- 5.4e | %5.1f | %3.2e |",
      0, primfeas, dualfeas, relgap, obj[1], obj[2], Sys.time() - tstart, sig
    ))
  }

  parNCG <- list(tolconst = 0.5, p = p)
  maxitersub <- 10
  prim_win <- 0
  dual_win <- 0

  ssncgop <- list(tol = stoptol, precond = 0, printsub = TRUE)

  sigmamax <- 1e6
  sigmamin <- 1e-4

  runhist <- list()

  PP <- P
  PP$G <- P$G - 1L
  PP$ind[1, ] <- P$ind[1, ] - 1L
  PP$ind[2, ] <- P$ind[2, ] - 1L
  for (iter in 1:maxit) {
    parNCG$sigma <- sig

    if (dualfeas < 1e-5) {
      maxitersub <- max(maxitersub, 30)
    } else if (dualfeas < 1e-3) {
      maxitersub <- max(maxitersub, 30)
    } else if (dualfeas < 1e-1) {
      maxitersub <- max(maxitersub, 20)
    }
    ssncgop$maxitersub <- maxitersub

    result <- sglssn_conjgrad_interface(
      y, as.numeric(Aty), x, as.numeric(Ax), A, as.numeric(b), lambda[1], lambda[2],
      PP$matrix, PP$G, PP$ind, PP$num_group, parNCG,
      ssncgop$printsub, ssncgop$maxitersub, ssncgop$tol
    )
    y <- result$y
    z <- result$z
    Aty <- result$Aty
    x <- result$x
    Ax <- result$Ax
    parNCG <- result$par
    runhist_NCG <- result$runhist
    info_NCG <- result$info

    if (info_NCG$breakyes < 0) {
      parNCG$tolconst <- max(parNCG$tolconst / 1.06, 1e-3)
    }

    Rd <- Aty + z
    Px <- P$matrix %*% x
    normRd <- norm(Rd, type = "2")
    dualfeas <- normRd / (1 + norm(z, type = "2"))
    Axb <- Ax - b
    Rp <- Axb - y
    normRp <- norm(Rp, type = "2")
    primfeas <- normRp / normb

    lasso <- c2 * P$Lasso_fz(Px) + c1 * sum(abs(x))
    dualobj <- -sum(y^2) / 2 - crossprod(b, y)
    primobj <- sum(Axb^2) / 2 + lasso

    gap <- primobj - dualobj
    relgap <- abs(gap) / (1 + abs(primobj) + abs(dualobj))

    eta_1 <- NULL
    if (stopopt == 1) {
      stop <- max(dualfeas, relgap) < stoptol
    } else if (stopopt == 2) {
      eta_tmp <- max(dualfeas, primfeas)
      if (eta_tmp < stoptol) {
        eta_1 <- norm(x - proximal_combo(x + z, lambda, P), type = "2") /
          (1 + norm(x, type = "2"))
        stop <- eta_1 < stoptol
      }
    } else if (stopopt == 3) {
      stop <- (dualfeas < stoptol) && (relgap < stoptol_gap)
    } else if (stopopt == 4) {
      stop <- (dualfeas < stoptol) && (gap < stoptol * sum(b^2))
    }

    ttime <- Sys.time() - tstart
    runhist$primfeas[iter] <- primfeas
    runhist$dualfeas[iter] <- dualfeas
    runhist$sigma[iter] <- sig
    runhist$primobj[iter] <- primobj
    runhist$dualobj[iter] <- dualobj
    runhist$gap[iter] <- gap
    runhist$relgap[iter] <- relgap
    runhist$ttime[iter] <- ttime
    runhist$itersub[iter] <- info_NCG$itersub

    if (printyes) {
      message(
        sprintf(
          "\n%5.0d| [%3.2e %3.2e]  %- 3.2e | %- 5.4e %- 5.4e | %5.1f | %3.2e |",
          iter, primfeas, dualfeas, relgap, primobj, dualobj, ttime, sig
        )
      )
    }

    if (stop || (iter == maxit)) {
      termination <- "converged"
      if (iter == maxit) termination <- "maxiter reached"
      runhist$termination <- termination
      runhist$iter <- iter
      runhist$nnz <- cardcal(x)
      obj[1] <- primobj
      obj[2] <- dualobj
      break
    }

    ratio <- primfeas / (dualfeas + .Machine$double.eps)
    runhist$ratio_seq[iter] <- ratio

    if (ratio < 1) {
      prim_win <- prim_win + 1
    } else {
      dual_win <- dual_win + 1
    }

    sigma_update_iter <- sigma_fun(iter)
    if (primfeas > 100 * stoptol) {
      if (runhist_NCG$av_findstep > 5) {
        sigmascale <- sqrt(3)
      } else {
        sigmascale <- 3
      }
    } else {
      if (runhist_NCG$av_findstep > 5) {
        sigmascale <- sqrt(5)
      } else {
        sigmascale <- 5
      }
    }

    if ((iter %% sigma_update_iter == 0) && info_NCG$breakyes < 0) {
      if (prim_win > max(1, 1.2 * dual_win)) {
        prim_win <- 0
        sig <- min(sigmamax, sig * sigmascale)
      } else if (dual_win > max(1, 1.2 * prim_win)) {
        dual_win <- 0
        sig <- max(sigmamin, sig / sigmascale)
      }
      if (info_NCG$breakyes >= 0 && iter >= 10) {
        sig <- max(sigmamin, 2 * sig / sigmascale)
      }
    }
  }

  info <- list()

  if (stopopt == 2) {
    kktres <- max(eta_tmp, eta_1)
    runhist$kktres <- kktres
    info$kktres <- kktres
  }

  info$maxfeas <- maxfeas
  info$iter <- iter
  info$ttime <- ttime
  info$termination <- termination
  info$relgap <- relgap
  info$msg <- termination

  return(list(obj = obj, y = y, z = z, x = x, info = info, runhist = runhist))
}

sigma_fun <- function(iter) {
  if (iter < 10) {
    return(2)
  } else if (iter < 20) {
    return(3)
  } else if (iter < 200) {
    return(3)
  } else if (iter < 500) {
    return(10)
  }
}
