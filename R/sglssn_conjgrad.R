sglssn_conjgrad <- function(y0, Aty0, x0, Ax0, A, b, lambda, P, par, options) {
  printsub <- 1
  breakyes <- 0
  maxitersub <- 50
  tiny <- 1e-10
  tol <- 1e-6
  maxitpsqmr <- 500
  precond <- 0

  if (!is.null(options$printsub)) printsub <- options$printsub
  if (!is.null(options$maxitersub)) maxitersub <- options$maxitersub
  if (!is.null(options$tol)) tol <- options$tol

  sig <- par$sigma
  normb <- 1 + norm(b, "2")

  y <- y0
  Aty <- Aty0
  u <- x0 / sig - Aty
  Prox_u <- proximal_combo(u, lambda, P)
  sigProx_u <- sig * Prox_u
  z <- u - Prox_u
  psi_y <- -(sum(b * y) + 0.5 * sum(y^2) + 0.5 * sig * sum(Prox_u^2))

  runhist <- list(
    priminf = numeric(maxitersub),
    dualinf = numeric(maxitersub),
    psi_y = numeric(maxitersub)
  )

  for (itersub in 1:maxitersub) {
    x <- sigProx_u
    Ax <- A %*% x
    Grad <- as.numeric(-y - b + Ax)
    normGrad <- norm(Grad, "2") / normb
    priminf_sub <- normGrad
    Rd <- Aty + z
    normRd <- norm(Rd, "2")
    dualinf_sub <- normRd / (1 + norm(z, "2"))

    if (max(priminf_sub, dualinf_sub) < tol) {
      tolsubconst <- 0.9
    } else {
      tolsubconst <- 1e-2
    }

    tolsub <- max(min(1, par$tolconst * dualinf_sub), tolsubconst * tol)
    runhist$priminf[itersub] <- priminf_sub
    runhist$dualinf[itersub] <- dualinf_sub
    runhist$psi_y[itersub] <- psi_y

    if (printsub) {
      message(sprintf(
        "\n      %2.0d  %- 11.10e  %3.2e   %3.2e %3.2e",
        itersub, psi_y, priminf_sub, dualinf_sub, par$tolconst
      ), appendLF = FALSE)
    }

    if (normGrad < tolsub && itersub > 1) {
      msg <- "good termination in subproblem:"
      if (printsub) {
        message(sprintf("\n       %s  ", msg), appendLF = FALSE)
        message(sprintf(
          " dualinfes = %3.2e, normGrad = %3.2e, tolsub = %3.2e",
          dualinf_sub, normGrad, tolsub
        ))
      }

      breakyes <- -1
      break
    }

    par$epsilon <- min(1e-3, 0.1 * normGrad)
    par$precond <- precond
    if (precond == 1) {
      par$invdiagM <- 1 / (1 + sig)
    }

    if (dualinf_sub > 1e-3 || itersub <= 5) {
      maxitpsqmr <- max(maxitpsqmr, 200)
    } else if (dualinf_sub > 1e-4) {
      maxitpsqmr <- max(maxitpsqmr, 300)
    } else if (dualinf_sub > 1e-5) {
      maxitpsqmr <- max(maxitpsqmr, 400)
    } else if (dualinf_sub > 5e-6) {
      maxitpsqmr <- max(maxitpsqmr, 500)
    }

    if (itersub > 1) {
      prim_ratio <- priminf_sub / runhist$priminf[itersub - 1]
      dual_ratio <- dualinf_sub / runhist$dualinf[itersub - 1]
    } else {
      prim_ratio <- 0
      dual_ratio <- 0
    }

    rhs <- Grad
    tolpsqmr <- min(5e-3, 0.001 * norm(rhs, "2"))
    const2 <- 1
    if (itersub > 1 && (prim_ratio > 0.5 || priminf_sub > 0.1 * runhist$priminf[1])) {
      const2 <- 0.5 * const2
    }
    if (dual_ratio > 1.1) const2 <- 0.5 * const2
    tolpsqmr <- const2 * tolpsqmr
    par$tol <- tolpsqmr
    par$maxit <- maxitpsqmr
    nnz <- sum(abs(x) > tol)
    par$nnz <- nnz

    result <- conjgrad_linsolver(A, rhs, u, lambda, P, par)
    dy <- result$dy
    resnrm <- result$resnrm
    solve_ok <- result$solve_ok

    Atdy <- crossprod(A, dy)
    iterpsqmr <- length(resnrm) - 1
    if (printsub) {
      message(sprintf(
        " | %3.1e %3.1e %3.0d  %2.1f %2.0d",
        par$tol, resnrm[length(resnrm)], iterpsqmr, const2, nnz
      ), appendLF = FALSE)
    }

    par$iter <- itersub

    step_opt <- list(stepop = ifelse((itersub <= 3 && dualinf_sub > 1e-4) || (itersub < 3), 1, 2))
    steptol <- 1e-5

    step_result <- findstep(b, sig, psi_y, u, Prox_u, sigProx_u, z, y, Aty, dy, Atdy, lambda, P, steptol, step_opt)
    psi_y <- step_result$psi_y
    u <- step_result$u
    Prox_u <- step_result$Prox_u
    sigProx_u <- step_result$sigProx_u
    z <- step_result$z
    y <- step_result$y
    Aty <- step_result$Aty
    alp <- step_result$alp
    iterstep <- step_result$iter

    runhist$solve_ok[itersub] <- solve_ok
    runhist$psqmr[itersub] <- iterpsqmr
    runhist$findstep[itersub] <- iterstep
    runhist$av_findstep <- mean(runhist$findstep)

    if (alp < tiny) breakyes <- 11
    psiy_ratio <- 1
    if (itersub > 1) {
      psiy_ratio <- (psi_y - runhist$psi_y[itersub - 1]) / (abs(psi_y) + .Machine$double.eps)
    }

    if (printsub) {
      message(sprintf(" %3.2e %2.0f", alp, iterstep), appendLF = FALSE)
      if (psiy_ratio < 0) message("-", appendLF = FALSE)
    }

    if (itersub > 4) {
      idx <- max(1, itersub - 3):itersub
      tmp <- runhist$priminf[idx]
      ratio <- min(tmp) / max(tmp)
      if (all(runhist$solve_ok[idx] <= -1) && ratio > 0.9 && min(runhist$psqmr[idx]) == max(runhist$psqmr[idx]) && max(tmp) < 5 * tol) {
        message("#", appendLF = FALSE)
        breakyes <- 1
      }

      const3 <- 0.7
      priminf_1half <- min(runhist$priminf[1:ceiling(itersub * const3)])
      priminf_2half <- min(runhist$priminf[ceiling(itersub * const3) + 1:itersub])
      priminf_best <- min(runhist$priminf[1:(itersub - 1)])
      priminf_ratio <- runhist$priminf[itersub] / runhist$priminf[itersub - 1]
      dualinf_ratio <- runhist$dualinf[itersub] / runhist$dualinf[itersub - 1]
      stagnate_idx <- which(runhist$solve_ok[1:itersub] <= -1)
      stagnate_count <- length(stagnate_idx)
      idx2 <- max(1, itersub - 7):itersub

      if (itersub >= 10 && all(runhist$solve_ok[idx2] == -1) && priminf_best < 1e-2 && dualinf_sub < 1e-3) {
        tmp <- runhist$priminf[idx2]
        ratio <- min(tmp) / max(tmp)
        if (ratio > 0.5) {
          if (printsub) message("##", appendLF = FALSE)
          breakyes <- 2
        }
      }

      if (itersub >= 15 && priminf_1half < min(2e-3, priminf_2half) && dualinf_sub < 0.8 * runhist$dualinf[1] && dualinf_sub < 1e-3 && stagnate_count >= 3) {
        if (printsub) message("###", appendLF = FALSE)
        breakyes <- 3
      }

      if (itersub >= 15 && priminf_ratio < 0.1 && priminf_sub < 0.8 * priminf_1half && dualinf_sub < min(1e-3, 2 * priminf_sub) && (priminf_sub < 2e-3 || (dualinf_sub < 1e-5 && priminf_sub < 5e-3)) && stagnate_count >= 3) {
        if (printsub) message(" $$", appendLF = FALSE)
        breakyes <- 4
      }

      if (itersub >= 10 && dualinf_sub > 5 * min(runhist$dualinf) && priminf_sub > 2 * min(runhist$priminf)) {
        if (printsub) message("$$$", appendLF = FALSE)
        breakyes <- 5
      }

      if (itersub >= 20) {
        dualinf_ratioall <- runhist$dualinf[2:itersub] / runhist$dualinf[1:(itersub - 1)]
        idx <- which(dualinf_ratioall > 1)
        if (length(idx) >= 3) {
          dualinf_increment <- mean(dualinf_ratioall[idx])
          if (dualinf_increment > 1.25) {
            if (printsub) message("^^", appendLF = FALSE)
            breakyes <- 6
          }
        }
      }

      if (breakyes > 0) {
        x <- sig * Prox_u
        Ax <- A %*% x
        break
      }
    }
  }

  info <- list(
    maxCG = max(runhist$psqmr),
    avgCG = sum(runhist$psqmr) / itersub,
    breakyes = breakyes,
    itersub = itersub,
    tolconst = par$tolconst
  )

  return(list(y = y, z = z, Aty = Aty, Prox_u = Prox_u, x = x, Ax = Ax, par = par, runhist = runhist, info = info))
}


findstep <- function(b, sig, psi_y0, u0, Prox_u0, sigProx_u0, z0, y0, Aty0, dy, Atdy, c, P, tol, options) {
  op <- options$stepop
  printlevel <- 1
  maxit <- ceiling(log(1 / (tol + .Machine$double.eps)) / log(2))
  c1 <- 1e-4
  c2 <- 0.9

  tmp1 <- crossprod(dy, y0 + b)
  tmp2 <- sum(dy^2)

  cross <- crossprod(Atdy, Prox_u0)
  g0 <- as.numeric(sig * cross - tmp1)
  if (g0 <= 0) {
    alp <- 0
    iter <- 0
    if (printlevel) {
      message(sprintf("\n Need an ascent direction, %2.1e  ", g0), appendLF = FALSE)
    }
    u <- u0
    y <- y0
    z <- z0
    Aty <- Aty0
    Prox_u <- Prox_u0
    sigProx_u <- sigProx_u0
    psi_y <- psi_y0
    return(list(psi_y = psi_y, u = u, Prox_u = Prox_u, sigProx_u = sigProx_u, z = z, y = y, Aty = Aty, alp = alp, iter = iter))
  }

  alp <- 1
  alpconst <- 0.5
  for (iter in 1:maxit) {
    if (iter == 1) {
      alp <- 1
      LB <- 0
      UB <- 1
    } else {
      alp <- alpconst * (LB + UB)
    }
    y <- y0 + alp * dy
    u <- u0 - alp * Atdy
    Prox_u <- proximal_combo(u, c, P)
    z <- u - Prox_u
    galp <- tmp1 + alp * tmp2
    sigProx_u <- sig * Prox_u
    galp <- as.numeric(crossprod(Atdy, sigProx_u) - galp)
    psi_y <- as.numeric(-(crossprod(b, y) + 0.5 * sum(y^2) + 0.5 * sig * sum(Prox_u^2)))

    if (printlevel > 1) {
      message(sprintf("\n --------------------------------------- "))
      message(sprintf(" alp = %2.2f, psi_y = %11.10e, psi_y0 = %11.10e", alp, psi_y, psi_y0))
      message(sprintf(" --------------------------------------- "))
    }

    if (iter == 1) {
      gLB <- g0
      gUB <- galp
      # message(sprintf("\n gLB: %.3f gUB: %.3f\n", gLB, gUB))
      if (sign(gLB) * sign(gUB) > 0) {
        if (printlevel) message("|", appendLF = FALSE)
        Aty <- Aty0 + alp * Atdy
        return(list(psi_y = psi_y, u = u, Prox_u = Prox_u, sigProx_u = sigProx_u, z = z, y = y, Aty = Aty, alp = alp, iter = iter))
      }
    }

    if (abs(galp) < c2 * abs(g0) && (psi_y - psi_y0 - c1 * alp * g0 > 1e-8 / max(1, abs(psi_y0)))) {
      if (op == 1 || (op == 2 && abs(galp) < tol)) {
        if (printlevel) message(":", appendLF = FALSE)
        Aty <- Aty0 + alp * Atdy
        return(list(psi_y = psi_y, u = u, Prox_u = Prox_u, sigProx_u = sigProx_u, z = z, y = y, Aty = Aty, alp = alp, iter = iter))
      }
    }

    if (sign(galp) * sign(gUB) < 0) {
      LB <- alp
      gLB <- galp
    } else if (sign(galp) * sign(gLB) < 0) {
      UB <- alp
      gUB <- galp
    }
  }

  if (iter == maxit) {
    Aty <- Aty0 + alp * Atdy
  }
  if (printlevel) message("m", appendLF = FALSE)

  return(list(psi_y = psi_y, u = u, Prox_u = Prox_u, sigProx_u = sigProx_u, z = z, y = y, Aty = Aty, alp = alp, iter = iter))
}
