proximal_l1 <- function(w, kappa) {
  if (kappa < 0) {
    stop("kappa needs to be non-negative")
  }
  v <- pmax(0, w - kappa) - pmax(0, -w - kappa)
  return(v)
}
