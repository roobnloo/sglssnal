proximal_combo <- function(v, c, P) {
  c1 <- c[1]
  c2 <- c[2]

  utmp <- proximal_l1(v, c1)
  u <- P$times(utmp)
  result <- P$ProjL2(u, c2)
  u <- result$Pz
  u <- utmp - P$trans(u)

  return(u)
}
