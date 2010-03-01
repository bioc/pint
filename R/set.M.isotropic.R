
set.M.isotropic <- function (wtw, sigma.sq, dx) {
  solve(wtw + diag(sigma.sq, dx, dx))
}

