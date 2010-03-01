
set.M.full <- function (W, phi.inv, dz) {
  # for full marginal covariance
  solve(t(W)%*%phi.inv%*%W + diag(dz))
}

