
set.M.full2 <- function (W, phi.inv, dz) {
  # for full marginal covariance
  solve(t(W$X)%*%phi.inv$X%*%W$X + t(W$Y)%*%phi.inv$Y%*%W$Y + diag(dz))
}

