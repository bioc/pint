get.W <-
function (vec, Dim) {

  # Convert parameter vector into matrices Wx and Wy
  W <- list()
  W$X <- array(vec[1:(Dim$X*Dim$Z)], dim = c(Dim$X,Dim$Z))
  T <- array(vec[-seq(Dim$X*Dim$Z)], dim = c(Dim$Y,Dim$X))
  W$Y <- T%*%W$X

  list(W = W, T = T)
}

