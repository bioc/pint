get.W2 <-
function (vec, Dim) {

  # Convert parameter vector into matrices Wx and Wy

  # Note that we always assume that W is positive Therefore remove the
  #sign here to speed up optimization
  vec <- abs(vec)
  
  # NOTE: assumes dx = dy  
  W <- list()
  nx <- (Dim$X*Dim$Z)
  W$X <- array(vec[1:nx], dim = c(Dim$X, Dim$Z))
  W$Y <- array(vec[(nx + 1):length(vec)], dim = c(Dim$Y, Dim$Z))
  W$total <- rbind(W$X, W$Y)
  
  W

}

