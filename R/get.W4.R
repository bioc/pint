get.W4 <-
function (vec, Dim) {

  # Convert parameter vector into matrices Wx and Wy

  # NOTE: assumes dx = dy  
  W <- list()
  W$X <- W$Y <- array(abs(vec), dim = c(Dim$X, Dim$Z))
  # Note that we always assume that W is positive
  # Therefore remove the sign here to speed up optimization
  # FIXME: later, allow individual scales for Wx, Wy?
  
  W

}

