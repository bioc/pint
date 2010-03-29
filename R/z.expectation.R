z.expectation <- function (model, X, Y) {

  W <- getW(model)
  phi <- getPhi(model)
  index <- which(dimnames(X$data)[[1]] == getGeneName(model))
	
  if (nrow(W$X) == 1)
    window <- sparse.window(X, Y, index, getWindowSize(model))
  else
    window <- fixed.window(X, Y, index, getWindowSize(model))
  X <- window$X
  Y <- window$Y
	
  S <- solve(t(W$X)%*%solve(phi$X)%*%W$X + t(W$Y)%*%solve(phi$Y)%*%W$Y + diag(ncol(W$X)))

  return(S%*%(t(W$X)%*%solve(phi$X)%*%X + t(W$Y)%*%solve(phi$Y)%*%Y))
}