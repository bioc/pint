z.expectation <- function (model, X, Y = NULL) {

  W <- getW(model)
  phi <- getPhi(model)
  
  if (!is.null(Y)){
    # Check if whole data is given instead window for this model
    if(class(X) == "list"){
      # Find correct window for this model
      index <- which(rownames(X$data) == getGeneName(model))
    
      # Check if model has only 1 variable from X data
      if (nrow(W$X) == 1)
        window <- sparse.window(X, Y, index, getWindowSize(model))
      else
        window <- fixed.window(X, Y, index, getWindowSize(model))
      X <- window$X
      Y <- window$Y
    }

	
    S <- solve(t(W$X)%*%solve(phi$X)%*%W$X + t(W$Y)%*%solve(phi$Y)%*%W$Y + diag(ncol(W$X)))

    return(S%*%(t(W$X)%*%solve(phi$X)%*%X + t(W$Y)%*%solve(phi$Y)%*%Y))
  }
  else {
    W <- W$total
    phi <- phi$total
    S <- solve(t(W) %*% solve(phi) %*% W + diag(ncol(W)))
    return(S %*% (t(W) %*% solve(phi) %*% X))
  }
}