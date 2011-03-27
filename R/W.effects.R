W.effects <- function(model, X, Y = NULL){

  # Check if whole data is given instead window for this model
  if (class(X) == "list"){
    # Find correct window for this model
    index <- which(rownames(X$data) == getGeneName(model))

    # Check if model has only 1 variable from X data
    if (nrow(getW(model)$X) == 1)
      window <- sparse.window(X, Y, index, getWindowSize(model))
    else
      window <- fixed.window(X, Y, index, getWindowSize(model))
    X <- window$X
    Y <- window$Y
  }  
  z <- z.expectation(model, X, Y)
  W <- getW(model)$total

  # Sort samples according to their absolute effect
  z.eff <- z.effects(model, X, Y)
  z.order <- order(abs(z.eff),decreasing=TRUE)

  # Calculate first component of PCA for W*z
  pca <- princomp(t(W%*%z))
  projvec <- pca$loadings[,1]


  if (!is.null(Y)){
    if (class(X) == "list"){
      # Find correct window for this model
      index <- which(rownames(X$data) == getGeneName(model))

      # Check if model has only 1 variable from X data
      if (nrow(getW(model)$X) == 1)
        window <- sparse.window(X, Y, index, getWindowSize(model))
      else
        window <- fixed.window(X, Y, index, getWindowSize(model))
      X <- window$X
      Y <- window$Y
    }
    # Check the sign of projvec so that it corresponds with the data
    if (cor(c(X[,z.order[1]],Y[,z.order[1]]), projvec) < 0)
      projvec <- -projvec
  
  } else {
    if (cor(X[,z.order[1]], projvec)  < 0)
      projvec <- -projvec
  }

  # for models with 2 data sets
  if (!is.null(Y)){
    Wx <- getW(model)$X
    # Divide to X and Y components
    projvecx <- projvec[(1:nrow(Wx))]
    projvecy <- projvec[-(1:nrow(Wx))]
   
    return(list(total = projvec, X = projvecx, Y = projvecy))
  }
  # for models with one data set
  else {
    return(list(total = projvec))
  }
}
