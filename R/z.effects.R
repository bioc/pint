z.effects <- function(model, X, Y = NULL){

  W <- getW(model)

  # for models from 2 data sets
  if (!is.null(Y)){
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
    
    # Check that data window is smaller than half the sample size
    if (getWindowSize(model) > ncol(X))
      stop("Contribution of samples cannot be calculated when the window size is more than half the number of samples")
    W <- W$total

    z <- z.expectation(model,X,Y)

    # Calculate first component of PCA for W*z
    pca <- princomp(t(W%*%z))
    projvec <- pca$loadings[,1]

    # Project data to this component
    data <- rbind(X,Y)
    proj <- t(data)%*%projvec
   
    return(proj)
  }
  # for models with one data set
  else {
    W <- W$total
    z <- z.expectation(model,X)
    
    # Calculate first component of PCA for W*z
    pca <- princomp(t(W%*%z))
    projvec <- pca$loadings[,1]

    # Project data to this component
    data <- X
    proj <- t(data)%*%projvec
   
    return(proj)
  
  }
}




