W.effects <- function(model,X,Y=NULL){

  z <- z.expectation(model,X,Y)
  W <- getW(model)$total

  # Calculate first component of PCA for W*z
  pca <- princomp(t(W%*%z))
  projvec <- pca$loadings[,1]

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

