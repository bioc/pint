W.effects <- function(model, X, Y = NULL){


  # Check if whole data is given instead window for this model
  if (class(X) == "list"){
    # Find correct window for this model
    index <- which(rownames(X$data) == getGeneName(model))

    # Check if model has only 1 variable from X data
    if (nrow(getW(model)$X) == 1) {
      window <- sparse.window(X, Y, index, getWindowSize(model))
      Xm <- window$X
      Ym <- window$Y
    } else {
      #window <- fixed.window(X, Y, index, getWindowSize(model))
      #X <- window$X
      #Y <- window$Y
      Xm <- t(centerData(t(X$data[rownames(model@W$X), ]), rm.na = TRUE))
      Ym <- t(centerData(t(Y$data[rownames(model@W$Y), ]), rm.na = TRUE))
    }
  }  
  z <- z.expectation(model, Xm, Ym)
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
      if (nrow(getW(model)$X) == 1) {
        window <- sparse.window(X, Y, index, getWindowSize(model))
        Xm <- window$X
        Ym <- window$Y
      } else {
        #window <- fixed.window(X, Y, index, getWindowSize(model))
        #X <- window$X
        #Y <- window$Y
        Xm <- t(centerData(t(X$data[rownames(model@W$X), ]), rm.na = TRUE))
        Ym <- t(centerData(t(Y$data[rownames(model@W$Y), ]), rm.na = TRUE))
      }
    }
    # Check the sign of projvec so that it corresponds with the data
    if (cor(c(Xm[,z.order[1]], Ym[,z.order[1]]), projvec) < 0)
      projvec <- -projvec

    Wx <- getW(model)$X
    # Divide to X and Y components
    projvecx <- projvec[(1:nrow(Wx))]
    projvecy <- projvec[-(1:nrow(Wx))]
   
    return(list(total = projvec, X = projvecx, Y = projvecy))
    
  } else {   # for models with one data set
    if (cor(Xm[,z.order[1]], projvec)  < 0)
      projvec <- -projvec

    return(list(total = projvec))
  }

}
