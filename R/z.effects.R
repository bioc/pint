z.effects <- function(model, X, Y = NULL){

  W <- getW(model)

  # for models from 2 data sets
  if (!is.null(Y)){
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
        Xm <- t(centerData(t(X$data[rownames(W$X), ]), rm.na = TRUE))
        Ym <- t(centerData(t(Y$data[rownames(W$Y), ]), rm.na = TRUE))
      }
    }
    
    # Check that data window is smaller than half the sample size
    # FIXME: not required in all models; loosen this where possible
    if (getWindowSize(model) > ncol(X$data)/2)
      stop("The number of samples must be at least two times higher than number of features (probes)")

    W <- W$total
    z <- z.expectation(model, Xm, Ym)
    
    # Calculate first component of PCA for W*z
    pca <- princomp(t(W%*%z))
    projvec <- pca$loadings[,1]

    # Project data to this component
    data <- rbind(Xm, Ym)
    proj <- t(data)%*%projvec
    
    # Make sure the highest value is always positive
    # FIXME: later adjust sign based on observed data to make directy interpretable
    # this is ok fix for now
    if (abs(min(proj)) > max(proj))
      proj <- -proj   

    return(proj)
  } else {   # for models with one data set
    W <- W$total
    z <- z.expectation(model, X$data[rownames(W$X), ])

    # FIXME: for clarity, make own function for this since same operation is
    # used also above
    # Calculate first component of PCA for W*z
    pca <- princomp(t(W%*%z))
    projvec <- pca$loadings[,1]
      
    # Project data to this component
    data <- t(centerData(t(X$data[rownames(W$X), ]), rm.na = TRUE))
    proj <- t(data)%*%projvec
   
    # Make sure the highest value is allways positive
    if (abs(min(proj)) > max(proj))
      proj <- -proj

    return(proj)
  
  }
}




