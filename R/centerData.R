centerData <-
function (X, rm.na = FALSE, meanvalue = NULL) {

  #remove col means from matrix X

  if (rm.na) {
    X2 <- array(NA, dim = dim(X), dimnames = dimnames(X))
    for (i in 1:ncol(X)) {
      x <- X[, i]
      nainds <- is.na(x)
      xmean <- mean(x[!nainds])
      X2[!nainds, i] <- x[!nainds] - xmean 
    }
  }
  if (!rm.na) {
    X2 <- apply(X, 2, function(x){ x - mean(x) })
  }
  if (length(meanvalue) > 0) {
    # Shift the data so that mean gets a specified value
    X2 <- X2 + meanvalue
  }
  
  X2
}

