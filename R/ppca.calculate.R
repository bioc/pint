ppca.calculate <- function (X, zDimension) {

 # FIXME: ensure that X is zero-mean

  # Probabilistic PCA
  # (See Tipping and Bishop 1999 / 3.2)

  # ML estimates W, sigma for probability model
  # X ~ N(Wz, sigma*I)
  # i.e. latent variable model with isotropic noise

  # Use full-rank if dimensionality is not specified
  zDimension <- ifelse(is.null(zDimension), nrow(X), zDimension)

  # eigenvalues D and eigenvectors U
    duv <- svd(X)
    U <- duv$u
    D <- sort(duv$d, decreasing=TRUE)

    # ML estimate of the variance given dimensionality zDimension for the latent
    # variable z in model X ~ Wz + N(0,sigma) where sigma is isotropic
    d <- nrow(U)
    phi <- sum(D[-seq(zDimension)])/(d-zDimension) 

    # ML estimate for W given variance
    # Here set R <- I (R is an arbitrary orthogonal rotation matrix)
    
    W <- as.matrix(U[,1:zDimension])%*%sqrt(diag(D)[seq(zDimension),seq(zDimension)]-phi*diag(1,zDimension,zDimension))


  if (zDimension == nrow(X)) {

      # If W is full-rank then the isotropic error term disappears assuming data X is gaussian
      # then X%*%t(X) i.e. cov(t(X)) approximates W%*%t(W) (since X ~ Wz and z ~ N(0,I))

    # Note rotational ambiguity for W, Z 
    cat("Full-rank PCA calculated, isotropic error term is zero. Consider checking the principal components.\n")
    W <- matrix.sqrt(cov(t(X)))
    phi <- 0
  }
  
    list(W = W, phi = phi)
}
