


ppca = function (X,Y=NULL,zDimension=1) {

  # Replaces function solve.CCA

  # if zDimension = NULL then full-rank solution is computed
  
  # Probabilistic PCA
  # (See Tipping and Bishop 1999)

  # ML estimates W, sigma for probability model
  # X ~ N(Wz, sigma*I)
  # i.e. latent variable model with isotropic noise

  # If only X is given in the argument, compute
  # pPCA for X

  # If both X and Y are given in the argument, compute
  # pPCA for concatenated [X; Y]
  # Assuming isotropic and identical marginal noise, the
  # principal subspace will capture the dependencies between X and Y.
  # This corresponds to the model (sigmax = sigmay = sigma)
  # X ~ N(Wx*z, sigma*I); Y ~ N(Wy*z, sigma*I)
  # This provides a simple comparison method for more
  # detailed dependency models.

  if (is.null(Y)) {
    res <- ppca.calculate(X, zDimension)
  } else {
    # If second argument (Y) given, compute
    # pPCA with two concatenated data sets
    res <- ppca.calculate(rbind(X,Y), zDimension)

    # Make phi diagonal matrix
    phitotal <- diag(res$phi,(nrow(X)+nrow(Y)))

    # Variable names to W and phi
    rownames(res$W) <- c(rownames(X),rownames(Y))
    rownames(phitotal) <- c(rownames(X),rownames(Y))
    colnames(phitotal) <- rownames(phitotal)

    # Divide W and phi to X and Y parts
    phi <- list(X = phitotal[(1:nrow(X)),(1:nrow(X))], Y = phitotal[-(1:nrow(X)),-(1:nrow(X))],
    	        total = phitotal)
    W <- list(X = as.matrix(res$W[(1:nrow(X)),]), Y = as.matrix(res$W[-(1:nrow(X)),]), total = res$W)

  }
  # Note that if X, Y given then phi$X = phi$Y in the pCCA model
  # Here W corresponds to W$total of other functions when X, Y both given
  list(W=W, phi=phi)
}
