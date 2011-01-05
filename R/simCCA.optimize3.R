simCCA.optimize3 <-
function (X, Y, zDimension = 1, mySeed = 123, epsilon = 1e-6, priors = NULL) {
  
  # Different from simCCA.optimize.R in that T is not optimized here
  # (not included in the model) and there is option to set prior on W
  # (W.prior)

  # zDimension: assumed dimensionality for the shared latent variable z

  # w.similarity parameter tuning Wx ~ Wy in the model where 
  # sigma.w <- priors$sigma.w
  # Wx - Wy ~ Nm(O, sigma.w*I, sigma.w*I)
  # sigma.w <- 0 : Wx = Wy
  # sigma.w <- Inf : Wx, Wy free
  
  set.seed(mySeed)
  
  #################################################
  # Initialize
  #################################################

  # samples are always matched
  Nsamples <- ncol(X)

  if ( length(priors) == 0 ) { priors <- list() }
  if ( is.null(priors$sigma.w) ) { priors$sigma.w <- Inf } # tunes similarity constraint Wx ~ Wy

  Dim <- list()
  Dim$X <- nrow(X)
  Dim$Y <- nrow(Y)
  Dim$Z <- zDimension

  Dcov <- list()
  Dcov$X <- cov(t(X))
  Dcov$Y <- cov(t(Y))
  Dcov$total <- cov(t(rbind(X, Y)))

  # scalar diagonal noise on the marginals (shared by all features) initially
  phi.init <- list(X = diag(var(as.vector(X)), Dim$X), Y = diag(var(as.vector(Y)), Dim$Y)) 
  
  # Initialize 
  W.init   <- list()
  W.init$X <- as.matrix(eigen(Dcov$X)$vectors[, 1:Dim$Z])
  W.init$Y <- as.matrix(eigen(Dcov$Y)$vectors[, 1:Dim$Z])
  W.init$total <- rbind(W.init$X, W.init$Y) # can be possibly removed in some special cases
  
  ##################################################
  # optimize until convergence
  ##################################################

  res <- optimize.W3(W.init, phi.init, Dim, Dcov, priors, epsilon, par.change = 1e6, cost.old = 1e6, mySeed = mySeed + 1)

  W <- list(X = res$W$X, Y = res$W$Y, total = rbind(res$W$X, res$W$Y))
  phi <- res$phi


  
  rownames(W$X) <- rownames(X)
  rownames(W$Y) <- rownames(Y)
  rownames(W$total) <- c(rownames(X), rownames(Y))
      
  #phiX <- diag(res$phi$X, nrow(X))
  rownames(phi$X) <- colnames(phi$X) <- rownames(X)

  #phiY <- diag(res$phi$Y, nrow(Y))
  rownames(phi$Y) <- colnames(phi$Y) <- rownames(Y)

  #phitotal <- diag(c(diag(phiX), diag(phiY)), nrow(X) + nrow(Y))
  rownames(phi$total) <- colnames(phi$total) <- c(rownames(X), rownames(Y))

  #phi <- list(X = phiX, Y = phiY, total = phitotal)

  return( list(W = W, phi = phi) )

}

