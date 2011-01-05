optimize.simCCA.W2 <- function(w, phi, Dim, Dcov, nullmat, epsilon =
                               1e-3, par.change = 1e6, cost.old = 1e6,
                               mySeed = 123, dz = NULL, priors = NULL) {

  # option to give prior on W
  # no analytical solution in general case

  # modified from optimize.fullcov function
  # input otherwise similar, except that w is a single matrix
  # dimX x dimZ (note that dimX = dimY as always with simcca)

  # if dimensionality not specified, use the dimensionality of
  # input w
  dz <- ifelse(is.null(dz), ncol(w), dz)

  # Ensure that the dimensionality of given w matches with given dz
  w <- w[, 1:dz]
  W <- list()
  W$X <- W$Y <- w
  W$total <- rbind(w, w)

  set.seed( mySeed )
  cost.new <- cost7(abs(as.vector(W$X)), phi, Dcov, Dim, priors)
  initcost <- cost.new

  phi.inv <- list()

  while (par.change > epsilon || par.change < 0) {

    cost.old <- cost.new
    
    phi.inv$X <- solve(phi$X)
    phi.inv$Y <- solve(phi$Y)
    phi.inv$total <- rbind(cbind(phi.inv$X,nullmat), cbind(nullmat,phi.inv$Y))
    
    # assuming Wx = Wy!
    # see Bach-Jordan 2005, sec. 4.1 for details
    phi.inv.sum <- phi.inv$X + phi.inv$Y
    M <- solve(t(W$X)%*%phi.inv.sum%*%W$X + diag(dz))

    # store
    W.old <- W
    
    # Update W
    opt <- optim(as.vector(W$X), cost7, method = "L-BFGS-B", phi = phi, priors = priors, Dim = Dim, Dcov = Dcov, control = list(maxit = 1e6), lower = -10*max(Dcov$total), upper = 10*max(Dcov$total))    
    w <- get.W4(abs(opt$par), Dim)$X
    W$X <- W$Y <- w
    W$total <- rbind(w, w)
    W.new <- W
    
    rownames(W$X) <- rownames(Dcov$X)
    rownames(W$Y) <- rownames(Dcov$Y)
    rownames(W$total) <- c(rownames(W$X), rownames(W$Y))
    
    # Update phi
    phi <- phi.EM.simcca(Dcov, W.new, phi.inv, W.old, M)
    rownames(phi$X) <- colnames(phi$X) <- rownames(Dcov$X)
    rownames(phi$Y) <- colnames(phi$Y) <- rownames(Dcov$Y)
    rownames(phi$total) <- colnames(phi$total) <- c(rownames(phi$X), rownames(phi$Y))

    # Check marginal likelihood (-logP) for data to keep eye on convergence correctness
    # the smaller, the better are the parameters
    cost.new <- cost7(abs(as.vector(W$X)), phi, Dcov, Dim, priors)

    par.change <- (cost.old - cost.new)

  }

  list(W = W, phi = phi, score = cost.new)

}
