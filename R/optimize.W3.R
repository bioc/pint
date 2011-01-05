optimize.W3 <-
function (W, phi, Dim, Dcov, priors, epsilon = 1e-6, par.change = 1e6, cost.old = 1e6, mySeed = 123) {
  
  set.seed(mySeed)
  cost.new <- cost.W2(c(as.vector(W$X), as.vector(W$Y)), phi, priors, Dim, Dcov)
  costs    <- c(cost.new)
  cnt      <- 1
  phi.inv  <- list()
  nullmat  <- matrix(0, nrow = Dim$X, ncol = Dim$Y)

  # if phi$Y is scalar (as in segmented/mir case) we can speed up here. Do later.
  phi.inv$X <- solve(phi$X)
  phi.inv$Y <- solve(phi$Y)  
  phi.inv$total <- rbind(cbind(phi.inv$X, nullmat), cbind(t(nullmat), phi.inv$Y))

  while (par.change > epsilon) {

    cost.old <- cost.new

    # Update W: initialize with previous W	
    
    # This optimizes Wx and Wy separately
    opt <- optim(c(as.vector(W$X), as.vector(W$Y)), cost.W2, method = "L-BFGS-B", phi = phi, priors = priors, Dim = Dim, Dcov = Dcov, control = list(maxit = 1e6), lower = -10*max(Dcov$total), upper = 10*max(Dcov$total))

    # Convert optimized W parameter vector to actual matrices
    # Note that here we always assume that W is positive

    W.old <- W
    W.new <- W <- get.W2(opt$par, Dim)

    # Update phi (full covariances)
    phi.inv$X <- solve(phi$X)
    phi.inv$Y <- solve(phi$Y)    
    phi.inv$total <- rbind(cbind(phi.inv$X, nullmat),
                           cbind(t(nullmat), phi.inv$Y))    

    # also check from optimize.fullcov.R
    #M1 <- set.M.full(W.old$total, phi.inv$total, dz = Dim$Z)  # for non-matched case
    M <- set.M.full2(W.old, phi.inv, dz = Dim$Z)
    # was Checked: both M and M1 give the same result
    
    phi <- update.phi.EM3(Dcov, W.new, phi.inv, W.old, M, nullmat)

    cnt <- cnt + 1

    # Check and print marginal likelihood (-logP) for the data
    # the smaller, the better are the parameters

    cost.new <- cost.W2(c(as.vector(W$X), as.vector(W$Y)), phi, priors, Dim, Dcov)
    
    par.change <- (cost.old - cost.new)
    costs[[cnt]] <- cost.new

  }

  list(W = W, phi = phi, costs = costs, score = costs[[length(costs)]])
}

