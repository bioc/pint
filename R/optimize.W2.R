optimize.W2 <-
function (W, phi, Dim, Dcov, priors, epsilon = 1e-6, par.change = 1e6, cost.old = 1e6, mySeed = 123) {
  
  set.seed(mySeed)
  
  cost.new <- cost.W2(c(as.vector(W$X), as.vector(W$Y)), phi, priors, Dim, Dcov)
  
  costs    <- c(cost.new)
  cnt      <- 1
  phi.inv  <- list()
  nullmat  <- matrix(0, nrow = Dim$X, ncol = Dim$Y)
  
  while (par.change > epsilon) {
    
    cost.old <- cost.new

    phi.inv$X <- solve(phi$X)
    phi.inv$Y <- solve(phi$Y)

    phi.inv$total <- rbind(cbind(phi.inv$X, nullmat),
                           cbind(t(nullmat), phi.inv$Y))    

    M <- list()
    M$X <- set.M.full(W$X, phi$X, dz = Dim$Z)
    M$Y <- set.M.full(W$Y, phi$Y, dz = Dim$Z)
    #M <- set.M.full(W$total, phi.inv$total, dz = Dim$Z)  # for non-matched case later
    W.old <- W

    # Update W: initialize with previous W	

    if ((priors$sigma.w == Inf || priors$sigma.w == 0)) {

      # This optimizes Wx and Wy separately
      opt <- optim(c(as.vector(W$X), as.vector(W$Y)), cost.W2, method = "L-BFGS-B", phi = phi, priors = priors, Dim = Dim, Dcov = Dcov, control = list(maxit = 1e6), lower = -10*max(Dcov$total), upper = 10*max(Dcov$total))
      # Convert optimized W parameter vector to actual matrices
      # Note that here we always assume that W is positive

      W <- get.W2(opt$par, Dim)
      
    } else {
      stop("Case sigma.w in (0, Inf) not implemented yet.")
    }

    W.new <- W

    ###########################################################################

    # Update phi (full covariances)

    phi.inv$X <- solve(phi$X)
    phi.inv$Y <- solve(phi$Y)    
    phi.inv$total <- rbind(cbind(phi.inv$X, nullmat),
                           cbind(t(nullmat), phi.inv$Y))    

    # check from optimize.fullcov.R

    M <- solve(t(W$total)%*%phi.inv$total%*%W$total + diag(ncol(W$total)))

    phi <- update.phi.EM2(Dcov, W.new, phi.inv, W.old, M, nullmat)

    cnt <- cnt + 1

    # Check and print marginal likelihood (-logP) for the data
    # the smaller, the better are the parameters

    cost.new <- cost.W2(c(as.vector(W$X), as.vector(W$Y)), phi, priors, Dim, Dcov)

    par.change <- (cost.old - cost.new)
    costs[[cnt]] <- cost.new
  }

  list(W = W, phi = phi, costs = costs, score = costs[[length(costs)]])
}

