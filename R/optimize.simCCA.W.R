optimize.simCCA.W <- function (w, phi, Dim, Dcov, nullmat, epsilon = 1e-3,
                               par.change = 1e6, cost.old = 1e6, mySeed = 123, dz = NULL) { 

  #optimize.simCCA.W(
  #w <- W.init$X; phi <- phi.init; epsilon <- 1e-3; par.change = 1e6; cost.old = 1e6; mySeed=124; dz = zDimension

  # assuming Wx = Wy!
 
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

  set.seed(mySeed)
  cost.new <- cost5(W$X, phi, Dcov)
  initcost <- cost.new

  phi.inv <- list()

  while (par.change>epsilon || par.change < 0) {

    cost.old <- cost.new
    
    phi.inv$X <- solve(phi$X)
    phi.inv$Y <- solve(phi$Y)
    phi.inv$total <- rbind(cbind(phi.inv$X,nullmat), cbind(nullmat,phi.inv$Y))
    
        # assuming Wx = Wy!
        # see Bach-Jordan 2005, sec. 4.1 for details
    phi.inv.sum <- phi.inv$X + phi.inv$Y
    M <- solve(t(W$X)%*%phi.inv.sum%*%W$X + diag(dz))
        #beta <- M%*%t(W$X)%*%phi.inv.sum # ct. set.beta.fullcov(M, W$total, phi.inv$total)

        # store
    W.old <- W

        ## EM update
        ##W$total = W.simcca.EM(Dcov, M, beta)
        #W$total = Dcov$total%*%t(beta)%*%solve(M + beta%*%Dcov$total%*%t(beta))
        #W$X = W$total[1:Dim$X,]
        #W$Y = W$total[-c(1:Dim$X),]
        ## robust solution: take average of the two estimates
        #W$X = W$Y = (W$X+W$Y)/2
        #W$total = rbind(W$X,W$Y)

        # alternative update
    phihat <- phi$X + phi$Y
    phihat.inv <- solve(phihat)
    what <- 2*W$X
    M.w <- solve(t(what)%*%phihat.inv%*%what + diag(dz))
    beta.w <- M.w%*%t(what)%*%phihat.inv
        # Update what
    what <- Dcov$sum%*%t(beta.w)%*%solve(M.w + beta.w%*%Dcov$sum%*%t(beta.w))
    w <- what/2
    W$X <- W$Y <- w
    W$total <- rbind(w,w)
    W.new <- W

        # Update phi
    phi <- phi.EM.simcca(Dcov, W.new, phi.inv, W.old, M)

        # Check marginal likelihood (-logP) for data to keep eye on convergence correctness
        # the smaller, the better are the parameters
    cost.new <- cost5(W$X, phi, Dcov)

    par.change <- (cost.old - cost.new)#/cost.old

  }

  list(W = W, phi = phi, score = cost.new)

}









