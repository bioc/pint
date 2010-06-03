simCCA.optimize.fullcov.EM <-
function (X, Y, zDimension = NULL, mySeed=123, epsilon = 1e-6) {

  set.seed(mySeed)

  #################################################
  # Initialize
  #################################################

  inits <- initialize2(X,Y)
  phi.init <- inits$phi
  W.init <- inits$W
  Dcov <- inits$Dcov
  Dim <- inits$Dim
  nullmat <- inits$nullmat
  Nsamples <- inits$Nsamples
  
  ##################################################
  # optimize until convergence
  ##################################################


  # use this for full W (EM algorithm, unstable for n ~ p)
  res <- optimize.simCCA.W(W.init$X, phi.init, Dim, Dcov,
                           nullmat, epsilon, par.change = 1e6,
                           mySeed=mySeed + 1, dz = zDimension)

  res

}

