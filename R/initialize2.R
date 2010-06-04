initialize2 <-
function (X, Y) {

  #################################################
  #  Initialization
  #################################################

  Nsamples <- ncol(X)
  
  Dim <- list()
  Dim$X <- nrow(X)
  Dim$Y <- nrow(Y)
  
  nullmat <- array(0, dim = c(Dim$X, Dim$Y))
  
  Dcov <- list()
  Dcov$X <- cov(t(X), use = "pairwise.complete.obs")
  Dcov$Y <- cov(t(Y), use = "pairwise.complete.obs")
  Dcov$xy <- cov(t(X), t(Y), use = "pairwise.complete.obs")
  Dcov$yx <- t(Dcov$xy)
  Dcov$total <- rbind(cbind(Dcov$X, Dcov$xy), cbind(Dcov$yx, Dcov$Y))
  Dcov$sum   <- Dcov$X + Dcov$Y + Dcov$xy + Dcov$yx
  #Dcov$yx <- cov(t(Y), t(X), use = "pairwise.complete.obs")
  #Dcov$total <- cov(t(rbind(X, Y)), use = "pairwise.complete.obs")
  #Dcov$sum <- cov(t(X + Y), use = "pairwise.complete.obs")
  
  # It is possible that covariances calculated with pairwise complete
  # observations are not positive semi-definite.
  # Check this. If not pos.sem.def, then replace with the closest
  # pos. semidefinite matrix
  if (any(eigen(Dcov$total)$values < 0)) {
    message("Covariance approximation used to avoid numerical instabilities.")
    Dcov$X  <- as.matrix(nearPD(Dcov$X)$mat)
    Dcov$Y  <- as.matrix(nearPD(Dcov$Y)$mat)
    Dcov$total <- rbind(cbind(Dcov$X, Dcov$xy), cbind(Dcov$yx, Dcov$Y))
    Dcov$sum   <- Dcov$X + Dcov$Y + Dcov$xy + Dcov$yx
  }

  # Initialize 
  W <- list()
  W$X <- W$Y <- eigen(Dcov$sum)$vectors
  W$total <- rbind(W$X,W$Y)
  
  phi.init <- list()
  phi.init$X <- Dcov$X
  phi.init$Y <- Dcov$Y
  phi.init$total <- rbind(cbind(phi.init$X,nullmat), cbind(nullmat, phi.init$Y))
  
  list(phi = phi.init, W = W, Dcov = Dcov, Dim = Dim, nullmat = nullmat, Nsamples = Nsamples)

}

