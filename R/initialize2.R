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
  Dcov$yx <- cov(t(Y), t(X), use = "pairwise.complete.obs")
  Dcov$total <- cov(t(rbind(X, Y)), use = "pairwise.complete.obs")
  Dcov$sum <- cov(t(X+Y), use = "pairwise.complete.obs")
  
  # Initialize 
  W <- list()
  W$X <- W$Y <- eigen(cov(t(X+Y), use = "pairwise.complete.obs"))$vectors
  W$total <- rbind(W$X,W$Y)
  
  phi.init <- list()
  phi.init$X <- Dcov$X
  phi.init$Y <- Dcov$Y
  phi.init$total <- rbind(cbind(phi.init$X,nullmat), cbind(nullmat, phi.init$Y))
  
  list(phi = phi.init, W = W, Dcov = Dcov, Dim = Dim, nullmat = nullmat, Nsamples = Nsamples)

}

