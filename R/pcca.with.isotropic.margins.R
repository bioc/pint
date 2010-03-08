############################################################################

# Example use
#data(chromosome17)
#X = t(centerData(t(geneExp$data[145:154,])))
#Y = t(centerData(t(geneCopyNum$data[145:154,])))
#res <- pcca.with.isotropic.margins(X, Y,1,1e-6)

############################################################################

pcca.with.isotropic.margins <- function (X, Y, dz = 1, epsilon = 1e-6, delta = 1e6) {

  # epsilon and delta are convergence parameters
  # dz determines the dimensionality of the shared latent variable Z

  dx <- nrow(X)
  dy <- nrow(Y)
  if (dz > min(dx,dy)) {stop("dz < min(nrow(X),nrow(Y)) required!")}
  
  #  Dependency model
  #  X ~ N(Wx*z, sigmax*I)
  #  y ~ N(Wy*z, sigmay*I)
  #  i.e. isotropic marginals but in general  sigmax != sigmay
  # This is a special case of probabilistic factor analysis model
  # X ~ N(W*z, Sigma), where Sigma is diagonal (but not necessarily isotropic).

  # FIXME: ensure that X, Y have zero-mean (shift if necessary);
  # alternatively add mean parameter in the model

  # initialize
     inits <- initialize2(X,Y)
  phi.init <- inits$phi
      Dcov <- inits$Dcov
       Dim <- inits$Dim
    W.init <- inits$W
         W <- as.matrix(W.init$total[,1:dz])
        Wx <- as.matrix(W[1:dx,1:dz])
        Wy <- as.matrix(W[-(1:dx),1:dz])
       phi <- list(X = 1, Y = 1)

  
  # iterate until convergence:
  while (delta>epsilon) {

          # print(delta)
          W.old <- W
          
          ##########################################

          # Update Phi
          phi$X <- update.phi.isotropic(Dcov$X, Wx, phi$X, Dim$X) 
          phi$Y <- update.phi.isotropic(Dcov$Y, Wy, phi$Y, Dim$Y)
          
          #######################################

          phi.inv.full <- diag(c(rep(1/phi$X, Dim$X), rep(1/phi$Y, Dim$Y)))
             M <- set.M.full(W, phi.inv.full, dz) # corresponds to G in Bishop's book
          beta <- set.beta.fullcov(M, W, phi.inv.full)
             W <- as.matrix(W.cca.EM(Dcov, M, beta))
            Wx <- as.matrix(W[1:dx,])
            Wy <- as.matrix(W[-(1:dx),])
           
          ########################################
          
          # check convergence (enough to check W)
          delta <- max(abs(as.vector(W-W.old)))
          
        }

  # retrieve W and phi as in other functions:
  W2 <- list()
  W2$total <- W
  W2$X <- as.matrix(W[1:dx,])
  W2$Y <- as.matrix(W[-(1:dx),])

  phiX <- diag(phi$X, nrow(X))
  rownames(phiX) <- rownames(X)
  colnames(phiX) <- rownames(phiX)
  phiY <- diag(phi$Y, nrow(Y))
  rownames(phiY) <- rownames(Y)
  colnames(phiY) <- rownames(phiY)
  phitotal <- diag(c(diag(phiX),diag(phiY)),(nrow(X)+nrow(Y)))
  rownames(phitotal) <- c(rownames(X),rownames(Y))
  colnames(phitotal) <- rownames(phitotal)
  phi <- list(X = phiX, Y = phiY, total = phitotal)

  list(W=W2, phi=phi)

}
