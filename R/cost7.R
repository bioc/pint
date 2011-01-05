cost7 <-
function (Wvec, phi, Dcov, Dim, priors) {

  if (!is.null(priors$W)) {
    Wvec <- abs(Wvec)
  }
    
  W <- get.W4(Wvec, Dim)$X
  
  # SimCCA - just use identical Wx = Wy
  wtw <- W%*%t(W)

  Sigma <- rbind(cbind(wtw + phi$X, wtw),
          cbind(wtw,wtw + phi$Y))

        # Marginal cost for the whole data set
        # integrated over z
        # given parameters W, phi
        # P(X,Y | W, phi) = integral N(X|Wx*z,phix)*N(Y|Wy*z,phiy)*N(z|0,I)
        # We report -logP here

        # restrict solutions to cases where det(Sigma)>=0

  # -logP for the data
  detsigma <- det( Sigma )
  if (detsigma > 0) {
    cost.data <- log(detsigma) + sum(diag(solve(Sigma)%*%Dcov$total))
  } else { cost.data <- Inf }
  
  # -logP for W prior
  # wcost <- sum((W$X)^2) * priors$W
  # Assuming exponential prior distribution with rate parameter priors$W
  vec <- abs(as.vector(W))
  if (!is.null(priors)) {
    #multiply by 2 to count for both wx and wy
    wcost <- -2*sum(dexp(vec, rate = priors$W, log = TRUE))

  } else {
    wcost <- 0 # no effect
  }

  #print(paste("wcost", wcost))
  #print(paste("cost.data", cost.data))
  
  cost.data + wcost

}

