cost.W2 <-
function (vec, phi, priors, Dim, Dcov) {

  # remove sign as we assume W always positive here
  vec <- abs(vec)
  
  # Retrieve W from the parameter vector
  W <- get.W2(vec, Dim)

  # Marginal cost for the whole data set
  # integrated over z
  # given parameters W, phi
  # P(X,Y | W, phi) = integral N(X|Wx*z,phix)*N(Y|Wy*z,phiy)*N(z|0,I)
  # We report -logP here
	
  # Data prob. Taken from probCCA paper, section 4, l1
  wtw.xy <- W$X%*%t(W$Y)

  Sigma <- rbind(cbind(W$X%*%t(W$X) + phi$X, wtw.xy),
                 cbind(t(wtw.xy), W$Y%*%t(W$Y) + phi$Y))
  
  # -logP for the data
  cost.data <- log(det(Sigma)) + sum(diag(solve(Sigma)%*%Dcov$total))
  
  # -logP for W prior
  # wcost <- sum((W$X)^2) * priors$W
  # Assuming exponential prior distribution with rate parameter priors$W
  if (!is.null(priors)) {
    wcost <- -sum(dexp(vec, rate = priors$W, log = TRUE))
  } else {
    wcost <- 0 # no effect
  }
  
  cost.data + wcost

}

