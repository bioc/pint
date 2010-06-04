cost5 <-
function (W, phi, Dcov) {

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
  detsigma <- det( Sigma )
  if (detsigma > 0) {
                #print(detsigma)
    log(detsigma) + sum(diag(solve(Sigma)%*%Dcov$total))
  } else { Inf }

}

