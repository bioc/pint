cost.W <-
function (vec, phi, priors, Dim, Dcov, H) {

	# Retrieve the actual W and T from the parameter vector
 	wt <- get.W(vec, Dim)
	W <- wt$W
	T <- wt$T

	# Marginal cost for the whole data set
	# integrated over z
	# given parameters W, phi
	# P(X,Y | W, phi) = integral N(X|Wx*z,phix)*N(Y|Wy*z,phiy)*N(z|0,I)
	# We report -logP here
	
	# Data prob. Taken from probCCA paper, section 4, l1

	wtw.xy <- W$X%*%t(W$Y)

	Sigma <- rbind(cbind(W$X%*%t(W$X) + phi$X*diag(Dim$X), wtw.xy),
		      cbind(t(wtw.xy),W$Y%*%t(W$Y) + phi$Y*diag(Dim$Y)))

	# -logP for the data
	cost.data <- log(det(Sigma)) + sum(diag(solve(Sigma)%*%Dcov$total))

	# -logP for T prior
        tcost <- sum((T - H)^2) * priors$T
		
	# -logP for W prior - skip since not used now
	#wcost <- sum((W$X)^2) * priors$W

	cost.data + tcost #+ wcost

}

