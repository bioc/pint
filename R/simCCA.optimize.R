simCCA.optimize <-
function (X, Y, zDimension=1, H, sigma2.T = 1e-12, sigma2.W = 1e12, mySeed=123, epsilon = 1e-6) {

	#
	# sigma2.T: tuning parameter for T prior: large is uninformative i.e. ordinary CCA, small values lead to SimCCA
	# sigma2.W: tuning parameter for W prior: large is uninformative 
	# zDimension: assumed dimensionality for the shared latent variable z
	# H: priori for T

	set.seed(mySeed)

	#################################################
	# Initialize
	#################################################

	Nsamples = ncol(X)
	
	priors = list() # small means weak prior
	priors$T = 1/(2 * Nsamples * sigma2.T)
	priors$W = 1/(2 * sigma2.W * Nsamples)

	Dim = list()
	Dim$X = nrow(X)
	Dim$Y = nrow(Y)
	Dim$Z = zDimension

	Dcov = list()
	Dcov$X = cov(t(X))
	Dcov$Y = cov(t(Y))
	#Dcov$xy = cov(t(X),t(Y))
	#Dcov$yx = cov(t(Y),t(X))
	Dcov$total = cov(t(rbind(X,Y)))

	# Initialize 
	W.init = list()
	W.init$X = as.matrix(eigen(Dcov$X)$vectors[,1:Dim$Z])
	T.init = array(rnorm(Dim$Y*Dim$X,0,sd(W.init$X)),dim=c(Dim$Y,Dim$X))
	#T.init = diag(rnorm(min(Dim$X,Dim$Y),0,sd(W.init$X)),Dim$Y,Dim$X)
	W.init$Y = T.init%*%W.init$X
	
	# Assume scalar diagonal noise on the marginals (shared by all features)
	phi.init = list(X=var(as.vector(X)),Y=var(as.vector(Y)))

	##################################################
	# optimize until convergence
	##################################################

	res = optimize.W(W.init, T.init, phi.init, Dim, Dcov, priors, H, epsilon, par.change = 1e6, cost.old = 1e6, mySeed=mySeed+1)

	W <- list(X = res$W$X, Y = res$W$Y, total = rbind(res$W$X,res$W$Y))
	rownames(W$X) <- rownames(X)
	rownames(W$Y) <- rownames(Y)
	rownames(W$total) <- c(rownames(X),rownames(Y))

	# Make phis diagonal matrices
  	phiX <- diag(res$phi$X, nrow(X))
  	rownames(phiX) <- rownames(X)
  	colnames(phiX) <- rownames(phiX)
  	phiY <- diag(res$phi$Y, nrow(Y))
  	rownames(phiY) <- rownames(Y)
  	colnames(phiY) <- rownames(phiY)
  	phitotal <- diag(c(diag(phiX),diag(phiY)),(nrow(X)+nrow(Y)))
  	rownames(phitotal) <- c(rownames(X),rownames(Y))
  	colnames(phitotal) <- rownames(phitotal)
  	phi <- list(X = phiX, Y = phiY, total = phitotal)

	return(list(W = W, phi = phi))

}

