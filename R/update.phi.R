update.phi <-
function (Dcov, M, beta, W, phi) {

	# This assumes simple marginal covariance: sigma2 * I
	# i.e. just one parameter.

	# Mean equals to dividing with dimension
	phix = mean(diag(Dcov$X - Dcov$X%*%t(beta$X)%*%t(W$X)))
	phiy = mean(diag(Dcov$Y - Dcov$Y%*%t(beta$Y)%*%t(W$Y)))

 	# Return phi
        list(X = phix, Y = phiy)
}

