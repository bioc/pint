set.beta <-
function (M, W, phi) {
	# assuming isotropic marginal covariance
	M%*%t(W)/phi
}

