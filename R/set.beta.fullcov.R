

set.beta.fullcov = function (M, W, phi.inv) {
	# assuming full marginal covariance
	# as in section 4.1 EM algorithm from BachJordan probCCA paper
	#phi.inv%*%W%*%M
	M%*%t(W)%*%phi.inv
}

