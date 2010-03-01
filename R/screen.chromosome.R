screen.chromosome <- function(X, Y, windowSize, chromosome, arm, method = "", params = list()){

	# Check ordering of samples
	if (any(colnames(X$data) != colnames(Y$data)))
		stop("Samples are not in the same order in both datas")


	### Check parameters and put defaults where needed ###

	# H (priori for T in Wy = T*Wx)
	if (method == "TPriorpSimCCA") {
		if (is.null(params$H)) {
                        params$H <- diag(1,windowSize,windowSize)
		}
	}
	else if (method == "pSimCCA") {
                params$H <- diag(1,windowSize,windowSize)
	}
	else {
		if(is.null(params$H))
			params$H <- NA
	}

	# Marginal covariances
	if (method == "pPCA") {		
		params$marginalCovariances <- "identical isotropic"
	}
	else if (method == "pFA") {		
		params$marginalCovariances <- "diagonal"
	}
	else if (method == "TPriorpSimCCA") {		
		params$marginalCovariances <- "isotropic"
	}
	else if (method == "pCCA") {
	     	params$marginalCovariances <- "full"
	}
	else {
		if(is.null(params$marginalCovariances)) 
			params$marginalCovariances <- "full"
	}

	# Dimension of z
	if (is.null(params$zDimension)) 
		params$zDimension <- 1

	# sigmas
	if (method == "TPriorpSimCCA") {
		if (is.null(params$sigmas)) {
			params$sigmas <- 1
		}
	}
	else {
		if(is.null(params$sigmas))
			params$sigmas <- 0
	}

	# Limit for convergence
	if (is.null(params$covLimit)) 
		params$covLimit <- 0
	if (is.null(params$mySeed)) 
		params$mySeed <- 566


	# Calculate dependency models
	if (missing(chromosome))
		models <- calculate.genome(X, Y, windowSize, method, params)
	
	else if (missing(arm))
		models <- calculate.chr(X, Y, windowSize, chromosome, method, params)
		
	else
		models <- calculate.arm(X, Y, windowSize, chromosome, arm, method, params)
		
	return(models)
}