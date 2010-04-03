screen.cgh.mrna <- function(X, Y, windowSize, chromosome, arm, method = "", params = list()){

	# Check ordering of samples
	if (any(colnames(X$data) != colnames(Y$data)))
		stop("Samples are not in the same order in both datas")


	### Check parameters and put defaults where needed ###

	# H (priori for T in Wy = T*Wx)
	else if (method == "pSimCCA") {
		if (is.null(params$H)) {
        	params$H <- diag(1,windowSize,windowSize)
		}
	}
	else if (method == "pPCA" || method == "pCCA" || method == "pFA") {
	        params$H <- NA
	}
	else {
		if(is.null(params$H))
			params$H <- diag(1,windowSize,windowSize)
	}

	# sigmas
	if(is.null(params$sigmas))
			params$sigmas <- 0

	# Marginal covariances
	if (method == "pPCA") {		
		params$marginalCovariances <- "identical isotropic"
	}
	else if (method == "pFA") {		
		params$marginalCovariances <- "diagonal"
	}
	else if (method == "pCCA") {
     	if (is.null(params$marginalCovariances)) {
			params$marginalCovariances <- "full"
		}
	}
	else {
		if (is.null(params$marginalCovariances)) {			
			if (params$sigmas == 0) {
				params$marginalCovariances <- "full"
			}
			else {
				params$marginalCovariances <- "isotropic"
			}
		}
	}

	# Dimension of z
	if (is.null(params$zDimension)) 
		params$zDimension <- 1

	# Limit for convergence
	if (is.null(params$covLimit)) 
		params$covLimit <- 0
	if (is.null(params$mySeed)) 
		params$mySeed <- 566

	# Set method name
	if (any(is.na(params$H))){
	    if (params$marginalCovariances == "full")
		   	method = "pCCA"
		if (params$marginalCovariances == "isotropic")
		   	method = "pCCA"
		if (params$marginalCovariances == "diagonal")
		   	method = "pFA"
		if (params$marginalCovariances == "identical isotropic")
		   	method = "pPCA"
	}	
	else {
	   	method = "pSimCCA"
	}

	# Calculate dependency models
	if (missing(chromosome))
		models <- calculate.genome(X, Y, windowSize, method, params)
	
	else if (missing(arm))
		models <- calculate.chr(X, Y, windowSize, chromosome, method, params)
		
	else
		models <- calculate.arm(X, Y, windowSize, chromosome, arm, method, params)
		
	
	return(models)
}