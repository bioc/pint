fit.dependency.model <-
function (X, Y, zDimension = 1,
          marginalCovariances = "full", H = 1,
          sigmas = 0, covLimit = 0, mySeed = 123){

	# Center data
	X <- t(centerData(t(X)))
	Y <- t(centerData(t(Y)))

	# Check if dimensionality is too big
	if(zDimension > ncol(X) || zDimension > ncol(Y))
		stop("Dimension of latent variable too big")

	# Strogae for results from particular function that calculates dependency model
	res = NA
	
	method <- ""
	if (any(is.na(H))) {
		method <- "pCCA"
		if (marginalCovariances == "full") {
			
			# ML solution for CCA
			# Analytical solution!
			res <- pcca(X, Y, zDimension)
		}
		else if (marginalCovariances == "diagonal"){
 		        # Probabilistic factorial analysis model as proposed in
     			# EM Algorithms for ML Factoral Analysis, Rubin D. and
     			# Thayer D. 1982
			res <- pfa(X, Y, zDimension)			
			method <- "pFA"
		}
		else if (marginalCovariances == "isotropic") {
			if (covLimit == 0) 
				covLimit <- 1e-6
			# pCCA assuming isotropic margins
			# with phiX != phiY
			res <- pcca.with.isotropic.margins(X, Y, zDimension, epsilon=covLimit)
		}
		else if(marginalCovariances == "identical isotropic"){
			method <- "pPCA"
			#pPCA
			res <- ppca(X, Y, zDimension) 
		}
	}

	else {

		method <- "pSimCCA"
		if (sigmas == 0 && marginalCovariances == "full") {
			# SimCCA with full covariances
			# with constraint Wx = Wy
			# EM algorithm
			# "probsimcca.full" = aucs.simcca.full
	
			#  Denoting Wy = T*Wx = TW; W=Wx this fits the case T = I with
			#  full-rank Wx, Wy, Sigmax, Sigmay: (dxd-matrices where d equals to
			#  number of features in X and Y)
			H <- diag(1,nrow(X),nrow(Y))			
			if (covLimit == 0) 
				covLimit <- 1e-3		
			res <- simCCA.optimize.fullcov.EM(X, Y, zDimension, mySeed=mySeed, epsilon=covLimit)
		}
		else if (marginalCovariances == 'isotropic' && sigmas != 0) {
		        # Make H indetity matrix if scalar is given
			if(length(H) == 1){
		      		H <- diag(1,nrow(X),nrow(Y))			
			}
		     	
			if(ncol(H) != nrow(X)){
			      stop("columns of H must match rows of X")
			}
			if(nrow(H) != nrow(Y)){
			      stop("rows of H must match rows of Y")
			}

			# SimCCA with isotropic covariances and possibility to tune prior for T
			if (covLimit == 0) 
				covLimit <- 1e-6
			res <- simCCA.optimize(X, Y, zDimension, H, sigma2.T=sigmas, sigma2.W = 1e12, mySeed,epsilon=covLimit)
		}
	}
	
	# Test whether model exists for given arguments
	if (any(is.na(res))) {
		stop("Error with model parameters")
	}
	else {
		params <- list(marginalCovariances = marginalCovariances, sigmas = sigmas, H = H, 
		       	       zDimension = zDimension, covLimit = covLimit)
		score <- dependency.score(res)
		
		geneName <- dimnames(X)[[1]][ trunc((nrow(X)+1)/2) ]
		if(is.null(geneName))
			geneName = ""

		model <- new("DependencyModel", W = res$W, phi = res$phi, score = score, 
					windowSize = dim(X)[1], method = method, params = params, geneName = geneName)	
	}
	model
}

