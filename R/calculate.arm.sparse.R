calculate.arm.sparse <- function(X, Y, windowSize, chromosome, arm, method = "pSimCCA", params = list()){

	# Get Indices to X and Y for chosen chromosome and arm
	Xm <- pick.chr.arm(X,chromosome,arm)
	Ym <- pick.chr.arm(Y,chromosome,arm)

	# Storage for dependency scores
	scores <- vector()
	# Storage for gene window location
	locs <- vector()
	# Stroage for gene names
	genes <- vector()

	# method name
	methodName <- method
	if (method == 'TPriorpSimCCA') methodName = 'pSimCCA with T prior'
	message(paste("Calculating dependency models for ",chromosome,arm," with method ",methodName, 
		", window size:",windowSize,sep=""))
	
	modelList <- list()
	
	# index for modelList
	k <- 1
	
	if (length(Xm$info$loc) > 0) {
		for (n in seq_along(Xm$info$loc)) {

			# Get window fo dependency modeling
			window <- sparse.window(Xm,Ym,n,windowSize)
		
			# Skip windows that overlaps chromosome arms
			if (!window$fail){
						
				model <- fit.dependency.model(window$X, window$Y,zDimension = params$zDimension, 
					marginalCovariances = params$marginalCovariances, H = params$H, sigmas = params$sigmas, 
					covLimit = params$covLimit, mySeed=params$mySeed)
				setLoc(model) <- window$loc
				#setGeneName(res) <- window$geneName
				modelList[[k]] <- model
				k <- k+1
			}
		}
	}
	# Change chromosome and arm factors and get levels from X
	chromosome <- factor(chromosome, levels = levels(Xm$info$chr))
	arm <- factor(arm, levels = levels(Xm$info$arm))
	
	return(new("ChromosomeArmModels", models = modelList, chromosome = chromosome, arm = arm, windowSize = windowSize, 
		method = method, params = params))
}

