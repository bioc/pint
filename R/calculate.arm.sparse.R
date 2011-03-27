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
						
				             model <- fit.dependency.model(window$X,         
                                           window$Y,   
                                           zDimension = params$zDimension,
                                           marginalCovariances = params$marginalCovariances,      
                                           priors = list(Nm.wxwy.mean = params$H,          
                                           Nm.wxwy.sigma = params$sigmas),includeData = FALSE, calculateZ = FALSE)        
				setLoc(model) <- window$loc
				setChromosome(model) <- as.character(chromosome)
        setGeneName(model) <- window$geneName
        if (!is.null(arm)) setArm(model)  <- arm   
				modelList[[k]] <- model
				k <- k+1
			}
		}
	}
 
 return(new("ChromosomeModels",                                  
            models = modelList,                                   
            chromosome = chromosome, 
            method = method,                                                
            params = params))                           
} 


