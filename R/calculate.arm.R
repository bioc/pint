calculate.arm <- function(X, Y, windowSize, chromosome, arm = NULL, method =
"pSimCCA", params = list(), match.probes = FALSE, priors = NULL, regularized = FALSE){

  # Get Indices to X and Y for chosen chromosome and arm
  Xm <- pick.chr.arm(X, chromosome, arm)
  Ym <- pick.chr.arm(Y, chromosome, arm)

  # Storage for dependency scores
  scores <- vector()

  # Storage for gene window location                                         
  locs <- vector()

  # Storage for gene names
  genes <- vector()

  # method name
  methodName <- method
  message(paste("Calculating dependency models for ",chromosome,arm," with method ",methodName, 
		", window size:",windowSize,sep=""))
	
  modelList <- list()
	
  # index for modelList
  k <- 1
	
 if (length(Ym$info$loc) > 0) {                           
   for (n in seq_along(Ym$info$loc)) {       
     message(chromosome, arm, "; window ", n, "/", length(Ym$info$loc))                  
     # Assuming Y is the copy number data (important when match.probes= TRUE)                                                              
        # Get window to dependency modeling                                
        if (!match.probes) {                          
	  #message("Using fixed chromosomal window size.")       
	  window <- fixed.window(Xm, Ym, n, windowSize)
	} else {                           
	  message("Selecting the closest expression probes for each copy number segment.")       
          #NOTE: this calculates for each in X, window in Y      
	  # therefore convert as we have copy number in Y        
	  tmp <- sparse.window(Ym, Xm, n, windowSize)  
	  window <- tmp          
	  window$X <- tmp$Y      
	  window$Y <- tmp$X      
	}
                                                   
	# Skip windows that overlaps chromosome arms                    
       if (!window$fail){                     
   
           if (!match.probes && !regularized) {               
	     #print(c("no-segment", "no-regu"))         
	     # still ensure that no regularization used here:                 
             priors$W <- NULL                           
             model <- fit.dependency.model(window$X,         
                                           window$Y,   
                                           zDimension = params$zDimension,
                                           marginalCovariances = params$marginalCovariances,      
                                           priors = list(Nm.wxwy.mean = params$H,          
                                           Nm.wxwy.sigma = params$sigmas),includeData = FALSE, calculateZ = FALSE)                                     
    
  
             } else if (!match.probes && regularized) {                     
	           #print(c("no-segment", "yes-regu"))             
		   # 1-dimensional cca (in general, Wx != Wy) with nonnegative W                 
                   # assuming matched probes                                
               model <- fit.dependency.model(window$X, window$Y, priors = priors)                 
                                                                 
            } else if (match.probes) {            
	      message("Note: regularization is used with segmented data.")                 
              regularized <- TRUE                          
	      # always use positive prior for W here     
	      if (is.null(priors$W)) {priors$W <- 1e-3} # uninformative                                  
	      	 model <- fit.dependency.model(window$X, window$Y,zDimension = params$zDimension, priors = priors)                                   
                                                                      
     	   }
           model <- as(model, "GeneDependencyModel")
           #setGeneName(model) <- rownames(window$X)[[ trunc((nrow(window$X) + 1)/2) ]]
           setGeneName(model) <- window$geneName
          #setLocs(model) <- window$loc        
          setLoc(model) <- window$loc
        setChromosome(model) <- chromosome
        if (!is.null(arm)) setArm(model)  <- arm          
        modelList[[k]] <- model
        k <- k + 1                                 
      }                
    }
  }
 #Change chromosome and arm factors and get levels from X          
 #chromosome <- factor(chromosome, levels = c(1:22,"X","Y"))     
 #arm <- factor(arm, levels = levels(X$info$arm))      
 
 return(new("ChromosomeModels",                                  
            models = modelList,                                   
            chromosome = chromosome,    
            method = method,                                                
            params = params))                           
}                     

