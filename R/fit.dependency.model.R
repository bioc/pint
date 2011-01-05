fit.dependency.model <-
function (X, Y,
          zDimension = 1,
          marginalCovariances = "full",
          H = 1, sigmas = 0, covLimit = 1e-3,
          mySeed = 123, priors = NULL, version = "standard")
{

  if (covLimit == 0)  {covLimit <- 1e-3}

  # Center data
  X <- t(centerData(t(X), rm.na = TRUE)) #moi
  Y <- t(centerData(t(Y), rm.na = TRUE))

  # Check dimensionality
  if(zDimension > nrow(X) || zDimension > nrow(Y)) {
    message("zDimension exceeds dimensionality of X or Y; using full dimensionality min(ncol(X), ncol(Y)).")
    zDimension <- min(nrow(X), nrow(Y))
  }
       
  if (nrow(X) < nrow(Y)) {stop("If the two data matrices do not have equal dimensionality, then place the smaller one in Y.")}
  
  # Storage for results
  res <- NA; method <- ""

  #####################################

  if (version == "devel") {

    message(paste("in fit.dependency.model using version", version))

  if (!nrow(X) == nrow(Y)) {
    message("Dependency model without assuming matched variables (dimX != dimY).")

    # Assuming full marginal covariances.    
    # If prior for W is given, we must optimize W (no analytical solution to EM)
    # (add the fast non-prior case a.k.a pCCA later)

    if (!is.null(priors$W)) {

      # Case Wx ~ Wy free and priors for W given    
      # priors$W is the rate parameter of the exponential.
      # The smaller, the flatter

      res <- simCCA.optimize3(X, Y,
                              zDimension,
                              mySeed  = mySeed,
                              epsilon = covLimit,
                              priors  = priors
                              )
    } else {
      stop("Dependency model for non-matched features implemented only given priors$W")
    }
  } else {
    
    # Matched case (e.g. for non-segmented data)
        
    if (!is.null(H) && any(is.na(H))) {
      stop("H cannot contain NA values!")
    } else {

    method <- "pSimCCA"

    if (sigmas == 0 && marginalCovariances == "full") {
      # SimCCA with full covariances
      # with constraint Wx = Wy
      # EM algorithm
      # "probsimcca.full" = aucs.simcca.full
      
      #  Denoting Wy = T*Wx = TW; W = Wx this fits the case T = I with
      #  full-rank Wx, Wy, Sigmax, Sigmay: (dxd-matrices where d equals to
      #  number of features in X and Y)


      # If prior for W is given, we must optimize W (no analytical solution to EM)
      if (!is.null(priors$W) && priors$sigma.w == Inf) {

        #print("CCA-type model i.e. no constraints on Wx ~ Wy with regularized W (W>=0)")
        
        # Case Wx ~ Wy free and priors for W given
        
        # Currently implemented exponential prior for W,
        # priors$W is the rate parameter of the exponential.
        # The smaller, the flatter

        # By default, does not constrain Wx~Wy (priors$sigma.w = Inf) but instead imposes
        # prior on W
        res <- simCCA.optimize2(X, Y,
                                zDimension,
                                mySeed  = mySeed,
                                epsilon = covLimit,
                                priors  = priors
                               )

      } else if (!is.null(priors$W) && priors$sigma.w == 0) {
        
        # SimCCA Wx = Wy with regularized W (W>=0)
        
        #message("Case Wx = Wy and priors for W given")
        # Initialize
           inits <- initialize2(X, Y)
        phi.init <- inits$phi
          W.init <- inits$W
            Dcov <- inits$Dcov
             Dim <- inits$Dim
           Dim$Z <- zDimension
         nullmat <- inits$nullmat
        Nsamples <- inits$Nsamples
 
        # use this for full W (EM algorithm, unstable for n ~ p)
        res <- optimize.simCCA.W2(W.init$X, phi.init, Dim = Dim, Dcov = Dcov,
                                 nullmat = nullmat, epsilon = covLimit, par.change = 1e6,
                                 mySeed = mySeed + 1, dz = zDimension, priors = priors)

      } else { # at least case is.null(priors)
        # non-regularized case with Wx = Wy        
        res <- simCCA.optimize.fullcov.EM(X, Y, zDimension, mySeed = mySeed, epsilon = covLimit)
      }

    } else if (marginalCovariances == 'isotropic' && sigmas != 0) {
      # Make H indentity matrix if scalar is given
      if(length(H) == 1){ H <- diag(1, nrow(X), nrow(Y)) }      
      if(ncol(H) != nrow(X)){stop("columns of H must match rows of X")}
      if(nrow(H) != nrow(Y)){stop("rows of H must match rows of Y")}
        
      # SimCCA with isotropic covariances and possibility to tune prior for T
      res <- simCCA.optimize(X, Y,
                             zDimension,
                             H,
                             sigma2.T = sigmas,
                             sigma2.W = 1e12,
                             mySeed,
                             epsilon = covLimit)
    }
   }
  }	
  
  } else if (version == "standard") {
     
    if (!is.null(H) && any(is.na(H))) { 
      method <- "pCCA"                       
      if (marginalCovariances == "full") {                           
        # ML solution for CCA  
        # Analytical solution!                      
	res <- calc.pcca(X, Y, zDimension)            
      } else if (marginalCovariances == "diagonal"){                                  
        # Probabilistic factor analysis model as proposed in          
        # EM Algorithms for ML Factoral Analysis, Rubin D. and      
	# Thayer D. 1982    
	res <- calc.pfa(X, Y, zDimension)    
	method <- "pFA"          
      } else if (marginalCovariances == "isotropic") {
	  # pCCA assuming isotropic margins with phiX != phiY 
	  res <- calc.pcca.with.isotropic.margins(X, Y, zDimension, epsilon = covLimit)  
      } else if(marginalCovariances == "identical isotropic"){        
      	method <- "pPCA"                                   
     	res <- calc.ppca(X, Y, zDimension)        
      }                                          
    } else {                                 
      method <- "pSimCCA"                    
      if (sigmas == 0 && marginalCovariances == "full") {      
        # SimCCA with full covariances             
        # with constraint Wx = Wy                       
        # EM algorithm
        # "probsimcca.full" = aucs.simcca.full                                                                                     
        #  Denoting Wy = T*Wx = TW; W=Wx this fits the case T = I with                                                         
        #  full-rank Wx, Wy, Sigmax, Sigmay: (dxd-matrices where d equals to         
        #  number of features in X and Y)      
	H <- diag(1, nrow(X), nrow(Y))   
	res <- simCCA.optimize.fullcov.EM(X, Y, zDimension, mySeed = mySeed, epsilon = covLimit)                                                                
     } else if (marginalCovariances == 'isotropic' && sigmas != 0) {
       # Make H indetity matrix if scalar is given
       if(length(H) == 1){ H <- diag(1, nrow(X), nrow(Y)) }                               
       if(ncol(H) != nrow(X)){ stop("columns of H must match rows of X") }
       if(nrow(H) != nrow(Y)){ stop("rows of H must match rows of Y") }
	         
       # SimCCA with isotropic covariances and possibility to tune prior for T                                      
       res <- simCCA.optimize(X, Y, zDimension, H, sigma2.T = sigmas, sigma2.W = 1e12,                                         
                             mySeed, epsilon = covLimit)                               
    }  
  }        
                                                            
  }

  ##################################################################

   # Test whether model exists for given arguments
  if (any(is.na(res))) {
    stop("Error with model parameters")
  } else {
    params <- list(marginalCovariances = marginalCovariances, sigmas = sigmas, H = H, 
                   zDimension = zDimension, covLimit = covLimit)
    score <- dependency.score(res)
    geneName <- rownames(X)[[ trunc((nrow(X) + 1)/2) ]]
  }
  
  if(is.null(geneName)) {geneName <- ""}
  model <- new("DependencyModel", W = res$W, phi = res$phi, score = score, chromosome = "", arm = "",
                windowSize = nrow(Y), method = method, params = params, geneName = geneName)	

  model
}


ppca <- function(X, Y = NULL, zDimension = 1){                
  if (!is.null(Y)) {                 
    fit.dependency.model(X,Y,zDimension,marginalCovariances = "identical isotropic", H = NA, sigmas = 0)        
  } else {         
    if (ncol(X) > 1)                         
      X <- t(centerData(t(X), rm.na = TRUE))               
           
    # Check if dimensionality is too big               
    if(zDimension > ncol(X))                         
      stop("Dimension of latent variable too big") 
                      
    res <- calc.ppca(X, Y, zDimension) 
    method <- "pPCA"    
                               
    params <- list(marginalCovariances = "isotropic", sigmas = 0, H = NA,                               
                           zDimension = zDimension, covLimit = 0)      
    score <- dependency.score(res) 
    geneName <- ""     
    model <- new("DependencyModel",W = res$W, phi = res$phi, score = score, 
    	     	 chromosome = "", arm = "",    
                 windowSize = dim(X)[1], method = method, 
		 params = params, geneName = geneName)   
    model      
  }        
}    

               
	       
pfa <- function(X, Y = NULL, zDimension = 1){       
  if (!is.null(Y)) {         
    fit.dependency.model(X, Y, zDimension, marginalCovariances = "diagonal", H = NA, sigmas = 0)       
  } else {         
    if (ncol(X) > 1)                         
      X <- t(centerData(t(X), rm.na = TRUE))               
           
    # Check if dimensionality is too big               
    if(zDimension > ncol(X))                         
      stop("Dimension of latent variable too big") 
                      
    res <- calc.pfa(X, Y, zDimension)        
    method <- "pFA"
                               
    params <- list(marginalCovariances = "diagonal", sigmas = 0, H = NA,                               
                   zDimension = zDimension, covLimit = 0)      
    score <- dependency.score(res) 
    geneName <- ""     
    model <- new("DependencyModel", W = res$W, phi = res$phi, score = score, 
    	     chromosome = "", arm = "", windowSize = dim(X)[1], method = method, 
	     params = params, geneName = geneName)                                                                    
  }                     
}                                                                     
                            
pcca.isotropic <- function(X, Y, zDimension = 1, covLimit = 1e-6){
  fit.dependency.model(X,Y,zDimension,marginalCovariances = "isotropic", H = NA, sigmas = 0, covLimit = 1e-6)          
}                                                                      
                                                                   
pcca <- function(X, Y, zDimension = 1){  
  fit.dependency.model(X, Y, zDimension, marginalCovariances = "full", H = NA, sigmas = 0)                                                     
}                                                                      
          
