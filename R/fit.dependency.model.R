fit.dependency.model <-
function (X, Y,
          zDimension = 1,
          marginalCovariances = "full",
          H = 1, sigmas = 0, covLimit = 1e-3,
          mySeed = 123, priors = list(), matched = TRUE)
{

  # (C) Olli-Pekka Huovilainen and Leo Lahti 
  # FreeBSD license (keep this notice).

  #zDimension = 1; marginalCovariances = "full"; H = 1; sigmas = 0; covLimit = 1e-3; mySeed = 123; priors = NULL
  
  # Fits the generative model
  # X = Wx * z + epsx
  # Y = Wy * z + epsy
  # with various modeling assumptions



  
  # FIXME using priors$sigma.w and the parameter sigmas for the same purpose
  # (in !matched case always sigma.w)
  # FIXME: later completely replace 'sigmas' by 'priors$sigma.w'
  # to make more explicit notation for prior parameters
  if (is.null(priors$sigma.w)) {
    if (exists("sigmas")) {
      priors$sigma.w <- sigmas
    } else {
      priors$sigma.w <- 0
    }
  } else if (!priors$sigma.w == sigmas) {
    warning("Inconsistent priors, priors$sigma.w is overriding sigmas")
    # priors$sigma.w and sigmas are alternative parameters, referring to the same prior.
    priors$sigma.w <- sigmas
  }    

  if (covLimit == 0)  {covLimit <- 1e-3} # avoid numerical overflows

  # Center data
  X <- t(centerData(t(X), rm.na = TRUE))
  Y <- t(centerData(t(Y), rm.na = TRUE))

  # Check dimensionality
  if(zDimension > nrow(X) || zDimension > nrow(Y)) {
    message("zDimension exceeds data dimensionality; using full dimensionality min(ncol(X), ncol(Y)).")
    zDimension <- min(nrow(X), nrow(Y))
  }       

  if (nrow(X) < nrow(Y)) {stop("If the two data matrices do not have equal dimensionality, then place the smaller one in Y.")} # FIXME automate
  
  # Storage for results
  res <- NA; method <- ""

  if (!nrow(X) == nrow(Y)) {
    message("The data sets have unequal dimensionality.")
    if (matched) {stop("Cannot use matched methods for nonmatched data.")}
    matched <- FALSE
    message("- Applying the methods for nonmatched data.")
  }

  # Matrix normal distribution mean matrix not specified -> no
  # effective priors here.
  if (is.null(H)) {
    H <- 1; priors$sigma.w <- Inf;
    warning("The matrix H is not specified. Assuming priors$sigma.w is infinite and Wx ~ Wy unconstrained.")
  }
  
  #####################################

  # FIXME: replace this if clause by user options of
  # matched/nonmatched data
  if (!matched) {
    #message("Assuming non-matched variables.")

    # no similarity prior for Wx ~ Wy, completely uncoupled
    # -> no matching, no similarities (TODO: some similarities might still apply, implement later?)
    priors$sigma.w <- Inf

    if (!zDimension == 1) {
      warning("For non-matched variables only 1-dimensional latent variable model has been implemented. Using zDimension = 1.") # FIXME: implement multidimensional cases for nonmatched variables. Some if it may rock already; test.
      zDimension <- 1
    }

    if (!is.null(priors$sigma.w)) {
      warning("Non-matched case (matched = FALSE). Similarity between Wx, Wy is not constrained by priors in this case. Ignoring priors for Wx~Wy relation (priors sigma.w).")
      priors$sigma.w <- NULL
    }

    # FIXME: implement here the fast non-prior case a.k.a pCCA when no
    # priors are provided?
    if (is.null(priors$W)) {
      stop("Provide priors$W with nonmatched variables!")
    }

    if (marginalCovariances == "full") {

      # Prior for W is given -> need to optimize W (no analytical
      # solution to EM)
      # Wx ~ Wy free and priors for W given    
      # priors$W is the rate parameter of the exponential.
      # The smaller, the flatter
      res <- simCCA.optimize3(X, Y, zDimension, mySeed = mySeed, epsilon = covLimit, priors = priors)
    } else {
      stop("Only the case marginalCovariances = full has been implemented for nonmatched variables.")
    }

    } else if (matched) {

      #message("Assuming matched variables.")
      # Matched case (for instance, for non-segmented data)
      
      if (any(is.na(H))) {

        warning("H cannot contain NAs! Using nonconstrained version with priors$sigma.w = Inf.")

        H <- 1
        priors$sigma.w <- Inf
      }

      if (priors$sigma.w == 0 && !marginalCovariances == "full") {
        stop("With priors$sigma.w = 0 the matched simcca model is implemented only with full marginal covariances.")
      }    

      # Case I : Relation Wx, Wy not constrained.
      # FIXME: define this in a more elegant manner through priors
      if (priors$sigma.w == Inf) {
        wx.wy.relation.constrained <- FALSE
      } else {
        wx.wy.relation.constrained <- TRUE
      }

      # Exponential prior for elements in W. Forces nonnegative solutions
      if (!is.null(priors$W)) {
        nonnegative.w <- TRUE
        #if (!priors$W > 0) { stop("Nonnegative prior for priors$W required!") }
      } else {
        nonnegative.w <- FALSE
      }
      
      if ( !wx.wy.relation.constrained ) {

        #message("Relation between Wx, Wy is unconstrained.")
        
        if (!nonnegative.w) {
         # message("Wx ~ Wy free. No regularization for W.")

          if (marginalCovariances == "full") {           
            # Analytical ML solution for CCA 
            res <- calc.pcca(X, Y, zDimension)            
            method <- "pCCA"
            message("Full marginal covariances.")
          } else if (marginalCovariances == "diagonal"){                                  
            # Probabilistic factor analysis model as proposed in          
            # EM Algorithms for ML Factoral Analysis, Rubin D. and      
	    # Thayer D. 1982    
            res <- calc.pfa(X, Y, zDimension)    
            method <- "pFA"
            message("Diagonal marginal covariances.")
          } else if (marginalCovariances == "isotropic") {
	    # pCCA assuming isotropic margins with phiX != phiY 
            res <- calc.pcca.with.isotropic.margins(X, Y, zDimension, epsilon = covLimit)  
            method <- "pFA with isotropic margins"
            message("Isotropic marginal covariances.")          
          } else if(marginalCovariances == "identical isotropic"){        
            res <- calc.ppca(X, Y, zDimension)
            message("Identical isotropic marginal covariances.")
            method <- "pPCA"
          }                         
        } else if (nonnegative.w) {

          message("Wx ~ Wy free. Regularized W (W>=0).")        

         if (!marginalCovariances == "full") {           
           warning("Only full marginal covariances implemented with uncontrained Wx/Wy and regularized W. Using full marginal covariances.")
           marginalCovariances <- "full"           
         }
          
          # Currently implemented exponential prior for W,
          # priors$W is the rate parameter of the exponential.
          # The smaller, the flatter. 
          # By default, this function does not constrain Wx~Wy
          # (priors$sigma.w = Inf) but imposes prior on W.  FIXME:
          # this should work also with constrained Wx ~ Wy, test and
          # compare FIXME: check and message how full/diag/isotropic
          # cov go here

          # any zDimension should work here
          res <- simCCA.optimize2(X, Y, zDimension, mySeed = mySeed,
                                 epsilon = covLimit, priors = priors)
        }

      } else if ( wx.wy.relation.constrained ){

        method <- "pSimCCA"
        #message("Regulating the relationship Wx ~ Wy.")
        
        # Case IIa: fully constrained case Wx = Wy
        if (priors$sigma.w == 0) { #Wx = Wy        
        
          #  SimCCA with full covariances with constraint Wx = Wy
          #  "probsimcca.full" = aucs.simcca.full      
          #  Denoting Wy = T*Wx = TW; W = Wx this fits the case T = I with
          #  full-rank Wx, Wy, Sigmax, Sigmay: (dxd-matrices where d equals to
          #  number of features in X and Y)

          # If prior for W is given, we must optimize W (no analytical
          # solution to EM)
          
          # Regularization (W > 0)
          if (nonnegative.w) {
            
            # SimCCA Wx = Wy with regularized W (W>=0)
            #message("Case Wx = Wy and regularized W.")

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
            res <- optimize.simCCA.W2(W.init$X,
                                      phi.init,
                                      Dim = Dim,
                                      Dcov = Dcov,
                                      nullmat = nullmat,
                                      epsilon = covLimit,
                                      par.change = 1e6,
                                      mySeed = mySeed + 1,
                                      dz = zDimension,
                                      priors = priors)

          } else if (!nonnegative.w) {
            # mlsp'09 simcca
            #message("Case Wx = Wy. No regularization for W.")
            res <- simCCA.optimize.fullcov.EM(X, Y, zDimension, mySeed = mySeed, epsilon = covLimit)
            # FIXME: priors for phi todo?            
            
          }
        } else if (priors$sigma.w != 0) {
          # Case IIb: partially constrained Wx ~ Wy
          if (nonnegative.w) {
            stop("Not implemented regularized (nonnegative) W with partially constrained Wx ~ Wy. Only special cases priors$sigma.w = 0 and priors$sigma.w = Inf are available.")            
          } else if (!nonnegative.w) {
            #message("Partially constrained Wx ~ Wy. No regularization for W.")
            
            if (marginalCovariances == 'isotropic') {

              # Make H identity matrix if scalar is given
              if(length(H) == 1){ H <- diag(1, nrow(X), nrow(Y)) }
              if(ncol(H) != nrow(X)){ stop("columns of H must match rows of X") }
              if(nrow(H) != nrow(Y)){ stop("rows of H must match rows of Y") }
                 
              #message("SimCCA with isotropic covariances and regularized H (through sigmas).")
              res <- simCCA.optimize(X, Y, zDimension, H, sigma2.T = sigmas, sigma2.W = 1e12, mySeed, epsilon = covLimit)
            } else if (!marginalCovariances == 'isotropic') {
              stop("With priors$sigma.w !=0, only isotropic marginal covariances implemented.")
            }
          }
        } 
      }
    }

  ##################################################################

  # Test whether model exists for given arguments
  if (any(is.na(res))) {
    stop("Error with model parameters.")
  } else {
    params <- list(marginalCovariances = marginalCovariances, sigmas = sigmas, H = H, zDimension = zDimension, covLimit = covLimit)
    
    score <- dependency.score(res)
    geneName <- rownames(X)[[ trunc((nrow(X) + 1)/2) ]]
  }
  
  if(is.null(geneName)) {geneName <- ""}
  model <- new("DependencyModel", W = res$W, phi = res$phi, score = score, chromosome = "", arm = "", windowSize = nrow(Y), method = method, params = params, geneName = geneName)	
  
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
          
