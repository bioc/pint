screen.cgh.mir <- function(X, Y, windowSize, chromosome, arm, method = "", params = list(), outputType = "models"){

  # Check ordering of samples
  if (any(colnames(X$data) != colnames(Y$data)))
    stop("Samples are not in the same order in both datas")


  ### Check parameters and put defaults where needed ###

  warningText <- FALSE

  if (method == "pSimCCA")
    warningText <- TRUE

  if(!is.null(params$H)){
    if(!is.na(params$H)){
      warningText <- TRUE
    }
  }
  params$H <- NA


  # sigmas
  if(!is.null(params$sigmas)){
    if(params$sigmas != 0){
      warningText <- TRUE
    }
  }
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
    params$marginalCovariances <- "full"
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
  if (params$marginalCovariances == "full")
    method = "pCCA"
  if (params$marginalCovariances == "isotropic")
    method = "pCCA"
  if (params$marginalCovariances == "diagonal")
    method = "pFA"
  if (params$marginalCovariances == "identical isotropic")
   	method = "pPCA"


  if(warningText)
    warning("Similarity constrained CCA cannot be used, method changed to CCA")

	# Convert chromosome argument to factor with correct levels
  if (!missing(chromosome)){
    if (chromosome == 23) chromosome <- "X"
    if (chromosome == 24) chromosome <- "Y"
    chromosome = factor(chromosome, levels = c(1:22,"X","Y"))
    if (is.na(chromosome))
      stop("Incorrect chromosome given.")
  }

  # Calculate dependency models
  if (missing(chromosome))
    models <- calculate.genome.sparse(X, Y, windowSize, method, params)
	
  else if (missing(arm))
    models <- calculate.chr.sparse(X, Y, windowSize, chromosome, method, params)
		
  else
    models <- calculate.arm.sparse(X, Y, windowSize, chromosome, arm, method, params)


  if(outputType == "data.frame"){
    return(as.data.frame(models))
  }
  else {
    return(models)
  }
}