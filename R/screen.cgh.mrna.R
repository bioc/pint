screen.cgh.mrna <- function(X, Y, windowSize = NULL,
                            chromosome,
                            arm,
                            method = "pSimCCA",
                            params = list(),
                            max.dist = 1e7, 
                            outputType = "models",
                            useSegmentedData = TRUE,
                            match.probes = TRUE,
                            regularized = FALSE)
{

#X <- ge; Y <- cn; windowSize = NULL; method = "pSimCCA"; params = list(); max.dist = 1e7; outputType = "models"; useSegmentedData = FALSE; match.probes = FALSE; regularized = FALSE 

  # FIXME: quick hack - later modify genomeModels class
  if (is.null(X$info$arm) && is.null(Y$info$arm)) {
    warning("Arm information missing, artificially adding p arm for all probes.")
    X$info$arm <- rep("p", nrow(X$info))
    Y$info$arm <- rep("p", nrow(Y$info))
  }

  if (is.null(X$info$arm) && !is.null(Y$info$arm) && nrow(X$info) == nrow(Y$info)) {
    warning("Arm information missing from X data, borrowing the arm info from Y data.")
    X$info$arm <- Y$info$arm
  }

  if (!is.null(X$info$arm) && is.null(Y$info$arm) && nrow(X$info) == nrow(Y$info)) {
    warning("Arm information missing from Y data, borrowing the arm info from X data.")
    Y$info$arm <- X$info$arm
  }    

  if (is.null(windowSize)) {
    windowSize <- min(floor(ncol(X$data)/3),15)
    cat(paste("Chromosomal window (windowSize) not specified. Using
default ratio of 1/3 between features and samples (with max window
size 15 probes). Using window size ", windowSize,"\n"))
  }
      
  #FIXME: move all preprocessing stuff to dedicated preprocessing functions
  # pint.data and pint.match
  
  # Check ordering of samples
  if (any(colnames(X$data) != colnames(Y$data))) {
    warning("Samples not in the same order in the two data sets. Using samples that are found in both data sets and reordering the samples.\n")

    commons <- intersect(colnames(X$data), colnames(Y$data))
    if (length(commons) > 1) {
      X$data <- X$data[, commons]
      Y$data <- Y$data[, commons]
    } else {
      stop("Not enough common samples found. Check that the corresponding samples in the two data sets have identical names.\n")
    }
  }

  ## Checks that segmented data is not used when not implicitely
  ## indicated by argument
  ## FIXME: this is independent of 'segmented' option, join these later
  if (!useSegmentedData){
    if (test.segmented(X$data) || test.segmented(Y$data)){
      warning("Segmented data found while method for non-segmented data is selected.\n", immediate. = TRUE)    
    }
  }

  if (match.probes) {   # FIXME: change name to match.probes or something 
    # Match probes
    tmp <- pint.match(X, Y, max.dist, useSegmentedData = useSegmentedData)
    X <- tmp$X
    Y <- tmp$Y
    match.probes <- FALSE # now the probes are matched.
  } else {
    # If user claims that no matching is needed
    # this will implicate that matching has already been
    # performed. Verify: the number of probes should be
    # identical in ge and cn data sets
    if (!nrow(X$data) == nrow(Y$data)) {
      stop("If match.probes == FALSE then the number of probes in ge and cn data sets (X$data, Y$data) need to match!")      
    }  
  }
  
  # Remove probes where observations are not available in either data set
  # TODO

  ############################################################################

  # Set priors

  # Currently imlement only pSimCCA for the screening
  if (method == "pSimCCA") {
        # similarity prior: Wx = Wy identically
	priors <- list(sigma.w = 0) 
  } else if (method == "nonmatched") {
    # no similarity prior for Wx ~ Wy, completely uncoupled 
    priors <- list(sigma.w = Inf)
  } else {
    message("Method not among pSimCCA, nonmatched. Using empty prior.")
    priors <- NULL
  }
	      
  if (regularized) {priors$W <- 1e-3} # W>=0 prior
		  
  ###########################################################################

  if (method == "pSimCCA") {
    if (is.null(params$H)) {
      params$H <- diag(1, windowSize, windowSize)
    }
  } else if (method == "pPCA" || method == "pCCA" || method == "pFA") {
    params$H <- NA
  } else {
    if (is.null(params$H))
      params$H <- diag(1, windowSize, windowSize)
  }

  if (is.null(params$sigmas))
    params$sigmas <- 0
  if (method == "pPCA") {
    params$marginalCovariances <- "identical isotropic"
  } else if (method == "pFA") {
    params$marginalCovariances <- "diagonal"
  } else if (method == "pCCA") {
    if (is.null(params$marginalCovariances)) {
      params$marginalCovariances <- "full"
    }
  } else {
    if (is.null(params$marginalCovariances)) {
      if (params$sigmas == 0) {
        params$marginalCovariances <- "full"
      } else {
        params$marginalCovariances <- "isotropic"
      }
    }
  }
  if (is.null(params$zDimension))
    params$zDimension <- 1
  if (!is.null(params$H) && any(is.na(params$H))) {
    if (params$marginalCovariances == "full")
      method = "pCCA"
    if (params$marginalCovariances == "isotropic")
      method = "pCCA"
    if (params$marginalCovariances == "diagonal")
      method = "pFA"
    if (params$marginalCovariances == "identical isotropic")
      method = "pPCA"
  } else {
    method = "pSimCCA"
  }

  # Convert X and Y chromosomes to numbers
  if (!missing(chromosome)){
    if (chromosome == 23) chromosome <- "X"
    if (chromosome == 24) chromosome <- "Y"
  }

  # give warning if arm param given and no arm data is available
  if (!missing(arm) && any(X$info$arm == "")){
    warning("No arm information in the data. Calculating thw whole chromosome")
    arm <- NULL
  }

  if (missing(chromosome)) {
    models <- calculate.genome(X, Y, windowSize, method, params, match.probes = match.probes, priors = priors, regularized = regularized)
  } else if (missing(arm) || is.null(arm)) {
    models <- calculate.chr(X, Y, windowSize, chromosome, method, params, match.probes = match.probes, priors = priors, regularized = regularized)
  } else {
    models <- calculate.arm(X, Y, windowSize, chromosome, arm, method, params, match.probes = match.probes, priors = priors, regularized = regularized)
  }

  # TODO: move this when segmented method is fully implemented (option 'match.probes')
  models@params <- c(models@params, segmentedData = useSegmentedData)

  if(outputType == "data.frame"){
    message("Convert model to data.frame")
    return(as.data.frame(models))
  }
  else {
    return(models)
  }
}

