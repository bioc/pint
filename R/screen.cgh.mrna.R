screen.cgh.mrna <- function(X, Y, windowSize, chromosome, arm, method = "", params = list(), max.dist = 1e7)
{

  # Check ordering of samples
  if (any(colnames(X$data) != colnames(Y$data))) {
    warning("Samples not in the same order in the two data sets. Using samples that are found in both data sets and reordering the samples..")

    commons <- intersect(colnames(X$data), colnames(Y$data))
    if (length(commons) > 1) {
      X$data <- X$data[, commons]
      Y$data <- Y$data[, commons]
    } else {
      stop("Not enough common samples found. Check that the corresponding samples in the two data sets have identical names.")
    }
  }

  ## Check probes. We should have observations in each data set for the probes.

  # Match probes by location
  tmp <- pint.match(X, Y, max.dist)
    X <- tmp$X
    Y <- tmp$Y

  # Remove probes where observations are not available in either data set
  # TODO


  ############################################################################
  
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
  if (is.null(params$covLimit))
    params$covLimit <- 0
  if (is.null(params$mySeed))
    params$mySeed <- 566
  if (any(is.na(params$H))) {
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

  if (missing(chromosome)) {
    models <- calculate.genome(X, Y, windowSize, method, params)
  } else if (missing(arm)) {
    models <- calculate.chr(X, Y, windowSize, chromosome, method, params)
  } else {
    models <- calculate.arm(X, Y, windowSize, chromosome, arm, method, params)
  }


  return(models)
}

