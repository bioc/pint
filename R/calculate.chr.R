calculate.chr <- function(X, Y, windowSize, chromosome, method = "pSimCCA", params = list(), segmented, priors, regularized){

  # Check if arm information is missing
  if (is.null(X$info$arm)){

    return(calculate.arm(X, Y, windowSize, chromosome, method=method, params=params, segmented=segmented, priors=priors, regularized=regularized))
    
  } else {
    pArm <- calculate.arm(X, Y, windowSize, chromosome, 'p', method, params, segmented, priors, regularized)
    qArm <- calculate.arm(X, Y, windowSize, chromosome, 'q', method, params, segmented, priors, regularized)
    
 
    return(new("ChromosomeModels",
               models = append(pArm@models, qArm@models),
               chromosome = chromosome,
               method = method,
               params = params))
  }
}