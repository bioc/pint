calculate.chr.sparse <- function(X, Y, windowSize, chromosome, method = "pSimCCA", params = list()){

  # Check if arm information is missing
  if (is.null(X$info$arm)){

    return(calculate.arm.sparse(X, Y, windowSize, chromosome, method=method, params=params))
    
  } else {
    pArm <- calculate.arm.sparse(X, Y, windowSize, chromosome, 'p', method, params)
    qArm <- calculate.arm.sparse(X, Y, windowSize, chromosome, 'q', method, params)
    
 
    return(new("ChromosomeModels",
               models = append(pArm@models, qArm@models),
               chromosome = chromosome,
               method = method,
               params = params))
  }
}