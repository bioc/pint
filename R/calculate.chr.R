calculate.chr <- function(X, Y, windowSize, chromosome, method = "pSimCCA", params = list(), segmented, priors, regularized){

  pArm <- calculate.arm(X, Y, windowSize, chromosome, 'p', method, params, segmented, priors, regularized)
  qArm <- calculate.arm(X, Y, windowSize, chromosome, 'q', method, params, segmented, priors, regularized)
  chromosome <- factor(chromosome, levels = levels(X$info$chr))

  return(new("ChromosomeModels",
             pArmModels = pArm,
             qArmModels = qArm,
             chromosome = chromosome,
             method = method,
             params = params))

}
