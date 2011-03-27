calculate.genome <- function(X, Y, windowSize, method = "pSimCCA", params = list(), segmented, priors, regularized){

  chromosomeModelList <- list()
  for (i in 1:24) {
    chromosomeModelList[i] <- calculate.chr(X, Y, windowSize, i, method, params, segmented, priors, regularized)
  }

  #chromosomeModelList[23] <- calculate.chr(X, Y, windowSize, 'X', method, params, segmented, priors, regularized)
  #chromosomeModelList[24] <- calculate.chr(X, Y, windowSize, 'Y', method, params, segmented, priors, regularized)

  return(new("GenomeModels", chromosomeModels = chromosomeModelList, method = method, params = params))

}
