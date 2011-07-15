calculate.genome <- function(X, Y, windowSize, method = "pSimCCA", params = list(), match.probes, priors, regularized){

  chromosomeModelList <- list()
  for (i in 1:24) {
    chromosomeModelList[i] <- calculate.chr(X, Y, windowSize, i, method, params, match.probes, priors, regularized)
  }

  return(new("GenomeModels", chromosomeModels = chromosomeModelList, method = method, params = params))

}
