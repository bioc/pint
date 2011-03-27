calculate.genome.sparse <- function(X, Y, windowSize, method = "pSimCCA", params = list()){
	chromosomeModelList = list()
	for (i in 1:24) {
		chromosomeModelList[i] = calculate.chr.sparse(X, Y, windowSize, i, method, params)
	}
	#chromosomeModelList[23] = calculate.chr.sparse(X, Y, windowSize, 'X', method, params)
	#chromosomeModelList[24] = calculate.chr.sparse(X, Y, windowSize, 'Y', method, params)
	return(new("GenomeModels", chromosomeModels = chromosomeModelList, method = method, params = params))
}
