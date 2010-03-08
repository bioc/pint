calculate.genome <- function(X, Y, windowSize, method = "pSimCCA", params = list()){
	chromosomeModelList = list()
	for (i in 1:22) {
		chromosomeModelList[i] = calculate.chr(X, Y, windowSize, i, method, params)
	}
	chromosomeModelList[23] = calculate.chr(X, Y, windowSize, 'X', method, params)
	chromosomeModelList[24] = calculate.chr(X, Y, windowSize, 'Y', method, params)
	return(new("GenomeModels", chromosomeModels = chromosomeModelList, method = method, params = params))
}
