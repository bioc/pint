setClass("DependencyModel", representation(W = "list", phi = "list", score = "numeric", loc = "numeric",  
		geneName = "character", windowSize = "numeric", 
		method = "character", params = "list"))

#setClass("SparseDependencyModel", contains = "DependencyModel")

setClass("ChromosomeArmModels", representation(models = "list", chromosome = "factor", arm = "factor", 
			windowSize = "numeric", method = "character", params = "list"))

setClass("ChromosomeModels", representation(pArmModels = "ChromosomeArmModels", qArmModels = "ChromosomeArmModels", chromosome = "factor",
			method = "character", params = "list"))

setClass("GenomeModels", representation(chromosomeModels = "list", method = "character", params = "list"))

#setClass("genomeArmModels", representation(ChromosomeArmModels = "list", arm, method = "character", parameters = "list"))

