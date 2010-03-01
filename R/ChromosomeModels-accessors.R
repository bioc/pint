setMethod("getChromosome","ChromosomeModels", 
	function(model) { 
		return(model@chromosome) 
	} 
) 

setMethod("getPArm","ChromosomeModels", 
	function(model) { 
		return(model@pArmModels) 
	} 
) 

setMethod("getQArm","ChromosomeModels", 
	function(model) { 
		return(model@qArmModels) 
	} 
) 

setMethod("getParams","ChromosomeModels", 
	function(model) { 
		return(model@params) 
	} 
) 

setMethod("getModelMethod","ChromosomeModels", 
	function(model) { 
		return(model@method) 
	} 
) 

setMethod("getWindowSize","ChromosomeModels", 
	function(model) { 
		return(getWindowSize(getPArm(model))) 
	} 
) 

setMethod("findHighestGenes","ChromosomeModels",
	function(model,num = 1) {
		pscores <- getScore(getPArm(model))
		qscores <- getScore(getQArm(model))
		scores <- c(pscores,qscores)
		
		pgenes <- getGeneName(getPArm(model))
		qgenes <- getGeneName(getQArm(model))
		genes <- c(pgenes,qgenes)

		data <- data.frame(scores,genes)
		#order dataframe and take num names of genes with highest scores 
		return(as.character(data[order(scores,decreasing=TRUE),]$genes[1:num]))
	}
)

setMethod("findHighestModels","ChromosomeModels",
	function(model,num = 1) {
		pscores <- getScore(getPArm(model))
		qscores <- getScore(getQArm(model))
		scores <- c(pscores,qscores)
		pindices <- 1:length(pscores)
		qindices <- 1:length(qscores)
		indices <- c(pindices,qindices)
		data <- data.frame(scores,indices,arm = c(rep('p',length(pscores)), rep('q',length(qscores))))
		#Order dataframe
		data <- data[order(scores,decreasing=TRUE),]
		returnList = list()
		for (i in 1:num) {
			if (data$arm[i] == 'p')
				returnList = c(returnList,getPArm(model)[[data$indices[i]]])
			else
				returnList = c(returnList,getQArm(model)[[data$indices[i]]])
		}
		return(returnList)
	}
)