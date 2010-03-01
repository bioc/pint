setMethod(f="[[", signature("ChromosomeArmModels"),
			definition=(function(x,i,j,drop) {
				return(x@models[[i]])
			} 
))

setReplaceMethod(f="[[",signature("ChromosomeArmModels"),
				definition=(function(x,i,j,value) {
					x@models[[i]] <- value
 					return(x)
 				}
))
	
setMethod("getScore","ChromosomeArmModels", 
	function(model) {
		scores <- vector()
		for (i in seq_along(model@models)) {
				scores[i] <- getScore(model[[i]])
		}
		return(scores) 
	} 
)			

setMethod("getLoc","ChromosomeArmModels", 
	function(model) {
		locs <- vector()
		for (i in seq_along(model@models)) {
			locs[i] <- getLoc(model[[i]])
		}
		return(locs) 
	} 
)			

setMethod("getGeneName","ChromosomeArmModels", 
	function(model) {
		geneNames = vector()
		for (i in seq_along(model@models)){
			geneNames[i] <- getGeneName(model[[i]])
		}
		return(geneNames) 
	} 
)	

setMethod("getChromosome","ChromosomeArmModels", 
	function(model) { 
		return(model@chromosome) 
	} 
) 

setMethod("getArm","ChromosomeArmModels", 
	function(model) { 
		return(model@arm) 
	} 
)

setMethod("getModelMethod","ChromosomeArmModels", 
	function(model) { 
		return(model@method) 	
	} 
)

setMethod("getParams","ChromosomeArmModels", 
	function(model) { 
		return(model@params) 
	} 
)

setMethod("getModelNumbers","ChromosomeArmModels", 
	function(model) { 
		return(length(model@models)) 
	} 
)

setMethod("getWindowSize","ChromosomeArmModels", 
	function(model) { 
		return(model@windowSize) 
	} 
)

setMethod("findHighestGenes","ChromosomeArmModels",
	function(model,num = 1) {
		scores <- getScore(model)
		genes <- getGeneName(model)
		data <- data.frame(scores,genes)
		# order dataframe and take num names of genes with highest scores 
		return(as.character(data[order(scores,decreasing=TRUE),]$genes[1:num]))	
	}
)

setMethod("findHighestModels","ChromosomeArmModels",
	function(model,num = 1) {
		scores <- getScore(model)
		indices <- 1:length(scores)
		data <- data.frame(scores,indices)
		# Order dataframe
		data <- data[order(scores,decreasing=TRUE),]
		returnList = list()
		for(i in 1:num){
			returnList = c(returnList,model[[data$indices[i]]])
		}
		return(returnList)
	}
)