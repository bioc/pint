
setMethod("getModelMethod","GenomeModels", 
	function(model) {
		return(model@method) 
	} 
)	

setMethod("getParams","GenomeModels", 
	function(model) {
		return(model@params) 
	} 
)	

setMethod(f="[[", signature("GenomeModels"),
			definition=(function(x,i,j,drop) {
				if(i == 'X') 
					i <- 23
				if(i == 'Y') 
					i <- 24
				return(x@chromosomeModels[[i]])
			} 
))

setReplaceMethod(f="[[",signature("GenomeModels"),
				definition=(function(x,i,j,value) {
					if(i == 'X') 
						i <- 23
					if(i == 'Y') 
						i <- 24
					
					x@chromosomeModels[[i]] <- value
					return(x) 
 				}
))

setMethod("getWindowSize","GenomeModels", 
	function(model) { 
		return(getWindowSize(getPArm(model[[1]]))) 
	} 
) 

setMethod("findHighestGenes",signature("GenomeModels"),
	function(model,num = 1) {
		scores <- vector()
		genes <- vector()
		for(i in 1:24){
			scores <- c(scores, getScore(getQArm(model[[i]])))
			scores <- c(scores, getScore(getPArm(model[[i]])))
			genes <- c(genes, getGeneName(getQArm(model[[i]])))
			genes <- c(genes, getGeneName(getPArm(model[[i]])))
		}
		data = data.frame(scores,genes)
		#order dataframe and take num names of genes with highest scores 
		return(as.character(data[order(scores,decreasing=TRUE),]$genes[1:num]))
	}
)