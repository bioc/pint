
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

setMethod("findHighestModels","GenomeModels",
        function(model,num = 1) {

		scores <- vector()
		arm <- vector()
		chr <- vector()
		indices <- vector()
		for(i in 1:24){
		        armpscores <- getScore(getPArm(model[[i]]))
			if (length(armpscores > 0)){
				arm <- c(arm,rep('p',length(armpscores)))
				indices <- c(indices,1:length(armpscores))
				chr <- c(chr, rep(i,length(armpscores)))
			}		
		        armqscores <- getScore(getQArm(model[[i]]))
			if (length(armqscores > 0)){
				arm <- c(arm,rep('q',length(armqscores)))
				indices <- c(indices,1:length(armqscores))
				chr <- c(chr, rep(i,length(armqscores)))
			}
			scores <- c(scores, armpscores, armqscores)
		}

                data <- data.frame(scores,indices,arm,chr)
                #Order dataframe
                data <- data[order(scores,decreasing=TRUE),]
                returnList = list()
                for (i in 1:num) {
                        if (data$arm[i] == 'p')
                                returnList = c(returnList,getPArm(model[[data$chr[i]]])[[data$indices[i]]])
                        else
                                returnList = c(returnList,getQArm(model[[data$chr[i]]])[[data$indices[i]]])               
		}
                return(returnList)
        }
)
