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
          geneNames <- vector()
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

setMethod("isEmpty","ChromosomeArmModels",
	function(model) {
		return(length(model@models) == 0)
	}
)

setMethod("topGenes","ChromosomeArmModels",
	function(model, num = NA) {
          scores <- getScore(model)
          genes <- getGeneName(model)
          data <- data.frame(scores,genes)
          if (is.na(num)){
		    num <- length(scores)
		  }
          # order dataframe and take num names of genes with highest scores 
          return(as.character(data[order(scores,decreasing=TRUE),]$genes[1:num]))	
	}
          
)

setMethod("topModels","ChromosomeArmModels",
	function(model,num = 1) {
		scores <- getScore(model)
		indices <- 1:length(scores)
		data <- data.frame(scores,indices)
		# Order dataframe
		data <- data[order(scores,decreasing=TRUE),]
		returnList <- list()

		for(i in 1:num){
			returnList <- c(returnList,model[[data$indices[i]]])
		}
		return(returnList)
	}
)

setMethod("orderGenes","ChromosomeArmModels",
  function(model){

    scores <- getScore(model)
    genes <- getGeneName(model)
    data <- data.frame(scores,genes)
    data <- data[order(scores,decreasing=TRUE),]
    return(as.character(data$genes))
  }	
)

setMethod("findModel","ChromosomeArmModels",
  function(model, name){
    index = which(getGeneName(model) == name)
    if (length(index) > 0)
      return(model[[which(getGeneName(model) == name)]])
    stop("No model found")
  }
)
    
setMethod(f="as.data.frame",signature("ChromosomeArmModels"),
          definition=(function (x, row.names = NULL, optional = FALSE, ...) {

      model <- x
      genes <- c(getGeneName(model))
      scores <- c(getScore(model))
      locs <- c(getLoc(model))
      arms <- as.character(rep(getArm(model),length(scores)))
      chrs <- as.character(rep(getChromosome(model),length(genes)))
      #indices <- c(pindices,qindices)
      data <- data.frame(geneName = genes, dependencyScore = scores, chr = chrs, 
                         arm = arms, loc = locs, stringsAsFactors = FALSE)
      return(data)
    }
))