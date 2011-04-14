setMethod(f="[[", signature("ChromosomeModels"),
  definition=(function(x,i,j,drop) {
    if (i == 'p') return(getPArm(x))
    if (i == 'q') return(getQArm(x)) 
    return(x@models[[i]])
  } 
))

setReplaceMethod(f="[[",signature("ChromosomeModels"),
  definition=(function(x,i,j,value) {
    x@models[[i]] <- value
    return(x)
  }
))

setMethod("getChromosome","ChromosomeModels", 
  function(model) { 
    return(model@chromosome) 
  } 
) 

setMethod("getPArm","ChromosomeModels", 
  function(model) {
    arms <- getArm(model)
	  if (all(arms == "") && length(arms) > 0){
	    stop("cannot return dependency models for an arm because no arm information was given")
	  }
	
	  return(new("ChromosomeModels",                                  
               models = model@models[arms == 'p'],                                   
               chromosome = getChromosome(model), 
               method = getModelMethod(model),                                                
               params = getParams(model)))       
  } 
) 

setMethod("getQArm","ChromosomeModels", 
  function(model) {
    arms <- getArm(model)
	if (all(arms == "") && length(arms) > 0){
	  stop("cannot return dependency models for an arm because no arm information was given")
	}
	
	return(new("ChromosomeModels",                                  
               models = model@models[arms == 'q'],                                   
               chromosome = getChromosome(model), 
               method = getModelMethod(model),                                                
               params = getParams(model)))       
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
    return(getWindowSize(model@models[[1]])) 
  } 
) 
setMethod("isEmpty","ChromosomeModels",
  function(model) {
    return(length(model@models) == 0)
  }
)

setMethod("getModelNumbers","ChromosomeModels",
  function(model) {
    return(length(model@models))
  }
)

setMethod("getScore","ChromosomeModels", 
  function(model) {
    scores <- vector()
    for (i in seq_along(model@models)) {
      scores[i] <- getScore(model[[i]])
    }
    return(scores) 
})	

setMethod("topGenes", "ChromosomeModels",
  function(model, num = NA) {

    scores <- getScore(model)
    genes <- getGeneName(model)        
    data <- data.frame(scores, genes)

    if (is.na(num)){
      num <- length(scores)
    }
    #order dataframe and take num names of genes with highest scores 
    return(as.character(data[order(scores, decreasing = TRUE),]$genes[1:num]))
})

setMethod("getGeneName", "ChromosomeModels",
  function(model){
    genes <- vector()
      for (i in seq_along(model@models)) {
        genes[i] <- getGeneName(model[[i]])
      }
    return(genes) 
	} 
)	

setMethod("getLoc","ChromosomeModels", 
	function(model) {
		locs <- vector()
		for (i in seq_along(model@models)) {
			locs[i] <- getLoc(model[[i]])
		}
		return(locs) 
	} 
)			

setMethod("getArm","ChromosomeModels", 
	function(model) {
    return(sapply(model@models,getArm))
	} 
)
	
setMethod("topModels","ChromosomeModels",
  function(model, num = 1) {
    scores <- getScore(model)
    genes <- getGeneName(model)        
    data <- data.frame(scores, genes)
    indices <- seq_along(scores)
    data <- data.frame(scores,indices)

    #Order dataframe
    data <- data[order(scores, decreasing = TRUE), ]
    returnList <- list()
    if (num > 1){
      for (i in 1:num) {
          returnList <- c(returnList,model[[data$indices[i]]])
      }
      return(returnList)
    } else {
      return(model[[data$indices[1]]])
    }
  }
)

setMethod("orderGenes","ChromosomeModels",
  function(model){

    return(topGenes(model))
  }	
)

setMethod("findModel","ChromosomeModels",
  function(model, name){
    index = which(getGeneName(model) == name)
    if (length(index) > 0)
      return(model[[index[1]]])
    stop("No model found")
  }
)

setMethod(f="as.data.frame",signature("ChromosomeModels"),
  definition=(function (x, row.names = NULL, optional = FALSE, ...) {

    model <- x
    genes <- getGeneName(model)
    scores <- getScore(model)
    locs <- getLoc(model)
    arms <- getArm(model)
    chrs <- rep(getChromosome(model),getModelNumbers(model))
    if (!any(arms == "")){  
      data <- data.frame(geneName = genes, dependencyScore = scores, chr = chrs, 
                         arm = arms, loc = locs, stringsAsFactors = FALSE)
    } else {
      data <- data.frame(geneName = genes, dependencyScore = scores, chr = chrs, 
                         loc = locs, stringsAsFactors = FALSE)
    } 
    return(data)
  }
))
