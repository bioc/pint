setReplaceMethod(f="setLoc", signature("DependencyModel"),
	definition=(function(model,value) {
		model@loc <- value
		return(model)
	}
))

setReplaceMethod(f="setGeneName", signature("DependencyModel"),
	definition=(function(model,value) {
		model@geneName <- value
		return(model)
	}
))

setReplaceMethod("setChromosome","DependencyModel", 
	function(model, value) { 
		model@chromosome <- as.character(value) 
		return(model)
	} 
) 

setReplaceMethod("setArm","DependencyModel", 
	function(model, value) { 
		model@arm <- as.character(value)
		return(model)
	} 
) 

setMethod("getW","DependencyModel", 
	function(model) { 
		return(model@W) 
	} 
) 

setMethod("getPhi","DependencyModel", 
	function(model) { 
		return(model@phi) 
	} 
) 

setMethod("getScore","DependencyModel", 
	function(model) { 
		return(model@score) 
	} 
) 

setMethod("getLoc","DependencyModel", 
	function(model) { 
		return(model@loc) 
	} 
) 

setMethod("getGeneName","DependencyModel", 
	function(model) { 
		return(model@geneName) 
	} 
) 

setMethod("getParams","DependencyModel", 
	function(model) { 
		return(model@params) 
	} 
) 

setMethod("getModelMethod","DependencyModel", 
	function(model) { 
		return(model@method) 
	} 
) 

setMethod("getWindowSize","DependencyModel", 
	function(model) { 
		return(model@windowSize) 
	} 
) 

setMethod("getChromosome","DependencyModel", 
	function(model) { 
		return(model@chromosome) 
	} 
) 

setMethod("getArm","DependencyModel", 
	function(model) { 
		return(model@arm) 
	} 
) 
