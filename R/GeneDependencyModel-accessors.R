setAs("GeneDependencyModel", "DependencyModel", function(from){
    res <- new("DependencyModel", W = from@W, phi = from@phi, 
      score = from@score, method = from@method, params = from@params)	
    return(res)    
  }
)

setReplaceMethod(f="setLoc", signature("GeneDependencyModel"),
                 definition=(function(model,value) {
                   model@loc <- value
                   return(model)
                 }
))

setReplaceMethod(f="setGeneName", signature("GeneDependencyModel"),
	definition=(function(model,value) {
		model@geneName <- value
		return(model)
	}
))

setReplaceMethod("setChromosome","GeneDependencyModel", 
	function(model, value) { 
		model@chromosome <- value
		return(model)
	} 
) 

setReplaceMethod("setArm","GeneDependencyModel", 
	function(model, value) { 
		model@arm <- as.character(value)
		return(model)
	} 
) 

setMethod("getLoc","GeneDependencyModel", 
	function(model) { 
		return(model@loc) 
	} 
) 

setMethod("getGeneName","GeneDependencyModel", 
          function(model) { 
            return(model@geneName) 
          } 
) 

setMethod("getChromosome","GeneDependencyModel", 
	function(model) { 
		return(model@chromosome) 
	} 
) 

setMethod("getArm","GeneDependencyModel", 
	function(model) { 
		return(model@arm) 
	} 
) 

setMethod("getWindowSize","GeneDependencyModel", 
  function(model) { 
    return(nrow(model@W$X)) 
  } 
) 

setMethod("getZ","GeneDependencyModel",
  function(model, X, Y) {
    if (missing(X) || missing(Y)) stop("Original data sets are needed as parameters")
    #window <- fixed.window(X,Y,which(rownames(X$data) == getGeneName(model)),getWindowSize(model))   
    Xm <- t(centerData(t(X$data[rownames(model@W$X), ]), rm.na = TRUE))
    Ym <- t(centerData(t(Y$data[rownames(model@W$Y), ]), rm.na = TRUE))
    return(z.expectation(model, Xm, Ym)) 
  }
)

