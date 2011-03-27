setMethod(f="show",signature("ChromosomeModels"),
  function(object){

    arms <- getArm(object)
    if (getModelNumbers(object) == 0){
      cat("No dependency models calculated for chromosome: ", as.character(object@chromosome), "\n", sep = "")
    } else {
      cat("*** Dependency models for chromosome:", as.character(object@chromosome))
      if (all(arms == 'p') || all(arms == 'q')){
        cat(as.data.frame(object)$arm[1])
      }          
      cat(" ***\n", sep = "")
      
      cat("Method used:", getModelMethod(object), "with window size", getWindowSize(object), "\n")
      if (any(arms == 'p') && any(arms == 'q')){
        cat("Number of models in ", as.character(object@chromosome), "p: ",  length(arms == 'p'), 
            ", ", as.character(object@chromosome), "q: ",  length(arms == 'q'), "\n", sep="")
      } else {
        cat("Number of models: ",  getModelNumbers(object), "\n", sep="")
      }
      #Printing parameters
      #if(length(object@params) > 0){
      #  cat("Method parameters: \n")
      #  names <- names(object@params)
      #  for (n in 1:length(names)){			
      #    if(names[n] == 'H'){
      #      if(identical(object@params$H,diag(1,getWindowSize(object),getWindowSize(object)))){
      #        cat("- H: identity matrix","\n")
      #      } else {
      #        cat("- ")
      #        cat(matrix.print(object@params$H,"H"),"\n", sep="")				
      #      }
      #    } else {
      #      cat("- ",names[n], ": ", object@params[[names[n]]], "\n", sep="")
      #    }
      #  }
      #}

      #Summary of score
      cat("Summary of dependency scores: \n")
      score <- getScore(object)
      print(summary(score))

      cat("******************************************\n")
    }
  }
)


setMethod(f="show",signature("GenomeModels"),
	function(object){
		if (getModelNumbers(object) == 0){
			cat("No dependency models calculated","\n")
		}
		cat("*** Dependency models for genome ***\n", sep = "")
	    cat("Number of models:", getModelNumbers(object),"\n")
		cat("Method used:", as.character(object@method), "; window size", getWindowSize(object), "\n")
		#Printing parameters
		#if(length(object@params) > 0){
		#	cat("Method parameters: \n")
		#	names <- names(object@params)
		#	for (n in 1:length(names)){			
		#	        if(names[n] == 'H'){
		#			if(identical(object@params$H,diag(1,getWindowSize(object),getWindowSize(object)))){
		#				cat("- H: identity matrix","\n")
		#			}
		#			else {
		#			     cat("- ")
		#			     cat(matrix.print(object@params$H,"H"),"\n", sep="")				
		#			}
		#		}
		#		else {
		#			cat("- ",names[n], ": ", object@params[[names[n]]], "\n", sep="")
		#		}
		#	}
		#}

		#Summary of score
		cat("Summary of dependency scores: \n")
		score <- vector()
		for(i in 1:24){
			score <- c(score, getScore(object[[i]]))
		}
		print(summary(score))

		cat("*************************************\n")
	}
)


setMethod(f="show",signature("GeneDependencyModel"),
  function(object){
    cat("***", object@method, "dependency model for window size:",getWindowSize(object),"*** \n")
    #Gene name and location
    if (length(object@geneName) > 0 | length(object@loc) > 0){   
      if (length(object@geneName) > 0) cat("Gene:",object@geneName)
      if(length(object@loc) > 0){
        cat("  Location: ")
        if (length(object@chromosome) > 0) cat(object@chromosome)
        if (length(object@arm) > 0) cat(object@arm)
        cat(", ")
        loc <- format((object@loc/1e6),digits=5)
        cat(loc,"Mbp",sep="")
      }
      cat("\n")
    }
    
    #Score
    cat("Dependency score:", object@score,"\n")

    if (is.null(object@W$X)){
      #W print
      cat("- ")
      cat(matrix.print(object@W$total,"W"),"\n", sep="")

      #Phi
      cat("- ")
      cat(matrix.print(object@phi$total,"phi"),"\n", sep="")
    }
    else {
      #WX print
      cat("- ")
      cat(matrix.print(object@W$X,"WX"),"\n", sep="")

      #WY print
      cat("- ")
      cat(matrix.print(object@W$Y,"WY"),"\n", sep="")
	
      #Phi X
      cat("- ")
      cat(matrix.print(object@phi$X,"phiX"),"\n", sep="")
    
      #Phi Y
      cat("- ")
      cat(matrix.print(object@phi$Y,"phiY"),"\n", sep="")
	}
    cat("************************************************\n")
  }
)


matrix.print <- function(matrix,name){
	# Prints size and 4 first values of matrix for show-methods

	if(any(is.na(matrix))){
		return(cat(name,": NA",rep=""))
	}
	else { 
		string1 <- paste(name, ": [1:", dim(matrix)[1], ", 1:",  dim(matrix)[2] ,"] ", sep = "")
		values <- format(matrix[1:min(4,length(matrix))],digits=3)
		string2 <- ""
		if(length(matrix) > 4) string2 <- "..."
		return(cat(string1,values,string2))
	}
}

