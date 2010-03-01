setMethod(f="show",signature("ChromosomeArmModels"),
	function(object){
		cat("*** Dependency models for chromosome: ", as.character(object@chromosome), 
			as.character(object@arm), " ***\n", sep = "")
		cat("Number of models:", length(object@models),"\n")
		cat("Method used:", as.character(object@method), "; window size", getWindowSize(object), "\n")
		
		#Printing parameters
		if(length(object@params) > 0){
			cat("Method parameters: \n")
			names <- names(object@params)
			for (n in 1:length(names)){			
			        if(names[n] == 'H'){
					if(identical(object@params$H,diag(1,getWindowSize(object),getWindowSize(object)))){
						cat("- H: identity matrix","\n")
					}
					else {
					     cat("- ")
					     cat(matrix.print(object@params$H,"H"),"\n", sep="")				
					}
				}
				else {
					cat("- ",names[n], ": ", object@params[[names[n]]], "\n", sep="")
				}
			}
		}
		
		#Summary of score
		cat("Summary of dependency scores: \n")
		print(summary(getScore(object)))
		
		cat("******************************************\n")
	}
)

setMethod(f="show",signature("ChromosomeModels"),
	function(object){
		cat("*** Dependency models for chromosome: ", as.character(object@chromosome)," ***\n", sep = "")
		cat("Method used:", as.character(object@method), "; window size", getWindowSize(object), "\n")
		cat("Number of models in ", as.character(object@chromosome), "p: ",  getModelNumbers(object@pArmModels), 
			", ", as.character(object@chromosome), "q: ",  getModelNumbers(object@qArmModels), "\n", sep="")
	
		#Printing parameters
		if(length(object@params) > 0){
			cat("Method parameters: \n")
			names <- names(object@params)
			for (n in 1:length(names)){			
			        if(names[n] == 'H'){
					if(identical(object@params$H,diag(1,getWindowSize(object),getWindowSize(object)))){
						cat("- H: identity matrix","\n")
					}
					else {
					     cat("- ")
					     cat(matrix.print(object@params$H,"H"),"\n", sep="")				
					}
				}
				else {
					cat("- ",names[n], ": ", object@params[[names[n]]], "\n", sep="")
				}
			}
		}

		#Summary of score
		cat("Summary of dependency scores: \n")
		score <- c(getScore(getPArm(object)),getScore(getQArm(object)))
		print(summary(score))

		cat("******************************************\n")
	}
)

setMethod(f="show",signature("GenomeModels"),
	function(object){
		cat("*** Dependency models for genome ***\n", sep = "")
		cat("Method used:", as.character(object@method), "; window size", getWindowSize(object), "\n")
	
		#Printing parameters
		if(length(object@params) > 0){
			cat("Method parameters: \n")
			names <- names(object@params)
			for (n in 1:length(names)){			
			        if(names[n] == 'H'){
					if(identical(object@params$H,diag(1,getWindowSize(object),getWindowSize(object)))){
						cat("- H: identity matrix","\n")
					}
					else {
					     cat("- ")
					     cat(matrix.print(object@params$H,"H"),"\n", sep="")				
					}
				}
				else {
					cat("- ",names[n], ": ", object@params[[names[n]]], "\n", sep="")
				}
			}
		}

		#Summary of score
		cat("Summary of dependency scores: \n")
		score <- vector()
		for(i in 1:24){
			score <- c(score, getScore(getPArm(object[[i]])), getScore(getQArm(object[[i]])))
		}
		print(summary(score))

		cat("*************************************\n")
	}
)


setMethod(f="show",signature("DependencyModel"),
	function(object){
		cat("***", object@method, "dependency model for window size:",object@windowSize,"*** \n")
		#Gene name and location
		cat("Gene:",object@geneName)
		if(length(object@loc) > 0){
			cat("  Location: ")
			loc <- format((object@loc/1e6),digits=5)
			cat(loc,"Mbp",sep="")
		}
		cat("\n")

		#Score
		cat("Dependency score:", object@score,"\n")

		#WX print
		cat("- ")
		cat(matrix.print(object@W$X,"WX"),"\n", sep="")
		#string1 <- paste("WX: [1:", dim(object@W$X)[1], ", 1:",  dim(object@W$X)[2] ,"] ", sep = "") 
		#values <- format(object@W$X[1:min(4,length(object@W$X))],digits=3)
		#string2 <- ""
		#if(length(object@W$X) > 4) string2 <- "..."
		#cat(string1,values,string2,"\n")
		
		#WY print
		cat("- ")
		cat(matrix.print(object@W$Y,"WY"),"\n", sep="")
		#string1 <- paste("WY: [1:", dim(object@W$Y)[1], ", 1:",  dim(object@W$Y)[2] ,"] ", sep = "") 
		#values <- format(object@W$Y[1:min(4,length(object@W$Y))],digits=3)
		#string2 <- ""
		#if(length(object@W$Y) > 4) string2 <- "..."
		#cat(string1,values,string2,"\n")
		
		#Phi X
		cat("- ")
		cat(matrix.print(object@phi$X,"phiX"),"\n", sep="")
		#cat("phiX:",object@phi$X,"  phiY:",object@phi$Y,"\n")		
		#string1 <- paste("phiX: [1:", dim(object@phi$X)[1], ", 1:",  dim(object@phi$X)[2] ,"] ", sep = "") 
		#values <- format(object@phi$X[1:min(4,length(object@phi$X))],digits=3)
		#string2 <- ""
		#if(length(object@phi$X) > 4) string2 <- "..."
		#cat(string1,values,string2,"\n")
		
		#Phi Y
		cat("- ")
		cat(matrix.print(object@phi$Y,"phiY"),"\n", sep="")
		#cat("phiX:",object@phi$X,"  phiY:",object@phi$Y,"\n")		
		#string1 <- paste("phiY: [1:", dim(object@phi$Y)[1], ", 1:",  dim(object@phi$Y)[2] ,"] ", sep = "") 
		#values <- format(object@phi$Y[1:min(4,length(object@phi$Y))],digits=3)
		#string2 <- ""
		#if(length(object@phi$Y) > 4) string2 <- "..."
		#cat(string1,values,string2,"\n")
		
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

