

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
            for(i in 1:length(model@chromosomeModels)) {
              if(getModelNumbers(model[[i]]) > 0)
              return(getWindowSize(model[[i]]))
            }
            return(NULL) 
          } 
) 

setMethod("topGenes",signature("GenomeModels"),
          function(model,num = NA) {
            scores <- vector()
            genes <- vector()
            for(i in 1:length(model@chromosomeModels)){
              scores <- c(scores, getScore(model[[i]]))
              genes <- c(genes, getGeneName(model[[i]]))
            }
            data = data.frame(scores,genes)
            if (is.na(num)){
              num = length(scores)
            }
            
            #order dataframe and take num names of genes with highest scores 
            if (num > getModelNumbers(model)){
              warning("Attempted to get more genes than there are models. Returning the genes from all models")
              num <- getModelNumbers(model)
            }
            return(as.character(data[order(scores,decreasing=TRUE),]$genes[1:num]))
          }
          )

setMethod("topModels","GenomeModels",
          function(model,num = 1) {
            
            scores <- vector()
            arm <- vector()
            chr <- vector()
            indices <- vector()
            for(i in 1:length(model@chromosomeModels)){
              chrScores <- getScore(model[[i]])
              indices <- c(indices,seq_along(chrScores))
              chr <- c(chr, rep(i,length(chrScores)))
              scores <- c(scores,chrScores)
              #armpscores <- getScore(getPArm(model[[i]]))
              #if (length(armpscores > 0)){
              #  arm <- c(arm,rep('p',length(armpscores)))
              #  indices <- c(indices,1:length(armpscores))
              #  chr <- c(chr, rep(i,length(armpscores)))
              #}		
              #armqscores <- getScore(getQArm(model[[i]]))
              #if (length(armqscores > 0)){
              #  arm <- c(arm,rep('q',length(armqscores)))
              #  indices <- c(indices,1:length(armqscores))
              #  chr <- c(chr, rep(i,length(armqscores)))
              #}
              #scores <- c(scores, armpscores, armqscores)
            }
            
            data <- data.frame(scores,indices,chr)
                                        #Order dataframe
            data <- data[order(scores,decreasing=TRUE),]
            returnList = list()

            if (num > getModelNumbers(model)){
              warning("Attempted to get more models than is calculated. Returning all models")
              num <- getModelNumbers(model)
            }

            for (i in 1:num) {
              #if (data$arm[i] == 'p')
        	  returnList = c(returnList,model[[data$chr[i]]][[data$indices[i]]])
              #else
              #returnList = c(returnList,getQArm(model[[data$chr[i]]])[[data$indices[i]]])               
            }
            return(returnList)
          }
          )

setMethod("getModelNumbers","GenomeModels",
  function(model){
    num <- 0
    for (i in 1:length(model@chromosomeModels)){
      num <- num + getModelNumbers(model[[i]])
    }
    return(num)
  }
)


setMethod("orderGenes","GenomeModels",

  function(model){

    scores <- vector()
    genes <- vector()
    for(i in 1:length(model@chromosomeModels)){
      scores <- c(scores, getScore(model[[i]]))
      genes <- c(genes, getGeneName(model[[i]]))
    }
    data <- data.frame(scores,genes,stringsAsFactors=FALSE)
    return(data[order(scores,decreasing=TRUE),])

  }	

)

setMethod("findModel","GenomeModels",
  function(model, name){

   for (i in 1:length(model@chromosomeModels)){
     index <- which(getGeneName(model[[i]]) == name)
     if(length(index) > 0) 
       #return(getPArm(model[[i]])[[pIndex[1]]])
       return(getPArm(model[[i]])[[1]]) # FIXME - was as above?
   }   
   stop("No model found")
  }
)

setMethod("as.data.frame","GenomeModels",
          function(x, ...){
    df <- data.frame(geneName = NULL, dependencyScore = NULL, chr = NULL,
                    arm = NULL, loc = NULL)
    for (i in 1:length(x@chromosomeModels)){
      df <- rbind(df, as.data.frame(x[[i]]))   
    }
    return(df)
  }
)
