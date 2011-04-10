

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

#setReplaceMethod(f="setSegmentedData", signature("GenomeModels"),
#                 definition=(function(model,value) {
#      if (!is.logical(value)){
#        stop("Incorrect value given. Use TRUE or FALSE")
#      }
#      if (is.null(model@params$segmentedData)){
#        model@params <- c(model@params, segmentedData = value)
#      }
#      else {
#        model@params$segmentedData <- value 
#      }
#      return(model)
#    }
#))

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
            for(i in 1:24) {
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
            for(i in 1:24){
              scores <- c(scores, getScore(getQArm(model[[i]])))
              scores <- c(scores, getScore(getPArm(model[[i]])))
              genes <- c(genes, getGeneName(getQArm(model[[i]])))
              genes <- c(genes, getGeneName(getPArm(model[[i]])))
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

            if (num > getModelNumbers(model)){
              warning("Attempted to get more models than is calculated. Returning all models")
              num <- getModelNumbers(model)
            }

            for (i in 1:num) {
              if (data$arm[i] == 'p')
                returnList = c(returnList,getPArm(model[[data$chr[i]]])[[data$indices[i]]])
              else
                returnList = c(returnList,getQArm(model[[data$chr[i]]])[[data$indices[i]]])               
            }
            return(returnList)
          }
          )

setMethod("getModelNumbers","GenomeModels",
  function(model){
    num <- 0
    for (i in 1:24){
      num <- num + getModelNumbers(model[[i]])
    }
    return(num)
  }
)


setMethod("orderGenes","GenomeModels",
  function(model){

    scores <- vector()
    genes <- vector()
    for(i in 1:24){
      scores <- c(scores, getScore(getQArm(model[[i]])))
      scores <- c(scores, getScore(getPArm(model[[i]])))
      genes <- c(genes, getGeneName(getQArm(model[[i]])))
      genes <- c(genes, getGeneName(getPArm(model[[i]])))
    }
    data <- data.frame(scores,genes,stringsAsFactors=FALSE)
    return(data[order(scores,decreasing=TRUE),])
  }	
)

setMethod("findModel","GenomeModels",
  function(model, name){

   for (i in 1:24){
     pIndex <- which(getGeneName(getPArm(model[[i]])) == name)
     if(length(pIndex) > 0) 
       return(getPArm(model[[i]])[[pIndex[1]]])
     qIndex <- which(getGeneName(getQArm(model[[i]])) == name)
     if(length(qIndex) > 0) 
       return(getQArm(model[[i]])[[qIndex[1]]])
   }   
   stop("No model found")
  }
)

setMethod("as.data.frame","GenomeModels",
          function(x, ...){
    df <- data.frame(geneName = NULL, dependencyScore = NULL, chr = NULL,
                    arm = NULL, loc = NULL)
    for (i in 1:24){
      df <- rbind(df, as.data.frame(x[[i]]))   
    }
    return(df)
  }
)