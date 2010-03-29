pint.match <- function(X,Y){
  # Match X probes with closest from Y
  inds <- vector()
  chr <- 0
  for(i in 1:length(X$info$loc)){
    chrinds <- which(as.character(Y$info$chr) == as.character(X$info$chr[i]))
    if (length(chrinds > 0)){
      inds[i] <- min(chrinds) - 1 + which.min(abs(X$info$loc[i] - Y$info$loc[chrinds]))
    }
  }
  yinds <- unique(inds)
  xinds <- (1:length(inds))[!duplicated(inds)]

  newY <- list(data = Y$data[unique(inds),], info = Y$info[unique(inds),])
  newX <- list(data = Y$data[unique(inds),], info = Y$info[unique(inds),])
  list(X=newX,Y=newY)
  #list(inds=inds,xinds=xinds,yinds=yinds)
  #inds
}