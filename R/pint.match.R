pint.match <- function(X, Y, max.dist = 1e7){

  message("Matching probes between the data sets..")
  
  chrs <- c(1:22, 'X', 'Y')
  xindices <- yindices <- vector()

  for (i in chrs){
    for (k in c('p', 'q')){
      ychrinds <- which(as.character(Y$info$chr) == i & Y$info$arm == k)
      xchrinds <- which(as.character(X$info$chr) == i & X$info$arm == k)
      if(length(xchrinds) > 0 && length(ychrinds) > 0){

        # Find indices of closest probe from Y for each from X
        inds <- sapply(X$info$loc[xchrinds], closest, vec = Y$info$loc[ychrinds])

        # Indices to X and Y with duplicates deleted
        yinds <- min(ychrinds-1) + unique(inds)

        if (any(duplicated(inds))) {message("Probes with identical location observed. Removing duplicates.")}

        # Remove duplicates
        xinds <- min(xchrinds-1) + (1:length(inds))[!duplicated(inds)]

        # delete indices which are further from each other than threshold
        near <- (abs(X$info$loc[xinds] - Y$info$loc[yinds]) < max.dist)      
        xindices <- c(xindices, xinds[near])
        yindices <- c(yindices, yinds[near])
      }
    }
  }

  newY <- list(data = Y$data[yindices,], info = Y$info[yindices,])
  newX <- list(data = X$data[xindices,], info = X$info[xindices,])

  return(list(X = newX, Y = newY))

}

closest <- function(a, vec){
  which.min(abs(a-vec))
}
