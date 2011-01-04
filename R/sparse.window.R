sparse.window <- function (X, Y, xIndex, windowSize) {
		
  # chromosome and arm of window
  chromosome <- as.character(X$info$chr[xIndex])
  arm <- as.character(X$info$arm[xIndex])
	
  if (chromosome == 'X') {
    if (!any(levels(X$info$chr) == 'X')) {
      chr <- 23
    }
  }
  if (chromosome == 'Y') {
    if (!any(levels(X$info$chr) == 'Y')) {
      chr <- 24
    }
  }
  
  # Indices of genes in window
  inds <- vector()
  loc <- X$info$loc[xIndex]
  # Indices for samples in same chromosome arm
  chrIndices <- which(Y$info$chr == chromosome & Y$info$arm == arm)

  # Test if Y has no samples with this chromosome arm
  if (length(chrIndices) == 0) {
    res <- list(fail = TRUE)
    return(res)
  } 

  # Index in that arm that is nearest to loc
  middleIndex <- which.min(abs(Y$info$loc[chrIndices] - loc)) + min(chrIndices) - 1


  inds[1] <- middleIndex
  geneName <- dimnames(X$data)[[1]][xIndex]
  
  # Indices for iterative finding
  lowerIndex <- middleIndex - 1
  higherIndex <- middleIndex + 1
	
  # Distances to lower and hgher genes in iteration
  lowerDist <- 1e12
  higherDist <- 1e12
  lowerEdge <- FALSE
  higherEdge <- FALSE
  
  armBig <- FALSE

  if(windowSize > 1) {
    # Pick windowSize-1 nearest genes to window
    for (i in 2:windowSize) {
		
      # Claculate distance to lower gene
      # Check if lower nearest gene is in the same chromosomal arm
      if(lowerIndex > 0) {
        if(Y$info$chr[lowerIndex] == chromosome && Y$info$arm[lowerIndex] == arm){
          lowerDist <- loc - Y$info$loc[lowerIndex]
        }
        else {
          lowerDist <- 1e12
          lowerEdge <- TRUE
        }
      }
      else {
        lowerDist <- 1e12
        lowerEdge <- TRUE
      }
      
         # Calculate distance to higher gene
         # Check if higher nearest gene is in the same chromosomal arm
      if (higherIndex <= length(Y$info$loc)) {
        if(Y$info$chr[higherIndex] == chromosome && Y$info$arm[higherIndex] == arm){
          higherDist <- Y$info$loc[higherIndex] - loc
        }
        else {
          higherDist <- 1e12
          higherEdge <- TRUE
        }
      }
      else {
        higherDist <- 1e12
        higherEdge <- TRUE
      }

      #test if lower gene is closer to middleIndex
      if (!lowerEdge && lowerDist < higherDist) {
        inds[i] <- lowerIndex
        lowerIndex <- lowerIndex - 1
      }
      #test if higher gene is closer to middleIndex
      else if(!higherEdge && higherDist < lowerDist) {
        inds[i] <- higherIndex
        higherIndex <- higherIndex + 1
      }
      #test if higher and lower genes are not in the same chromosomal arm
      else {
        armBig <- TRUE
        break;
      }
    }
  }
  
  inds <- sort(inds, decreasing = FALSE)
  
  if(!armBig){
    Xm <- X$data[xIndex,,drop=FALSE]
    Ym <- Y$data[inds,,drop=FALSE]
    
    res <- list(X = Xm,Y = Ym,loc = loc, geneName = geneName, edge = (higherEdge || lowerEdge), fail = FALSE)
  }
  else {
    res <- list(fail = TRUE)
    warning("Chromosome contains less than windowSize genes.")
  }
  res
}

