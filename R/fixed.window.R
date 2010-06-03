fixed.window <-
function (X, Y, middleIndex, windowSize){
		
  # chromosome and arm of window
  chr <- X$info$chr[middleIndex]
  arm <- X$info$arm[middleIndex]

  # Location and name of the gene
  loc <- X$info$loc[middleIndex]
  geneName <- rownames(X$data)[middleIndex]
	
  # Indices for the window
  inds <- (middleIndex - (trunc((windowSize - 1)/2))) : (middleIndex + (trunc(windowSize/2)))	
	
  # Check that indices don't get out of bounds
  indsOutBounds <- (min(inds) < 0 || max(inds) > length(X$info$chr))

  # Check that chromosome and arm are the same
  sameArm <- (identical(X$info$chr[min(inds)], X$info$chr[max(inds)]) && 
              identical(X$info$arm[min(inds)], X$info$arm[max(inds)]))

  if(!indsOutBounds && sameArm){

    # TODO: check if rm.na is needed elsewhere in the code with centerData
    Xm <- t(matrix(centerData(t(X$data), rm.na = TRUE)[, inds], ncol = length(inds)))
    Ym <- t(matrix(centerData(t(Y$data), rm.na = TRUE)[, inds], ncol = length(inds)))
    
    res <- list(X = Xm, Y = Ym, loc = loc, geneName = geneName, fail = FALSE)
  }
  else {
    res <- list(fail = TRUE)
  }
  res
}

