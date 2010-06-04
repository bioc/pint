fixed.window <-
function (X, Y, middleIndex, windowSize){
		
  # chromosome and arm of window
  chr <- X$info$chr[middleIndex]
  arm <- X$info$arm[middleIndex]

  # Location
  loc <- X$info$loc[middleIndex]
	
  # Indices for the window
  inds <- (middleIndex - (trunc((windowSize - 1)/2))) : (middleIndex + (trunc(windowSize/2)))	
	
  # Check that indices don't get out of bounds
  indsOutBounds <- (min(inds) < 0 || max(inds) > length(X$info$chr))

  # Check that chromosome and arm are the same
  sameArm <- (identical(X$info$chr[min(inds)], X$info$chr[max(inds)]) && 
              identical(X$info$arm[min(inds)], X$info$arm[max(inds)]))

  if(!indsOutBounds && sameArm){

    Xm <- t(centerData(t(X$data), rm.na = TRUE)[, inds])
    Ym <- t(centerData(t(Y$data), rm.na = TRUE)[, inds])

    #name of the gene
    geneName <- rownames(X$data)[[middleIndex]]

    res <- list(X = Xm, Y = Ym, loc = loc, geneName = geneName, fail = FALSE)
  }
  else {
    res <- list(fail = TRUE)
  }
  res
}

