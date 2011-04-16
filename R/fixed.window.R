fixed.window <-
function (X, Y, middleIndex, windowSize){
		
  # chromosome and arm of window
  chr <- X$info$chr[middleIndex]

  if ("arm" %in% colnames(X$info)) {
    arm <- X$info$arm[middleIndex]
  } else if ("arm" %in% colnames(Y$info)) {
    arm <- Y$info$arm[middleIndex]
  } else {arm = NULL}
  
  # Location
  loc <- X$info$loc[middleIndex]
  arm <- X$info$arm[middleIndex]  
	
  # Indices for the window. If indices exceed the limits, then use the
  # last window that fits (all windows must have same size, and it is
  # safer not to exclude any gene since it is in the ends of the arms)  
  #inds <- (middleIndex - (trunc((windowSize - 1)/2))) : (middleIndex + (trunc(windowSize/2)))	

  if (is.null(arm)) {
    arminds <- which(X$info$chr == chr)
  } else {
    arminds <- which(X$info$chr == chr & X$info$arm == arm)
  }

  nmin <- min(arminds)
  nmax <- max(arminds)
  if ((nmax - nmin + 1) < windowSize) {stop(paste("windowSize cannot exceed number of genes on chromosomal arm ", chr, arm, sep = ""))}
  wstart <- middleIndex - trunc((windowSize - 1)/2)
  wstop  <- middleIndex + trunc(windowSize/2)
  if (wstart < nmin)   { wstart <- nmin;  wstop <- wstart + windowSize - 1 }
  if (wstop > nmax) { wstop <- nmax; wstart <- wstop - windowSize + 1 }
  inds <- wstart:wstop
  
  # Check that indices don't get out of bounds
  indsOutBounds <- (min(inds) < nmin || max(inds) > nmax)

  # Check that chromosome and arm are the same in both ends
  if (!is.null(arm)) {
    sameArm <- (identical(X$info$chr[min(inds)], X$info$chr[max(inds)]) && 
                identical(X$info$arm[min(inds)], X$info$arm[max(inds)]))
  } else {
    sameArm <- (identical(X$info$chr[min(inds)], X$info$chr[max(inds)]))
  }
 
  if(!indsOutBounds && sameArm){
  #if( sameArm ) {

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

