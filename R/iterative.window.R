iterative.window <- function (X, Y, middleIndex, windowSize) {
		
	# Flag to show if window is in the edge of chromosomal arm
	higherEdge <- FALSE
	lowerEdge <- FALSE
	
	# Flag to show if chromosomal arm is shorter than window size
	armBig <- FALSE
	
	# chromosome and arm of window
	chromosome <- X$info$chr[middleIndex]
	arm <- X$info$arm[middleIndex]

	# Indices of genes in window
	inds <- vector()
	inds[1] = middleIndex
	loc <- X$info$loc[middleIndex]
	geneName <- dimnames(X$data)[[1]][middleIndex]

	# Indices for iterative finding
	lowerIndex <- middleIndex-1
	higherIndex <- middleIndex+1
	
	# Distances to lower and hgher genes in iteration
	lowerDist <- 1e12
	higherDist <- 1e12

	if(windowSize > 1) {
		# Pick windowSize-1 nearest genes to window
		for (i in 2:windowSize) {
		
			# Claculate distance to lower gene
			# Check if lower nearest gene is in the same chromosomal arm
			if(lowerIndex > 0) {
				if(X$info$chr[lowerIndex] == chromosome && X$info$arm[lowerIndex] == arm){
					lowerDist <- loc - X$info$loc[lowerIndex]
				}
				else {
					lowerEdge <- TRUE
					lowerDist <- 1e12
				}
			}
			else {
				lowerEdge <- TRUE
				lowerDist <- 1e12
			}
		
			# Claculate distance to higher gene
			# Check if higher nearest gene is in the same chromosomal arm
			if (higherIndex <= length(X$info$loc)) {
				if(X$info$chr[higherIndex] == chromosome && X$info$arm[higherIndex] == arm){
					higherDist <- X$info$loc[higherIndex] - loc
				}
				else {
					higherEdge <- TRUE
					higherDist <- 1e12
				}
			}
			else {
				higherEdge <- TRUE
				higherDist <- 1e12
			}
			#print(paste("dists:",lowerDist,higherDist,"edges:",lowerEdge,higherEdge))
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

	inds = sort(inds, decreasing = FALSE)

	if(!armBig){
		Xm = X$data[inds,]
		Ym = Y$data[inds,]
	
		res = list(X = Xm,Y = Ym,loc = loc, geneName = geneName, edge = (higherEdge || lowerEdge), fail = FALSE)
	}
	else {
		res = list(fail = TRUE)
		warning("Chromosome contains less than windowSize genes.")
	}
	res
}

