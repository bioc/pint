pick.chr.arm <- function(X,chr,arm){

	chromindices = which(X$info$chr == chr)
	armindices = which(X$info$arm == arm)
	indices = intersect(chromindices,armindices)
	
	chrarmX = list(data = X$data[indices,], info = list(chr = X$info$chr[indices], 
					arm = X$info$arm[indices], loc = X$info$loc[indices])) 
	
	chrarmX
}