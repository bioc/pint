pick.chr.arm <- function(X,chr,arm){

	if (chr == 'X') {
	   	if (!any(levels(X$info$chr) == 'X')) {
		   	chr = 23
		}
	}
	if (chr == 'Y') {
	   	if (!any(levels(X$info$chr) == 'Y')) {
		   	chr = 24
		}
	}

	chromindices = which(X$info$chr == chr)
	armindices = which(X$info$arm == arm)
	indices = intersect(chromindices,armindices)
	
	chrarmX = list(data = X$data[indices,,drop=FALSE], info = list(chr = X$info$chr[indices], 
					arm = X$info$arm[indices], loc = X$info$loc[indices])) 
	
	chrarmX
}