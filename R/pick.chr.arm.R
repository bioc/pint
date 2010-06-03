pick.chr.arm <- function(X, chr, arm){

  if (chr == 'X') {
    if (!any(levels(X$info$chr) == 'X')) {
      chr <- 23
    }
  }
  if (chr == 'Y') {
    if (!any(levels(X$info$chr) == 'Y')) {
      chr <- 24
    }
  }

  # pick probes for this arm
  indices <- which(X$info$chr == chr & X$info$arm == arm)
  
  chrarmX <- list(data = X$data[indices,,drop=FALSE],
                  info = X$info[indices,])
  #chrarmX <- list(data = X$data[indices, ],
  #                info = X$info[indices, ])

  
  chrarmX
}
