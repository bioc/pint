pick.chr.arm <- function(X, chr, arm = NULL){

  if (chr == 'X') { chr <- 23 }
  if (chr == 'Y') { chr <- 24 }

  if (!is.null(arm)){
    # pick probes for this arm
    indices <- which(X$info$chr == as.numeric(chr) & X$info$arm == arm)
  } else {
    #pick probes for whole chromosome
    warning("arm information not given, modeling the whole chromosome")
    indices <- which(X$info$chr == as.numeric(chr))  
  }

  chrarmX <- list(data = X$data[indices,,drop=FALSE],
                  info = X$info[indices,])

  chrarmX
}
