pint.data <- function(data, info, impute = TRUE, replace.inf = TRUE){

  # Probe data
  if (class(data) == "data.frame"){
    data <- as.matrix(data)
  }


  # Impute missing values
  if (impute) {
    message("Imputing missing values..")
    data <- impute(data)
  }
  
  # Replace infinite values by highest possible seen in the data
  if (replace.inf) {
     message("Replacing infinite values with highest non-infinite values seen in the data")
     inds <- is.infinite(data)
     data[inds] <- sign(data[inds])*max(abs(data[!inds])) # note the sign
     message(paste("...", 100*mean(inds), "percent of the values replaced."))
  }
  
  
  ## Dealing with NaNs
  ## variables (rows) with NaNs
  #rows <- which(rowSums(is.nan(data)) > 0)
  #if (length(rows) > 0) {
  #  nans <- 0;
  #  # Replace NaNs with samples from normal distribution
  ##  for(i in rows) {
  #    mean <- mean(data[i,],na.rm=TRUE)
  ##    sd <- sd(data[i,],na.rm=TRUE)
  ##    inds <- which(is.nan(data[i,]))
  #    for(j in inds){
  #        data[i,j] <- rnorm(1,mean,sd)
  #        nans <- nans + 1
  #    }
  #  }
  #  warning(paste(nans,"NaNs in probe data imputed with samples from normal distribution"))
  #}
  
  
 
  # Location
  if (is.null(info$loc) && is.null(info$bp)) {
    loc <- (info$start + info$end)/2
  } 
  else {
    if (!is.null(info$bp))
      loc <- info$bp
    else
      loc <- info$loc
  }
  
  # Arm
  if (is.null(info$arm)) {
    arm <- factor(rep("p",length(loc)),levels=c("p","q"))
  }
  else {
    arm <- info$arm
  }
  
  info <- data.frame(chr = factor(info$chr, levels = c(1:22,"X","Y")), arm = arm, loc = as.numeric(loc))

  # Order data by chr, arm and loc
  ord <- order(info$chr,info$arm,info$loc)
  data <- data[ord,]
  info <- info[ord,]

  list(data = data, info = info)
}
