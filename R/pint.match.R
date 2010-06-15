pint.match <- function(X, Y, max.dist = 1e7, chrs = NULL){

  X$info[["chr"]] <- as.character(X$info[["chr"]])
  Y$info[["chr"]] <- as.character(Y$info[["chr"]])
  X$info[X$info[["chr"]] == "X", "chr"] <- "23"
  X$info[X$info[["chr"]] == "Y", "chr"] <- "24"
  Y$info[Y$info[["chr"]] == "X", "chr"] <- "23"
  Y$info[Y$info[["chr"]] == "Y", "chr"] <- "24"

  # First order chromosomes 1...22, X, Y, then chromosomes with other names
  if (is.null(chrs)) {chrs <- c(as.character(1:24), sort(setdiff(unique(X$info[["chr"]]), as.character(1:24))))}

  message("Matching probes between the data sets..")
  
  xindices <- yindices <- vector()

  # If location information ('loc') for the probe is missing
  # but start and end positions are available, use them to calculate
  # probe middle bp location
  if (!"loc" %in% colnames(Y$info)) {
    # sometimes factors, sometimes numerics given; this should handle both cases correctly
    Y$info[["loc"]] <- (as.numeric(as.character(Y$info[, "start"])) + as.numeric(as.character(Y$info[, "end"])))/2
      #rowMeans(Y$info[, c("start","end")])
  }
  if (!"loc" %in% colnames(X$info)) {
    X$info[["loc"]] <- (as.numeric(as.character(X$info[, "start"])) + as.numeric(as.character(X$info[, "end"])))/2
    #X$info[["loc"]] <- rowMeans(X$info[, c("start","end")])
  }
  
  for (chr in chrs){
    # Note: chromosome arm information is used in the matching if it is available
    tmp <- get.neighboring.probes(X, Y, chr, max.dist)
    xindices <- c(xindices, tmp$xinds)
    yindices <- c(yindices, tmp$yinds)
  }

  # TODO: remove duplicates in the data matrix, unless segmented data is used
  # (which should be explicitly indicated: add the option to function call)
  
  newY <- list(data = as.matrix(Y$data[yindices,]), info = Y$info[yindices,])
  newX <- list(data = as.matrix(X$data[xindices,]), info = X$info[xindices,])

  return(list(X = newX, Y = newY))

}

closest <- function(a, vec){which.min(abs(a - vec))}

get.neighboring.probes <- function (X, Y, chr, max.dist, control.arms = TRUE) {

  # Use arm information if it is available and not blocked
  if (("arm" %in% names(X$info)) && ("arm" %in% names(Y$info)) && control.arms) {
    
    for (arm in c('p', 'q')){      
      # Investigate specified arm
      xchrinds <- which(as.character(X$info$chr) == chr & X$info$arm == arm)
      ychrinds <- which(as.character(Y$info$chr) == chr & Y$info$arm == arm)
      inds <- get.neighs(X, Y, xchrinds, ychrinds, max.dist) 
    }
  } else {
    # Investigate the whole chromosome
    xchrinds <- which(as.character(X$info$chr) == chr)
    ychrinds <- which(as.character(Y$info$chr) == chr)
    inds <- get.neighs(X, Y, xchrinds, ychrinds, max.dist)
  }

  # return the indices
  list(xinds = inds$xinds, yinds = inds$yinds)
  
}


get.neighs <- function (X, Y, xchrinds, ychrinds, max.dist) {

  xinds <- yinds <- NULL

  if ( length(xchrinds) > 0 && length(ychrinds) > 0 ){
    
    #Find indices of closest probe from Y for each from X
    xi <- 1:length(xchrinds)
    yi <- sapply(X$info$loc[xchrinds], closest, vec = Y$info$loc[ychrinds])

    # Remove duplicates
    keep <- !duplicated(yi)
    xi <- xi[keep]
    yi <- yi[keep]
      
    # Corresponding indices between X and Y
    xinds <- xchrinds[xi]
    yinds <- ychrinds[yi]

    # delete indices which are further from each other than threshold
    near <- (abs(X$info$loc[xinds] - Y$info$loc[yinds]) < max.dist)
    xinds <- xinds[near]
    yinds <- yinds[near]
  
    # calculate mean location for each pair for ordering of the pairs
    xy.loc <- (X$info$loc[xinds] + Y$info$loc[yinds])/2

    # ensure probes are ordered by location
    ord <- order(xy.loc)
    xinds <- xinds[ord]
    yinds <- yinds[ord]

  }

  list(xinds = xinds, yinds = yinds)
  
}
