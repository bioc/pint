pint.match <- function(X, Y, max.dist = 1e7, chrs = NULL,
useSegmentedData = FALSE, impute = TRUE, replace.inf = TRUE){

  if (all(is.na(X$data))) {stop("X data is empty/NA.")}
  if (all(is.na(Y$data))) {stop("Y data is empty/NA.")}
  
  # Match genes/probes/clones in X/Y data sets based on location information
  #X <- geneExp; Y <- geneCopyNum

  # If same number of rows and columns for data and info fields, 
  # assume that they match between data and info fields and give both
  # same rownames
  if (nrow(X$data) == nrow(X$info) && !rownames(X$data) == rownames(X$info)) {
     warning("The data and info fields in the X data set had different
rownames; using the rownames from data field for both.")
     rownames(X$info) <- rownames(X$data)
  }
  if (nrow(Y$data) == nrow(Y$info)) {
     warning("The data and info fields in the Y data set had different
rownames; using the rownames from data field for both.")  
    rownames(Y$info) <- rownames(Y$data)
  }

  # First, provide information in a form that has corresponding rows in
  # data and info matrices (only take those available in both)
  coms <- intersect(rownames(X$data), rownames(X$info))
  coms <- setdiff(coms, c(""))
  X$data <- X$data[coms,]
  X$info <- X$info[coms,]

  coms <- intersect(rownames(Y$data), rownames(Y$info))
  coms <- setdiff(coms, c(""))
  Y$data <- Y$data[coms,]
  Y$info <- Y$info[coms,]

  X$info[["chr"]] <- as.character(X$info[["chr"]])
  Y$info[["chr"]] <- as.character(Y$info[["chr"]])

  # X/Y chromosome for data set X
  if ("X" %in% X$info[["chr"]]) {X$info[X$info[["chr"]] == "X", "chr"] <- "23"}
  if ("Y" %in% X$info[["chr"]]) {X$info[X$info[["chr"]] == "Y", "chr"] <- "24"}

  # X/Y chromosome for data set X
  if ("X" %in% Y$info[["chr"]]) {Y$info[Y$info[["chr"]] == "X", "chr"] <- "23"}
  if ("Y" %in% Y$info[["chr"]]) {Y$info[Y$info[["chr"]] == "Y", "chr"] <- "24"}

  # Quarantee that there are no duplicated rows (probes) in the data
  if (!useSegmentedData){
    dupl <- duplicated(X$data)
    if (any(dupl)) {
      cat("Removing duplicate probe signals on X data..\n")
      X$data <- X$data[!dupl, ]
      X$info <- X$info[!dupl, ]
    }

    dupl <- duplicated(Y$data)
    if (any(dupl)) {
      cat("Removing duplicate probe signals on Y data..\n")
      Y$data <- Y$data[!dupl, ]
      Y$info <- Y$info[!dupl, ]
    }
  }


  # Impute missing values, replace infinite values etc.
  X <- pint.data(X$data, X$info, impute, replace.inf)
  Y <- pint.data(X$data, X$info, impute, replace.inf)

  
  # First order chromosomes 1...22, X, Y, then chromosomes with other names
  if (is.null(chrs)) {
    chrs <- c(as.character(1:24), sort(setdiff(unique(X$info[["chr"]]), as.character(1:24))))
  } else {
    chrs <- as.character(chrs)
  }

  # If location information ('loc') for the probe is missing
  # but start and end positions are available, use them to calculate
  # probe middle bp location
  if (!"loc" %in% colnames(Y$info)) {
    # sometimes factors, sometimes numerics given; this should handle both cases correctly
    Y$info[["loc"]] <- (as.numeric(as.character(Y$info[, "start"])) + as.numeric(as.character(Y$info[, "end"])))/2
  }
  if (!"loc" %in% colnames(X$info)) {
    X$info[["loc"]] <- (as.numeric(as.character(X$info[, "start"])) + as.numeric(as.character(X$info[, "end"])))/2
  }
  

  message("Matching probes between the data sets..")
  xindices <- yindices <- vector()
  for (chr in chrs){
    # Note: chromosome arm information is used in the matching if it is available
    tmp <- get.neighboring.probes(X, Y, chr, max.dist)
    xindices <- c(xindices, tmp$xinds)
    yindices <- c(yindices, tmp$yinds)
  }

  xdat <- as.matrix(X$data[xindices,], length(xindices))
  ydat <- as.matrix(Y$data[yindices,], length(yindices))

  # Sometimes file reading (with csv at least) leads to situation where the last column is NA.
  # To avoid this and other cases, remove 'NA samples'.
  nainds <- (colMeans(is.na(xdat)) == 1 | colMeans(is.na(ydat)) == 1)
  if (sum(nainds) > 0) {
    xdat <- xdat[, !nainds]
    ydat <- ydat[, !nainds]
    warning(paste("Samples ", colnames(X$data)[nainds], " contained exclusively NA's; removed."))
  }
  
  # Convert chr to factor
  if ("23" %in% X$info[["chr"]]) {X$info[X$info[["chr"]] == "23", "chr"] <- "X"}
  if ("24" %in% X$info[["chr"]]) {X$info[X$info[["chr"]] == "24", "chr"] <- "Y"}

  if ("23" %in% Y$info[["chr"]]) {Y$info[Y$info[["chr"]] == "23", "chr"] <- "X"}
  if ("24" %in% Y$info[["chr"]]) {Y$info[Y$info[["chr"]] == "24", "chr"] <- "Y"}

  X$info$chr = factor(X$info$chr, levels = c(1:22, "X", "Y"))
  Y$info$chr = factor(Y$info$chr, levels = c(1:22, "X", "Y"))

  newX <- list(data = xdat, info = X$info[xindices,])
  newY <- list(data = ydat, info = Y$info[yindices,])

  # Check if arm info is missing

  list(X = newX, Y = newY)

}

closest <- function(a, vec){which.min(abs(a - vec))}

get.neighboring.probes <- function (X, Y, chr, max.dist, control.arms = TRUE) {

  xinds <- yinds <- c()
  
  # Use arm information if it is available and not blocked
  if (("arm" %in% names(X$info)) && ("arm" %in% names(Y$info)) && control.arms) {    

    for (arm in c('p', 'q')){      
      # Investigate specified arm
      xchrinds <- which(as.character(X$info$chr) == chr & X$info$arm == arm)
      ychrinds <- which(as.character(Y$info$chr) == chr & Y$info$arm == arm)
      inds <- get.neighs(X, Y, xchrinds, ychrinds, max.dist)
      xinds <- c(xinds, inds$xinds)
      yinds <- c(yinds, inds$yinds)
    }
  } else {
    # Investigate the whole chromosome
    xchrinds <- which(as.character(X$info$chr) == chr)
    ychrinds <- which(as.character(Y$info$chr) == chr)
    inds <- get.neighs(X, Y, xchrinds, ychrinds, max.dist)
    xinds <- c(xinds, inds$xinds)
    yinds <- c(yinds, inds$yinds)
  }

  # return the indices
  list(xinds = xinds, yinds = yinds)
  
}


get.neighs <- function (X, Y, xchrinds, ychrinds, max.dist) {

  xinds <- yinds <- NULL

  if ( length(xchrinds) > 0 && length(ychrinds) > 0 ){
    
    #Find indices of closest probe from Y for each from X
    xi <- 1:length(xchrinds)
    yi <- sapply(as.numeric(as.character(X$info$loc[xchrinds])), closest, vec = as.numeric(as.character(Y$info$loc[ychrinds])))

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


