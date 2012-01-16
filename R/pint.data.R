pint.data <- function(data, info, impute = TRUE, replace.inf = TRUE, remove.duplicates = TRUE){

  data2 <- as.data.frame(data)
  info2 <- as.data.frame(info)

  if (is.null(colnames(data))) {
    warning("No column names for data matrix given, naming by indices..")
    colnames(data) <- as.character(1:ncol(data))
  }
  
  if (is.null(rownames(data))) {
    warning("No row names for data matrix given, naming by indices..")
    rownames(data) <- as.character(1:nrow(data))
  }  

  colnames(data2) <- colnames(data) 
  colnames(info2) <- colnames(info)
  
  if (any(duplicated(rownames(data)))) {
    apply(cbind(rownames(data), 1:nrow(data)), 1, function (x) {paste(x, collapse = "-row")})
  }
  
  if (any(duplicated(rownames(info)))) {
    apply(cbind(rownames(data), 1:nrow(data)), 1, function (x) {paste(x, collapse = "-row")})
  }
  
  data <- data2
  info <- info2

  ## Replace synonymous info fields
  # chromosome -> chr
  if ("chromosome" %in% colnames(info) && !"chr" %in% colnames(info)) {
    colnames(info)[which(colnames(info) == "chromosome")] <- "chr"
  }
  if ("stop" %in% colnames(info) && !"end" %in% colnames(info)) {
    colnames(info)[which(colnames(info) == "stop")] <- "end"
  }

  # If same number of rows and columns for data and info fields, 
  # assume that they match between data and info fields and give both
  # same rownames
  if (nrow(data) == nrow(info) && !rownames(data) == rownames(info)) {
     warning("The data and info fields had different rownames; using the rownames from data field for both.")
     rownames(info) <- rownames(data)
  }

  # First, provide information in a form that has corresponding rows in
  # data and info matrices (only take those available in both)
  coms <- intersect(rownames(data), rownames(info))
  coms <- setdiff(coms, c(""))
  data <- data[coms,]
  info <- info[coms,]
  info[["chr"]] <- as.character(info[["chr"]])

  # X/Y chromosome naming
  if ("X" %in% info[["chr"]]) {
    message("Changed chromosome name X to 23 for compatibility.")
    info[info[["chr"]] == "X", "chr"] <- "23"
  }
  if ("Y" %in% info[["chr"]]) {
    message("Changed chromosome name Y to 24 for compatibility.")
    info[info[["chr"]] == "Y", "chr"] <- "24"
  }

  # Quarantee that there are no duplicated rows (probes) in the data
  if (remove.duplicates){
    dupl <- duplicated(data)
    if (any(dupl)) {
      message("Removing duplicate probe signals..\n")
      data <- data[!dupl, ]
      info <- info[!dupl, ]
    }
  }
  
  # Impute missing values
  if (impute) {
    message("Imputing missing values..")
    data <- imputation(data)
  }
  
  # Replace infinite values by highest possible seen in the data
  inds <- is.infinite(as.matrix(data))
  if (replace.inf && any(inds)) {     message("Replacing infinite values with highest non-infinite values  in the data")
     data[inds] <- sign(data[inds])*max(abs(data[!inds])) # note the sign
     message(paste("...", 100*mean(inds), "percent of the values replaced."))
  }
  
  # Location
  if (is.null(info$loc) && is.null(info$bp)) {
    loc <- (info$start + info$end)/2
  } else if (!is.null(info$bp)) {
    loc <- info$bp
  } else {
    loc <- info$loc
  }
 
  # Arm
  if (!is.null(info$arm)) {
    info2 <- data.frame(chr = info$chr, arm = info$arm, loc = as.numeric(loc))
    rownames(info2) <- rownames(info)
    info <- info2
    # Order data by chr, arm and loc
    ord <- order(as.numeric(info$chr),info$arm,info$loc)
  } else {  # Arm info missing
    info2 <- data.frame(chr = as.numeric(info$chr), loc = as.numeric(loc))
    rownames(info2) <- rownames(info)
    info <- info2
    # Order data by chr and loc
    ord <- order(info$chr,info$loc)
  }  
  
  data <- data[ord,]
  info <- info[ord,]


  # If location information ('loc') for the probe is missing
  # but start and end positions are available, use them to calculate
  # probe middle bp location
  if (!"loc" %in% colnames(info)) {
    # sometimes factors, sometimes numerics given;
    # this should handle both cases correctly
    info[["loc"]] <- (as.numeric(as.character(info[, "start"])) + as.numeric(as.character(info[, "end"])))/2
  }

  list(data = as.matrix(data), info = info)
}
