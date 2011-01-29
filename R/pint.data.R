pint.data <- function(data, info, impute = TRUE, replace.inf = TRUE, remove.duplicates = TRUE){

  data2 <- as.data.frame(data)
  info2 <- as.data.frame(info)

  dimnames(data2) <- dimnames(data) 
  dimnames(info2) <- dimnames(info)
  
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
  inds <- is.infinite(data)
  if (replace.inf && any(inds)) {
     message("Replacing infinite values with highest non-infinite values seen in the data")
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
  if (is.null(info$arm)) {
    arm <- factor(rep("p",length(loc)),levels=c("p","q"))
  } else {
    arm <- info$arm
  }
  
  info2 <- data.frame(chr = factor(info$chr, levels = c(1:24)), arm = arm, loc = as.numeric(loc))
  rownames(info2) <- rownames(info)
  info <- info2

  # Order data by chr, arm and loc
  ord <- order(info$chr,info$arm,info$loc)
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

  list(data = data, info = info)
}
