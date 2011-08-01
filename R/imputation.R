imputation <- function (X) {

  # Impute missing values from a Gaussian as a
  # standard simple imputation approach.

  nas <- (is.na(X) | is.nan(as.matrix(X)))      # indices of missing values

  if (sum(nas)>0) {
    for (i in 1:ncol(X)) {
      x <- X[, i]
      nas <- (is.na(x) | is.nan(x))      # indices of missing values
      if (mean(nas) < 1) {
        X[nas, i] <- rnorm(sum(nas), mean(x[!nas]), # replace missing
                         sd(x[!nas]))
      } else {
      	stop(paste("Column ", i, "has all values missing; check the data!"))
      }
    }
  } else {} # No missing vals
  
  X
				     
}
				         
