test.segmented <- function(X){
  # How many identical values is needed before the function
  # return TRUE
  threshold.fraction = 0.33 
  threshold = floor(threshold.fraction*nrow(X)*ncol(X)) 

  # Number of identical values in the same column in adjacent rows 
  identicals = (X[1:(nrow(X)-1),] == X[2:nrow(X),])

  return(sum(identicals,na.rm=TRUE) > threshold)
}