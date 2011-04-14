fit.cgh.mrna.byname <- function(X, Y, geneName, windowSize, ...){
                    
  index <- which(rownames(X$data) == geneName)
  window <- fixed.window(X,Y,index,windowSize)
  fit.dependency.model(window$X,window$Y, ...)

}
