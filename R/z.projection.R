z.projection <- function(model,X,Y){

   # Find original data
   index <- which(dimnames(X$data)[[1]] == getGeneName(model))
   window <- fixed.window(X, Y, index, getWindowSize(model))

   z <- z.expectation(model,X,Y)
   W <- getW(model)$total

   # Calculate first component of PCA for W*z
   pca = princomp(t(W%*%z))
   projvec = pca$loadings[,1]

   # Project data to this component
   data = rbind(window$X,window$Y)
   proj = t(data)%*%projvec
   
   return(proj)

   #barplot(t(proj),main="Probe data projection to first pca loading vector, zDim = 1 (pCCA)")
}

