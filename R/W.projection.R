W.projection <- function(model,X,Y){

   z <- z.expectation(model,X,Y)
   W <- getW(model)$total
   Wx <- getW(model)$X

   # Calculate first component of PCA for W*z
   pca <- princomp(t(W%*%z))
   projvec <- pca$loadings[,1]

   # Divide to X and Y components
   projvecx <- projvec[(1:nrow(Wx))]
   projvecy <- projvec[-(1:nrow(Wx))]
   
   return(list(total = projvec, X = projvecx, Y = projvecy))

}

