update.phi.EM2 <-
function (Dcov, W.new, phi.inv, W.old, M, nullmat) {

  # From BachJordan sec. 4.1	

  #mat <- W.old$total%*%M%*%t(W.new$total)
  #mat.x <- phi.inv$X%*%mat
  #mat.y <- phi.inv$Y%*%mat
  #foo <- rbind(cbind(mat.x, mat.x), cbind(mat.y, mat.y))
  #mat2 <- Dcov$total - Dcov$total%*%foo

  mat2 <- Dcov$total - Dcov$total%*%phi.inv$total%*%W.old$total%*%M%*%t(W.new$total)

  # ensure that diagonals are nonnegative by adding a small positive constant
  dd <- diag(mat2)
  dd[dd<1e-6] <- 1e-6
  diag(mat2) <- dd

  # for large dimensionality, regularize more the diagonal
  # the more dimensions, more regularization
  if (ncol(Dcov$X) > 5) {mat2 <- mat2 + (1e-2)*diag(ncol(Dcov$total))}

  phi <- list()	
  phi$total <- mat2
  phi$X <- as.matrix(phi$total[1:nrow(Dcov$X),1:nrow(Dcov$X)], nrow = nrow(Dcov$X))
  phi$Y <- as.matrix(phi$total[-seq(nrow(Dcov$X)), -seq(nrow(Dcov$X))], nrow = nrow(Dcov$Y))

  phi

}

