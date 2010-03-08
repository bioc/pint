# Perhaps this could be modified later?
update.phi.EM.simcca = function (Dcov, W.new, phi.inv, W.old, M) {

  # From BachJordan sec. 4.1

  dx <- ncol(Dcov$X)

  # Reduces to this when Wx = Wy:
  mat <- W.old$X%*%M%*%t(W.new$X)
  nullmat <- matrix(0,nrow=dx,ncol=dx)
  mat2 <- Dcov$total - Dcov$total%*%rbind(cbind(phi.inv$X%*%mat,nullmat), cbind(nullmat,phi.inv$Y%*%mat))


  # ensure that diagonals are nonnegative by adding a small positive constant
  dd = diag(mat2)
  dd[dd<1e-100] = 1e-100
  diag(mat2) = dd

  # Diagonal is regularized to avoid singluar matrix
  mat2 <- mat2 + (1e-3)*diag(ncol(Dcov$total))

  phi <- list()
  phi$total <- mat2
  phi$X <- phi$total[1:dx,1:dx]
  phi$Y <- phi$total[(dx+1):(2*dx),(dx+1):(2*dx)]

  phi

}
