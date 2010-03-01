update.phi.EM <-
function (Dcov, W.new, phi.inv, W.old, M,nullmat) {

  # From BachJordan sec. 4.1	

  mat = phi.inv$total%*%W.old$total%*%M%*%t(W.new$total)
  mat2 = Dcov$total - Dcov$total%*%mat

  # ensure that diagonals are nonnegative by adding a small positive constant
  dd = diag(mat2)
  dd[dd<1e-100] = 1e-100
  diag(mat2) = dd

  # for large dimensionality, also regularize more the diagonal
  # 25 selected here since after that regularization was needed
  # the more dimensions, more regularization
  if (ncol(Dcov$X)>10) {
	#mat2 = mat2 + (ncol(Dcov$X)/(1000*51))*diag(ncol(Dcov$total))
	mat2 = mat2 + (1e-3)*diag(ncol(Dcov$total))
  }

  phi = list()	
  phi$total = mat2
  phi$X = phi$total[1:ncol(nullmat),1:ncol(nullmat)]
  phi$Y = phi$total[(ncol(nullmat)+1):(2*ncol(nullmat)),(ncol(nullmat)+1):(2*ncol(nullmat))]

  phi

}

