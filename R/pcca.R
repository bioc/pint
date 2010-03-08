
pcca = function (X, Y, zDimension=NULL) {

  # replaces: solve.CCA.full
  
  # If zDimension given, then
  # only pick zDimension first principal components
  # and estimate marginals accordingly
  # relies on the fact that the principal components
  # can be estimated consecutively in pCCA
  
  # Add here centering of the data matrices X, Y
  # (center the dimensions to 0)
  
	Dcov = list()
	Dcov$X = cov(t(X))
	Dcov$Y = cov(t(Y))

	# Solve W
	W = solve.w(t(X),t(Y),Dcov$X,Dcov$Y,zDimension)

        # Then take only the zDimension first components if defined
        #ifelse(is.null(zDimension), zDimension <- ncol(W$X), zDimension <- zDimension)
	W$X = as.matrix(W$X)
	W$Y = as.matrix(W$Y)
	W$total = rbind(W$X,W$Y)
        
	# estimate
	phi = list()
	phi$X = Dcov$X - W$X%*%t(W$X)
	phi$Y = Dcov$Y - W$Y%*%t(W$Y)

        # This would retrieve principal canonical components from the prob.CCA model
        # assuming that W and phi are known (at least in the full-rank case)
        # U = solve.archambeau(X, Y, W$X, W$Y, phi$X, phi$Y)
        
	list(W=W, phi=phi)
}


