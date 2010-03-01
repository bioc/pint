solve.w <-
function (X,Y) {

	# Solve W for CCA case using Archambeau06 equations

	cca <- cancor(t(X),t(Y))

	# center row means to zero
	# source("matrixOperations.R")
	X = t(centerData(t(X))) 
	Y = t(centerData(t(Y))) 

	Cxx = cov(t(X))
	Cyy = cov(t(Y))
	
	Ux = cca$xcoef
	Uy = cca$ycoef
	
	#Qx = diag(nrow(X)) # Note: only requirement is that
	#Qy = diag(cca$cor) # Qx%*%t(Qy) = canonical correlations
	
	if (nrow(X)>nrow(Y)) {
		Qx <- diag(1,nrow(X),nrow(Y))
		Qy <- diag(cca$cor)
	} else {
		Qx <- diag(cca$cor)
		Qy <- diag(1,nrow(Y),nrow(X))
	}

	# ML estimates for model W:
	Wx = as.matrix(Cxx%*%Ux%*%Qx)
	Wy = as.matrix(Cyy%*%Uy%*%Qy)

	list(X = Wx, Y = Wy)
}

