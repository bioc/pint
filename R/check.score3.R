check.score3 <-
function (W,phi) {

	# Expectation of the hidden variable(s) for each sample
	#zexp = z.expectation.fullcov(W, phi, X, Y)

	# this equals to the trace of the full Phi
	noise = sum(diag(as.matrix(phi$X))) + sum(diag(as.matrix(phi$Y))) 


	#if(is.null(W$WWt)){
	#	wtw = W$total%*%t(W$total) # shared.covariance
	#}
	#else{
	#	wtw = W$WWt
	#}
	wtw = rbind(W$X,W$Y)%*%t(rbind(W$X,W$Y))


	#signal = sum(eigen(wtw)$values)
	signal = sum(diag(wtw)) # trace of full WWt covariance 

	cost = signal/noise

	# Trace(W * t(W)) / Trace(Phi)
	#cost = sum(diag(wtw))/sum(diag(phi$total))

	#cost = sum(eigen(wtw)$values)/noise
	#cost = eigen(wtw)$values[[1]]/noise

	cost
}

