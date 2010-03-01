
W.cca.EM = function (Dcov, M, beta) {
  
	# EM update for W when phi is assumed known
    	#beta <-M%*%t(W$total)%*%phi.inv$total
        #M <- solve(t(W)%*%phi.inv%*%W + I)
  Dcov$total%*%t(beta)%*%solve(M + beta%*%Dcov$total%*%t(beta))

}

