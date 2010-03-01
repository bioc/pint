optimize.fullcov <-
function (W, phi, Dim, Dcov, nullmat, epsilon = 1e-3, par.change = 1e6, cost.old = 1e6, mySeed=123) {

  set.seed(mySeed)

  cost.new = cost5(W$X, phi, Dcov)
  initcost = cost.new

  while (par.change>epsilon || par.change < 0) {

	cost.old = cost.new

	##########################################################################

	phi.inv = list()
  	phi.inv$X = solve(phi$X)
  	phi.inv$Y = solve(phi$Y)
  	phi.inv$total = rbind(cbind(phi.inv$X,nullmat),cbind(nullmat,phi.inv$Y))
        
  	# assuming Wx = Wy!
  	#M = solve(t(W$X)%*%(phi.inv$X + phi.inv$Y)%*%W$X + diag(ncol(W$X)))
	M = solve(t(W$total)%*%phi.inv$total%*%W$total + diag(ncol(W$total)))
  	beta = M%*%t(W$total)%*%phi.inv$total #set.beta.simcca(M, W, phi.inv) 

	###########################################################################

	# store 
	W.old = W

	## EM update
	##W$total = W.simcca.EM(Dcov, M, beta)
	#W$total = Dcov$total%*%t(beta)%*%solve(M + beta%*%Dcov$total%*%t(beta))
	#W$X = W$total[1:Dim$X,]
	#W$Y = W$total[-c(1:Dim$X),]
	## robust solution: take average of the two estimates
	#W$X = W$Y = (W$X+W$Y)/2
	#W$total = rbind(W$X,W$Y)


	# alternative update
	phihat = phi$X + phi$Y
	phihat.inv = solve(phihat)
	what = 2*W$X
	M.w = solve(t(what)%*%phihat.inv%*%what + diag(ncol(what)))
  	beta.w = M.w%*%t(what)%*%phihat.inv 

	# Update what. what = 2*w
	what = Dcov$sum%*%t(beta.w)%*%solve(M.w + beta.w%*%Dcov$sum%*%t(beta.w))
	w = what/2
	wlist = list(); wlist$X = wlist$Y = w; wlist$total = rbind(w,w)
	W = wlist

	W.new = W

	#cost.w = cost5(W$X, phi, Dcov)
	#print(paste("cost.w",cost.w))

	#############################################
		
	# Update phi 
	#phi = update.phi.fullcov(Dcov, wtw) #ML
	phi = update.phi.EM(Dcov, W.new, phi.inv, W.old, M, nullmat)

	#cost.phi = cost(W, phi, Dim, Dcov)
	#cost.phi = cost5(W$X, phi, Dcov)
	#print(paste("cost.phi",cost.phi))
	#print(phi.inv)

	###########################################################################
	
	# Check marginal likelihood (-logP) for data to keep eye on convergence correctness
	# the smaller, the better are the parameters
	cost.new = cost5(W$X, phi, Dcov)

        #print(cost.new)

	par.change = (cost.old - cost.new)#/cost.old			

   }
	
   #list(W=W, phi=phi, score=cost.new, walt=w, phialt=phihat)
   list(W=W, phi=phi, score=cost.new)
   
}

