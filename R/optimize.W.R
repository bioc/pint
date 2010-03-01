optimize.W <-
function (W, T, phi, Dim, Dcov, priors, H, epsilon = 1e-6, par.change = 1e6, cost.old = 1e6, mySeed=123) {

  #We assume here that Wy = T%*%Wx. Optimize also T. 

  # Ensure that Wy = T * Wx:
  W$Y = T%*%W$X

  #par.old = c(unlist(W$X),unlist(T),unlist(phi))
  cost.new = cost.W(c(as.vector(W$X),as.vector(T)), phi, priors, Dim, Dcov, H)

  set.seed(mySeed)
  costs = c(cost.new)
  cnt = 1

  while (par.change>epsilon) {

	cost.old = cost.new

	##########################################################################

	  M = list()
	  M$X = set.M(W$X, phi$X)
	  M$Y = set.M(W$Y, phi$Y)

	  beta = list()
	  beta$X = set.beta(M$X, W$X, phi$X)
	  beta$Y = set.beta(M$Y, W$Y, phi$Y)

	###########################################################################

	  # Update W: initialize with previous W	

	  opt = optim(c(as.vector(W$X),as.vector(T)), cost.W, method = "L-BFGS-B", phi=phi, priors=priors, Dim=Dim, Dcov=Dcov, H = H, control = list(maxit=1e6),lower = -10*max(Dcov$total), upper = 10*max(Dcov$total))
	  #opt = optim(c(as.vector(W$X),as.vector(T)), cost.W, method = "BFGS", X=X, Y=Y, phi=phi, priors=priors, Dim=Dim, Dcov=Dcov, control = list(maxit=1e6))

	  # Convert optimized W parameter vector to actual matrices
	  wt = get.W(opt$par, Dim)
	  W = wt$W
	  T = wt$T


	###########################################################################
		
	  # Update phi

	  phi = update.phi(Dcov, M, beta, W, phi)
	  
	###########################################################################

	  cnt = cnt + 1
	
	  # Check and print marginal likelihood (-logP) for the data
	  # the smaller, the better are the parameters
	  cost.new = cost.W(c(as.vector(W$X),as.vector(T)), phi, priors, Dim, Dcov, H)
	  #print(log(cost.new))
	  #print(T)
	  par.change = (cost.old - cost.new)#/cost.old
          costs[[cnt]] = cost.new
	}
	
 	list(W=W, phi=phi, T=T, costs=costs, score=costs[[length(costs)]])
}

