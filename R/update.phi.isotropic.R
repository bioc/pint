
update.phi.isotropic <- function (Xcov, W, epsilon, dx) {

  # ML estimate of the variance given W
  # See Roweis, 'EM Algorithms for PCA and SPCA'
  # section 4 (SPCA; M-step for epsilon)

  # epsilon is the variance from previous
  # iteration step
  
  # auxiliary variables
  wtw <- W%*%t(W)
  M <- set.M.isotropic(wtw, epsilon, dx)

  # Calculate updated phi (= epsilon) and return
  sum(diag(Xcov - wtw%*%M%*%Xcov))/dx

}
