solve.w <- function (Xc, Yc, Cxx, Cyy, dz = NULL) {

  # assumes Xc, Yc : samples x features, zero-mean features
  # Cxx and Cyy are covariances from cov(Xc) and cov(Yc)
  # dz shows the desired rank of latent Z


  # NOTE: here the dimensions of Xc and Yc do not need to match
  # Note: in previous solve.w the input data was features x samples

  # Traditional CCA solution (modified from cancor function):
  nr <- nrow(Xc)

  qx <- qr(Xc)
  qy <- qr(Yc)
  dx <- qx$rank
  dy <- qy$rank

  z <- svd(qr.qty(qx, qr.qy(qy, diag(1, nr, dy)))[1L:dx, , drop = FALSE], dx, dy)

  xcoef <- backsolve((qx$qr)[1L:dx, 1L:dx, drop = FALSE], z$u)
  ycoef <- backsolve((qy$qr)[1L:dy, 1L:dy, drop = FALSE], z$v)

  #rownames(xcoef) <- colnames(Xc)[qx$pivot][1L:dx]
  #rownames(ycoef) <- colnames(Yc)[qy$pivot][1L:dy]
  #cca <- list(cor = z$d, xcoef = xcoef, ycoef = ycoef)

  # Solve W using Archambeau06 equations
  # Note: only requirement for Q is that Qx%*%t(Qy) = canonical correlations
  # Q corresponds to M (dz x dz) in Bach-Jordan 2005, p.8 (before sec 4.1)
  Qx <- diag(z$d[1:dz],dz,dz) # dz x dz matrix
  #Qy <- diag(1, nrow(Qx)) # also a dz x dz matrix: identity matrix -> omit

  # ML estimates for the prob. model W:
  dz <- ifelse(is.null(dz), length(z$d), dz)
  Wx <- Cxx%*%xcoef[,1:dz]%*%Qx
  Wy <- Cyy%*%ycoef[,1:dz]#%*%Qy # Qy is identity matrix -> omit

  list(X = Wx, Y = Wy)
}
