dependency.score <- function (res) {

  W <- res$W
  phi <- res$phi

  if (!is.null(W$X)){
    # this equals to the trace of the full Phi
    noise <- sum(diag(as.matrix(phi$X))) + sum(diag(as.matrix(phi$Y))) 

    wtw <- rbind(W$X,W$Y)%*%t(rbind(W$X,W$Y))
  }
  else {
    noise <- sum(diag(as.matrix(phi$total)))

    wtw <- W$total %*% t(W$total)
  }

  signal = sum(diag(wtw)) # trace of full WWt covariance 

  cost = signal/noise
  cost
}

