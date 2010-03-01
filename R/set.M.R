set.M <-
function (W, phi) {solve(t(W)%*%W/phi + diag(ncol(W)))}

