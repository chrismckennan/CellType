##This function estimates beta using GEE
#We assume C.1 is K x n.1, X.i is d x n.i, Omega is K x d
#theta = (beta, l)
Estimator.2 <- function(var, C.1, X.1, X.2, Omega, V, y.1, y.2, theta.0, tol) {
	K <- nrow(C.1)
	n.1 <- ncol(C.1)
	n.2 <- ncol(X.2)
	d <- nrow(X.1)
	
	beta.0 <- theta.0[1:d]
	l.0 <- theta.0[(d+1):(d+K)]
}