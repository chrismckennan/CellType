##Do we underestimate Omega with GLS?

CompareEstimate <- function(X, alpha, Gamma, Sigma) {     #It is assumed X is n x d where the first column is a vector of 1's
	r <- ncol(Gamma)
	p <- nrow(Gamma)
	d <- ncol(X)
	n <- nrow(X)
	
	Q <- qr.Q(qr(X), complete=T)
	R <- qr.R(qr(X))
	Q1 <- Q[,1:d]
	Q2 <- Q[,(d+1):n]
	W <- matrix(rnorm(r*n), nrow=r, ncol=n)
	W1 <- W %*% Q1     #r x d
	E <- matrix(rnorm(n*p), nrow=p, ncol=n) * sqrt(Sigma)
	E1 <- E %*% Q1
	
	Y <- Gamma %*% (W + alpha %*% t(X[,2:d])) + E
	rm(E)
	out.em <- fa.em( t( Y %*% Q2 ), r=r )
	Gamma.hat <- out.em$Gamma
	Sigma.hat <- out.em$Sigma
	Info <- t(Gamma.hat / Sigma.hat) %*% Gamma.hat / p
	Rest <- t(t((Gamma) / Sigma.hat) %*% Gamma.hat / p)
	rm(out.em)
	
	alpha.hat <- solve(t(Gamma.hat / Sigma.hat) %*% Gamma.hat, t(Gamma.hat / Sigma.hat) %*% Y %*% Q1 %*% solve(t(R)))[,2:d]
	alpha.true <- alpha + (W1 %*% solve(t(R)))[,2:d]
	return(list(alpha.true=alpha.true, alpha.hat=alpha.hat, resid=solve(t(Gamma.hat / Sigma.hat) %*% Gamma.hat, t(Gamma.hat / Sigma.hat) %*% E1 %*% solve(t(R)))[,2:d], Info=Info, Rest=Rest))
}