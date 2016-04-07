###This script contains functions necessary to compute the first, estimator and the second quasi-likelihood estimator for beta_g

#This is the first estimator
#C.1 is n.1 x K, orthog.X.1 is n.1 x n.1, X.2 is d x n.2, Omega is d x K, l.g is a K-vector estimate for l_g, V.C.1t is the orthonormal basis for C.1' and is n.1 x n.1 - K
Estimator.1 <- function(var, C.1, V.C.1t, orthog.X.1, X.2, Omega, V, l.g, y.1, y.2) {
	n.1 <- nrow(C.1)
	K <- ncol(C.1)
	n.2 <- ncol(X.2)
	m <- n.1 - K + n.2
	U.g <- array(NA, dim=c(m, m))    #An m x m matrix
	l.g <- cbind(l.g)
	
	XO.2 <- t(X.2) %*% Omega		    #An n.2 x K matrix
	kron.l <- kronecker(diag(n.2), t(l.g))
	
	U.g[1:(n.1-K), 1:(n.1-K)] <- var * diag(n.1-K)
	U.g[1:(n.1-K), (n.1-K+1):m] <- -var * t(V.C.1t) %*% orthog.X.1 %*% C.1 %*% solve( t(C.1) %*% orthog.X.1 %*% C.1, t(XO.2) )
	U.g[(n.1-K+1):m, 1:(n.1-K)] <- t(U.g[n.1-K, (n.1-K+1):m])
	U.g[(n.1-K+1):m, (n.1-K+1):m] <- var * diag(n.2) + kron.l %*% V %*% t(kron.l) + var * XO.2 %*% solve( t(C.1) %*% orthog.X.1 %*% C.1, t(XO.2) )
	
	Y <- cbind( c(as.vector(t(V.C.1t) %*% cbind(y.1)), y.2 - as.vector(XO.2 %*% l.g)) )
	X <- rbind( t(V.C.1t) %*% t(X.1), t(X.2) )
	beta <- solve( t(X) %*% solve( U.g, X ), t(X) ) %*% solve( U.g, Y )
	var.beta <- solve( t(X) %*% solve( U.g, X ) )
	return(list(beta=beta, var.beta=var.beta))
}

#This is the second estimator
#The dimension of the input are the same as above
#theta.0 = (beta.0, l.0)
Estimator.2 <- function(var, C.1, X.1, X.2, Omega, V, y.1, y.2, theta.0, tol) {
	n.1 <- ncol(X.1)
	n.2 <- ncol(X.2)
	K <- ncol(C.1)
	d <- nrow(X.1)
	XO.2 <- t(X.2) %*% Omega   #X.2' Omega = An n.2 x K matrix
	SX1X1 <- X.1 %*% t(X.1)      #A d x d matrix
	SX1C1 <- X.1 %*% C.1         #A d x K matrix
	SC1C1 <- t(C.1) %*% C.1     #A K x K matrix
	
	U.mat <- rbind( cbind( X.1, X.2 ), cbind( t(C.1), t(XO.2) ) )      #U.k = U.mat stand.y
	T.var <- rbind( cbind( 1/var * SX1X1, 1/var * SX1C1 ), cbind( 1/var * t(SX1C1), 1/var * SC1C1 ) )    #T.k = T.var + T.G.k
	
	count <- 1
	grad.norm <- tol + 1
	while(count < 1e4 && grad.norm > tol) {
		beta.0 <- theta.0[1:d]
		l.0 <- theta.0[(d+1):(d+K)]
		mu1.0 <- t(X.1) %*% beta.0 + C.1 %*% l.0
		mu2.0 <- t(X.2) %*% beta.0 + XO.2 %*% l.0
		
		kron.l <- kronecker(diag(n.2), rbind(l.0))
		G.0 <- kron.l %*% V %*% t(kron.l) + var * diag(n.2)
		U.0 <- U.mat %*% c( 1/var * (y.1 - mu1.0), solve(G.0, y.2 - mu2.0) )
		grad.norm <- max(abs(U.0))
		
		X2GinvX2 <- X.2 %*% solve(G.0, t(X.2))
		X2GinvC2 <- X.2 %*% solve(G.0, XO.2)
		T.0 <- T.var + rbind( cbind( X2GinvX2, X2GinvC2 ), cbind( t(X2GinvC2), t(Omega) %*% X2GinvX2 %*% Omega ) )
		
		theta.0 <- theta.0 + solve(T.0, U.0)
		count <- count + 1
	}
	
	return(list( beta=theta.0[1:d], l=theta.0[(d+1):(d+K)], n.iter=count, inf.norm=grad.norm, var.theta=solve(T.0) ))
}

#Compute alpha and gamma for gamma distribution using ML
#Input is starting value and p-vector s = (y_1'y_1, ..., y_p'y_p)
#theta = (alpha, gamma) (here gamma := beta in normal Gamma distribution parametrization)
ML.Gamma <- function(theta.0, s, n, tol) {
	p <- length(s)
	tol.boundary <- 1e-8
	count <- 1
	norm.grad <- 10
	rho <- 3/4
	c <- 0.01
	while (norm.grad > tol && count < 1e4) {
		alpha.0 <- theta.0[1]
		beta.0 <- theta.0[2]
		grad.0 <- -p * c( log(beta.0) - digamma(alpha.0) + digamma(alpha.0 + n/2) - mean( log(1/2*s + beta.0) ), alpha.0/beta.0 - (alpha.0 + n/2) * mean( 1/(1/2*s + beta.0) ) )
		norm.grad <- max(abs(grad.0))
		H.0 <- -p * matrix( c( trigamma(alpha.0 + n/2) - trigamma(alpha.0), 1/beta.0 - mean(1/(1/2*s + beta.0)), 1/beta.0 - mean(1/(1/2*s + beta.0)), -alpha.0/beta.0^2 + (alpha.0 + n/2) * mean( 1/(1/2*s + beta.0)^2 ) ), nrow=2, ncol=2 )
		eig.0 <- eigen(H.0)$values
		if (eig.0[2] < 0) {
			u.0 <- -grad.0/sqrt(sum(grad.0 * grad.0))
		} else {
			u.0 <- -solve(H.0, grad.0)
		}
		
		theta.1 <- theta.0 + u.0
		while(theta.1[1] < tol.boundary || theta.1[2] < tol.boundary || -LogLike.Gamma(theta.1, s, n) > -LogLike.Gamma(theta.0, s, n) + c*sum(grad.0 * u.0)) {
			u.0 <- rho * u.0
			theta.1 <- theta.0 + u.0
		}
		theta.0 <- theta.1
		count <- count + 1
	}
	return(list( alpha=theta.0[1], gamma=theta.0[2], hessian=H.0, inf.norm=norm.grad, n.iter=count ))
}

#Compute the log-likelihood function
LogLike.Gamma <- function(theta, s, n) {
	alpha <- theta[1]
	beta <- theta[2]
	p <- length(s)
	return( p*alpha*log(beta) - p*lgamma(alpha) + p*lgamma(alpha + n/2) - alpha*sum( log(1/2*s + beta) ) )
}

#Compute s needed to maximize the Gamma log-likelihood
#X is n x d, Y is p x n
Compute.s <- function(X, Y) {
	p <- nrow(Y)
	n <- nrow(X)
	Q <- diag(n) - X %*% solve(t(X) %*% X, t(X))
	Y.tilde <- Y %*% Q
	s <- rep(0, p)
	for (g in 1:p) {
		y.tilde.g <- Y.tilde[g,]
		s[g] <- sum(y.tilde.g * y.tilde.g)
	}
	return(s)
}