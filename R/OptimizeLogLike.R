###This script contains functions necessary to optimize the log-likelihood of a N(mean, V(mean)), where V(mean) is the multinomial asymptotic variance
###The constraints are that the mean must lie WITHIN the simplex, NOT on the boundary. This is satisfied, since the log-like is -inf on the boundary

#This computes the gradient of the log-likelhood
#X in n x d, where d is the number of covariates
#Omega is d x K, where K + 1 is the number of cell types
#vec(Omega) = (w_1, ..., w_K)^T, where w_k is the ith column of Omega => grad(Omega) = (grad(w_1), ..., grad(w_K))^T
#In R, as.vector(Omega) = c(w_1, ..., w_K)
#C is n x k
#v is the variance (sigma^2)

CompLogLike <- function(X, C, Omega, v) {
	n <- nrow(X)
	d <- ncol(X)
	K <- ncol(Omega)
	v <- as.numeric(v)
	LL <- -n/2 * log(v) - 1/2*sum( log(X %*% Omega) ) - 1/2*sum( log(1 - apply(X %*% Omega, 1, sum)) )   #Determinant portion of log-likelihood
	for (i in 1:n) {
		x.i <- X[i,]
		c.i <- cbind(C[i,])
		p.i <- t(Omega) %*% cbind(x.i)
		LL <- LL - 1/(2*v) * t(c.i - p.i) %*% (diag(p.i) - p.i %*% t(p.i)) %*% (c.i - p.i)
	}
	return(LL)
}

CompGrad <- function(X, C, Omega, v) {
	n <- nrow(X)
	d <- ncol(X)
	K <- ncol(Omega)
	v <- as.numeric(v)
	
	s1 <- 0
	s2 <- 0
	grad <- rep(0, d*K)
	for (i in 1:n) {
		x.i <- X[i,]
		c.i <- C[i,]
		denom.i <- 1 - sum( t(Omega) %*% cbind(x.i) )
		s1 <- s1 + x.i / denom.i / 2
		s2 <- s2 - 1/(2*v) * ( sum(c.i) - 1 )^2 / denom.i^2 * x.i
		
		for (k in 1:K) {
			w.k <- Omega[,k]
			grad[(1 + (k-1)*d):(k*d)] <- grad[(1 + (k-1)*d):(k*d)] - ( 1/2 / sum(x.i*w.k) * x.i ) + ( 1/(2*v) * C[i,k]^2 / sum(x.i*w.k)^2 * x.i )
		}
	}
	grad <- grad + rep(s1 + s2, K)
	return(-grad)    #I want to MINIMIZE the -log-likelihood function
}

#Compute full Hessian of Omega (but not v)
#If H is NOT pd, I will update the eigen values of H so that it is pd
CompHessian <- function(X, C, Omega, v) {
	n <- nrow(X)
	d <- ncol(X)
	K <- ncol(Omega)
	v <- as.numeric(v)
	one.vec <- cbind(rep(1, K))
	delta <- 1e-8    #Make smallest eigenvalue of H delta if any eigenvalues of H are smaller than this
	
	tmp.mat <- X %*% Omega  #An n x K matrix
	row.sums.XO <- apply(tmp.mat, 1, sum)
	row.sums.C <- apply(C, 1, sum)
	
	R <- t(X) %*% diag( 1/(1 - row.sums.XO)^2/2 - (1-row.sums.C)^2 / (1 - row.sums.XO)^3/v ) %*% X
	H <- kronecker(one.vec %*% t(one.vec), R)
	for (k in 1:K) {
		ind.k <- (1 + (k-1)*d):(k*d)
		dot.k <- as.vector(X %*% cbind(Omega[,k]))    #x_i' w_k for i = 1, ..., n
		H[ind.k, ind.k] <- H[ind.k, ind.k] + t(X) %*% diag(1/2 * 1/dot.k^2 - 1/v * C[,k]^2 / dot.k^3) %*% X
	}
	H <- -H
	
	eigen.H <- eigen(H, symmetric = T)
	if (eigen.H$values[d*K] < delta) {
		Q <- eigen.H$vectors
		eigs <- eigen.H$values
		eigs[eigs < delta] <- delta
		H <- Q %*% diag(eigs) %*% t(Q)
	}
	return(H)
}

#Compute the Fisher information matrix, including the indices corresponding to sigma^2
CompFisher <- function(X, Omega, v) {
	n <- nrow(X)
	d <- ncol(X)
	K <- ncol(Omega)
	one.vec <- cbind(rep(1, K))
	
	denom <- 1 - as.vector(t(one.vec) %*% t(Omega) %*% t(X))
	rest <- 1/v - 1/2 + (1 - 1/v) * as.vector(t(one.vec) %*% t(Omega) %*% t(X))
	R <- t(X) %*% diag( 1/denom^2 * rest ) %*% X   #E d2l/dw_k dw_r', r \neq k, a d x d matrix
	H.0 <- kronecker(one.vec %*% t(one.vec), R)
	H <- matrix(0, nrow=d*K+1, ncol=d*K+1)
	H[1:(d*K), 1:(d*K)] = H.0
	
	for (k in 1:K) {
		tmp.vec.k <- as.vector(X %*% cbind(Omega[,k]))
		diag.k <- 1/tmp.vec.k^2 * (1/2 + (1/v - 1) * tmp.vec.k)
		H[(1 + (k-1)*d):(k*d), (1 + (k-1)*d):(k*d)] = H[(1 + (k-1)*d):(k*d), (1 + (k-1)*d):(k*d)] + t(X) %*% diag(diag.k) %*% X
		H[(d*K+1),1:(d*K)] = 1/v/2 * as.vector(t(X) %*% cbind( 1/as.vector(X %*% cbind(Omega[,k])) )) - 1/v/2 * as.vector( t(X) %*% cbind(1/denom) )
		H[1:(d*K), (d*K+1)] = H[(d*K+1),1:(d*K)]
	}
	H[d*K+1, d*K+1] <- n*K/(2*v^2)
	
	return(H)
}

#Check boundary conditions
#Return 1 if they are satisfied, 0 if not
CheckBound <- function(X, Omega, dir) {    ##It's assumed Omega is a matrix and direction is a vector
	tol <- 1e-3
	n <- nrow(X)
	K <- ncol(Omega)
	d <- nrow(Omega)
	dir.mat <- matrix(dir, nrow=d, ncol=K)
	
	mean.mat <- X %*% (Omega + dir.mat)    #All of the entries of this matrix should be between 0 and 1
	mean.mat.vec <- as.vector(mean.mat)      #All of the entries of this vector should be between 0 and 1
	mean.vec <- apply(mean.mat, 1, sum)      #All of the entries of this vector should be between 0 and 1
	
	if (sum( mean.mat.vec >= (1-tol) ) + sum(mean.mat.vec <= tol) + sum( mean.vec >= (1-tol) ) + sum(mean.vec <= tol)) {
		return(0)
	} else {
		return(1)
	}
}

#Compute the maximum likelihood estimate for Omega and v
#The inputs are X, C, Omega_start, v_start and a tolerance for the gradient
#The algorithm first updates Omega using the Fisher information insted of the exact Hessian and then updates v with an exact solution
#It performs a line search in the modified Newton direction
#If sigma^2 is known/fixed, set vfixed=T
MaxLike.NewtonLS <- function(X, C, Omega_start, v_start, tol, vfixed=F) {
	max.iter <- 1e5
	rho <- 1/2           ##This shrinks the direction until it is within the boundary. Note that if B_k is pd, the direction will always be a search direction
	n <- nrow(X)
	d <- ncol(X)
	K <- ncol(Omega_start)
	count <- 1
	norm.grad <- 1 + tol
	Omega.0 <- Omega_start
	v.0 <- v_start
	grad.vec <- rep(NA, max.iter)
	
	while (norm.grad > tol && count < max.iter) {
		
		##Update Omega##
		grad.0 <- CompGrad(X, C, Omega.0, v.0)
		norm.grad <- max(abs(grad.0))
		grad.vec[count] <- norm.grad
		H.0 <- CompHessian(X, C, Omega.0, v.0)
		dir.1 <- - solve(H.0, grad.0)       ##Omega is in this direction. I need to make sure the solution is within the boundary
		count.check <- 0
		
		while(!CheckBound(X, Omega.0, dir.1) && count.check < 1000) {
			dir.1 <- rho*dir.1
			count.check <- count.check + 1
		}
		Omega.1 <- Omega.0 + matrix(dir.1, nrow=d, ncol=K)
		
		##Update v##
		if (!vfixed) {
			A = 0
			for (i in 1:n) {
				x.i <- X[i,]
				c.i <- C[i,]
				p.i <- as.vector(t(Omega.1) %*% cbind(x.i))
				Sigma.i <- diag(p.i) - cbind(p.i) %*% rbind(p.i)
				A = A + rbind(c.i - p.i) %*% cbind( solve(Sigma.i, c.i - p.i) )
			}
			v.1 <- as.numeric(A/n/K)
			v.0 <- v.1
		}
		
		count <- count + 1
		Omega.0 <- Omega.1
	}
	
	return(list(Omega = Omega.0, I = CompFisher(X, Omega.0, v.0), v = v.0, infnorm.grad = norm.grad))
}

