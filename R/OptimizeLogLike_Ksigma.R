###This script contains functions necessary to optimize the log-likelihood of a N(mean, V(mean)), where V(mean) is the multinomial asymptotic variance
###The constraints are that the mean must lie WITHIN the simplex, NOT on the boundary. This is satisfied, since the log-like is -inf on the boundary

#This computes the gradient of the log-likelhood
#X in n x d, where d is the number of covariates
#Omega is d x K, where K + 1 is the number of cell types
#vec(Omega) = (w_1, ..., w_K)^T, where w_k is the ith column of Omega => grad(Omega) = (grad(w_1), ..., grad(w_K))^T
#In R, as.vector(Omega) = c(w_1, ..., w_K)
#C is n x K
#sigma = (sigma_1, ..., sigma_K) is the standard deviation
#The difference between this script and OptimizeLogLike.R is this script allows for multiple variances (i.e K different dispersion parameters)

CompLogLike_Ksigma <- function(X, C, Omega, sigma) {   #sigma is a K-vector. This script returns the log-likelihood
	n <- nrow(X)
	d <- ncol(X)
	K <- ncol(Omega)
	sigma <- as.vector(sigma)
	LL <- -n * sum(log(sigma)) - 1/2 * sum(log(X %*% Omega)) - 1/2 * sum(log( 1 - apply(X %*% Omega, 1, sum) ))
	for (i in 1:n) {
		x.i <- X[i,]
		c.i <- C[i,]
		p.i <- as.vector(t(Omega) %*% cbind(x.i))
		LL <- LL - 1/2 * ( sum( c.i^2/sigma^2/p.i ) - 2*sum( c.i/sigma^2 ) + sum(p.i / sigma^2) + sum( (c.i - p.i)/sigma )^2/(1 - sum(p.i)) )
	}
	return(LL)
}

CompGrad_Ksigma <- function(X, C, Omega, sigma) {         #Returns the gradient of the -LL
	n <- nrow(X)
	d <- ncol(X)
	K <- ncol(Omega)
	sigma <- as.vector(sigma)
	XO <- X %*% Omega   		#A n x K matrix of expected cell types
	XO.sum <- apply(XO, 1, sum)			#An n-vector with components 1'Omega'x_i
	A.vec <- as.vector(( C - XO ) %*% cbind( 1/sigma ))/(1-XO.sum)      #An n-vector
	grad <- rep(0, d*K)
	for (k in 1:K) {
		sigma.k <- sigma[k]
		ind <- (1 + (k-1)*d):(k*d)
		XO.k <- XO[,k]
		c.k <- C[,k]
		grad[ind] <- -1/2 * as.vector(rbind(1/XO.k) %*% X) + 1/2 * as.vector(rbind(1/(1 - XO.sum)) %*% X) + 1/2 * as.vector(rbind( c.k^2/sigma.k^2/XO.k^2 ) %*% X) - 1/2/sigma.k^2 * as.vector(rbind( rep(1,n) ) %*% X) - 1/2 * as.vector(rbind(A.vec^2) %*% X) + 1/sigma.k * as.vector(rbind(A.vec) %*% X)
	}	
	return(-grad)
}

CompHessian_Ksigma <- function(X, C, Omega, sigma) {          #Compute the Hessian for -LL, only including the Omega terms. The output is a Kd x Kd matrix
	n <- nrow(X)
	d <- ncol(X)
	K <- ncol(Omega)
	sigma <- as.vector(sigma)
	one.vec <- cbind(rep(1,K))
	delta <- 1e-8
	
	XO <- X %*% Omega			#An n x K matrix
	XO.sum <- apply(XO, 1, sum)			#An n-vector with components 1'Omega'x_i
	func.C <- as.vector((C - XO) %*% cbind(1/sigma))
	R <- t(X) %*% diag( 1/2/(1-XO.sum)^2 - func.C^2/(1-XO.sum)^3 ) %*% X
	H <- kronecker(one.vec %*% t(one.vec), R)
	for (s in 1:K) {
		ind.s <- (1 + (s-1)*d):(s*d)
		sigma.s <- sigma[s]
		for (r in s:K) {
			ind.r <- (1 + (r-1)*d):(r*d)
			sigma.r <- sigma[r]
			if (s == r) {
				c.s <- C[,s]
				XO.s <- XO[,s]
				H[ind.s, ind.r] <- H[ind.s, ind.r] + t(X) %*% diag( 1/2/XO.s^2 - c.s^2/sigma.s^2/XO.s^3 + 2*func.C/sigma.s/(1-XO.sum)^2 - 1/sigma.s^2/(1-XO.sum) ) %*% X
			} else {
				add.rs <- t(X) %*% diag( (1/sigma.s + 1/sigma.r) * func.C/(1-XO.sum)^2 - 1/sigma.s/sigma.r/(1-XO.sum) ) %*% X
				H[ind.s, ind.r] <- H[ind.s, ind.r] + add.rs
				H[ind.r, ind.s] <- H[ind.r, ind.s] + add.rs
			}
		}
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

UpdateSigma <- function(X, Omega, C, sigma) {				#Update all indices of sigma and compute gradient after update
	n <- nrow(X)
	d <- ncol(X)
	K <- ncol(Omega)
	sigma <- as.vector(sigma)

	XO <- X %*% Omega			#An n x K matrix
	XO.sum <- apply(XO, 1, sum)			#An n-vector with components 1'Omega'x_i
	
	for (s in 1:K) {
		XO.s <- XO[,s]
		diff.s <- C[,s] - XO.s
		func.cs <- as.vector((C[,-s] - XO[,-s]) %*% cbind(1/sigma[-s]))
		
		sigma[s] <- quad.solve(n, -sum( diff.s/(1-XO.sum)*func.cs ), -sum( (1/XO.s + 1/(1-XO.sum)) * diff.s^2 ))
	}
	
	grad <- rep(0, K)
	for (s in 1:K) {
		sigma.s <- sigma[s]
		XO.s <- XO[,s]
		diff.s <- C[,s] - XO.s
		func.cs <- as.vector((C[,-s] - XO[,-s]) %*% cbind(1/sigma[-s]))
		grad[s] <- 1/sigma.s * ( -n + sum( diff.s^2 * (1/sigma.s^2/XO.s + 1/sigma.s^2/(1-XO.sum)) + diff.s/sigma.s/(1-XO.sum) * func.cs ) )		
	}
	return(list(sigma=sigma, grad=grad))
}

quad.solve <- function(a,b,c) {
	return( (-b + sqrt(b^2 - 4*a*c))/2/a )
}


CompFI_Ksigma <- function(X, C, Omega, sigma) {   #Computes the (d+1)K x (d+1)K Fisher Information matrix
	n <- nrow(X)
	d <- ncol(X)
	K <- ncol(Omega)
	sigma <- as.vector(sigma)
	XO <- X %*% Omega			#An n x K matrix
	XO.sum <- apply(XO, 1, sum)			#An n-vector with components 1'Omega'x_i	
	
	H <- matrix(0, nrow=(d+1)*K, ncol=(d+1)*K)
	for (s in 1:K) {
		sigma.s <- sigma[s]
		ind.s <- (1 + (s-1)*d):(s*d)
		for (r in s:K) {
			sigma.r <- sigma[r]
			ind.r <- (1 + (r-1)*d):(r*d)
			if (s == r) {
				##d2l/dw_r dw_s##
				H[ind.s, ind.r] <- t(X) %*% diag( 1/(1-XO.sum)^2 * ( 1/sigma.s/sigma.r - 1/2 + (1-1/sigma.s/sigma.r)*XO.sum ) + 1/XO[,s]^2 * (1/2 + (1/sigma.s^2 - 1)*XO[,s]) ) %*% X
				
				##d2l/dw_r dsigma_s##
				H[ind.s, d*K+r] <- as.vector( rbind( (1-XO[,s])/sigma.s/XO[,s] - XO[,s]/sigma.s/(1-XO.sum) ) %*% X )
				H[d*K+r, ind.s] <- H[ind.s, d*K+r]
				
				##d2l/dsigma_r dsigma_s##
				H[d*K+r, d*K+r] <- -n/sigma.s^2 + 1/sigma.s^2 * sum( 3*(1-XO[,s]) + XO[,s]*(3 - XO[,s] - 2*XO.sum)/(1-XO.sum) )
			} else {
				##d2l/dw_r dw_s##
				H[ind.s, ind.r] <- t(X) %*% diag( 1/(1-XO.sum)^2 * ( 1/sigma.s/sigma.r - 1/2 + (1-1/sigma.s/sigma.r)*XO.sum ) ) %*% X
				H[ind.r, ind.s] <- H[ind.s, ind.r]
				
				##d2l/dw_r dsigma_s##
				H[d*K+s, ind.r] <- -1/sigma.s * as.vector( rbind( XO[,s]/(1-XO.sum) ) %*% X )
				H[ind.r, d*K+s] <- H[d*K+s, ind.r]
				
				H[d*K+r, ind.s] <- -1/sigma.r * as.vector( rbind( XO[,r]/(1-XO.sum) ) %*% X )
				H[ind.s, d*K+r] <- H[d*K+r, ind.s]
				
				##d2l/dsigma_r dsigma_s##
				H[d*K+r, d*K+s] <- - sum( XO[,s]*XO[,r]/(1-XO.sum)/sigma.s/sigma.r )
				H[d*K+s, d*K+r] <- H[d*K+r, d*K+s]
			}
		}
	}
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


#Compute the ML estiamte for Omega and sigma = (sigma_1, ..., sigma_K)
MaxLike.NewtonLS_Ksigma <- function(X, C, Omega_start, sigma_start, tol) {
	max.iter <- 1e5
	rho <- 1/2           ##This shrinks the direction until it is within the boundary. Note that if B_k is pd, the direction will always be a search direction
	n <- nrow(X)
	d <- ncol(X)
	K <- ncol(Omega_start)
	count <- 1
	norm.grad <- 1 + tol
	Omega.0 <- Omega_start
	sigma.0 <- sigma_start
	grad.0 <- 0
	
	while (norm.grad > tol && count < max.iter) {
		
		##Update Omega##
		if (length(grad.0) < 2) {
			grad.0 <- CompGrad_Ksigma(X, C, Omega.0, sigma.0)
		}
		H.0 <- CompHessian_Ksigma(X, C, Omega.0, sigma.0)
		dir.1 <- - solve(H.0, grad.0)
		count.check <- 0
		
		while(!CheckBound(X, Omega.0, dir.1) && count.check < 1000) {
			dir.1 <- rho*dir.1
			count.check <- count.check + 1
		}
		Omega.1 <- Omega.0 + matrix(dir.1, nrow=d, ncol=K)	
		
		##Update sigma##
		up.sigma <- UpdateSigma(X, Omega.1, C, sigma.0)
		sigma.1 <- up.sigma$sigma
		
		##Check size of the gradient##
		grad.1 <- CompGrad_Ksigma(X, C, Omega.1, sigma.1)
		norm.grad <- max(abs(up.sigma$grad), abs(grad.1))
		
		##Re-initialize##
		Omega.0 <- Omega.1
		sigma.0 <- sigma.1
		grad.0 <- grad.1
		
		count <- count + 1
	}
	
	return(list(Omega=Omega.0, sigma=sigma.0, infnorm.grad=norm.grad, I=CompFI_Ksigma(X,C,Omega.0,sigma.0), n.iter = count))	
}