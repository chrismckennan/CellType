Lambda.tmp <- Lambda
	A <- Z1 - L %*% F.mat
	SinvL <- L / Sigma
	SinvG <- Gamma / Sigma
	SinvA <- A / Sigma
	SinvZ <- Z2 / Sigma
	
	LtSinvZ <- t(SinvL) %*% Z2
	GtSinvZ <- t(SinvG) %*% Z2
	LtSinvG <- t(SinvL) %*% Gamma
	GtSinvG <- t(Gamma) %*% SinvG
	LtSinvL <- t(L) %*% SinvL
	AtSinvG <- t(A) %*% SinvG
	mid.r <- diag(r) + GtSinvG
	
	##Compute H matrices##
	H1 <- LtSinvZ - LtSinvG %*% solve(mid.r, GtSinvZ)           #K x m2
	H2 <- LtSinvL - LtSinvG %*% solve(mid.r, t(LtSinvG))        #K x K
	H3 <- AtSinvG - AtSinvG %*% solve(mid.r, GtSinvG)         #m1 x r
	H4 <- t(GtSinvZ) - t(GtSinvZ) %*% solve(mid.r, GtSinvG)  #m2 x r
	H5 <- LtSinvG - LtSinvG %*% solve(mid.r, GtSinvG)         #K x r
	H6 <- GtSinvG - GtSinvG %*% solve(mid.r, GtSinvG)       #r x r 
	
	##Compute G matrices##
	mid.Lambda <- solve(Lambda) + H2
	G1 <- H1 - H2 %*% solve(mid.Lambda, H1)         #K x m2
	G2 <- H2 - H2 %*% solve(mid.Lambda, H2)         #K x K
	G3 <- H3															      #m1 x r
	G4 <- H4 - t(H1) %*% solve(mid.Lambda, H5)     #m2 x r
	G5 <- H5 - H2 %*% solve(mid.Lambda, H5)        #K x r
	G6 <- H6																  #r x r
	G7 <- H6 - t(H5) %*% solve(mid.Lambda, H5)     #r x r

rho <- 3/4
c <- 0.01
for (i in 1:100) {
	SSG1inG2 <- solve(G2, G1)
	tmp.grad <- G2 - 1/m2 * G1 %*% t(G1)
	s.dir <- 1/m2 * SSG1inG2 %*% t(SSG1inG2) - solve(G2)
	mll.0 <- mll.lambda(Lambda.tmp, m2, H1, H2)
	deriv.0 <- trace(s.dir %*% tmp.grad)
	if (max(abs(tmp.grad)) < 1e-6 || deriv.0 > 0) {
		break
	}
	alpha <- 1
	count <- 1
	if ( mll.lambda(Lambda.tmp + alpha * s.dir, m2, H1, H2) - mll.0 > c * alpha * deriv.0 ) {
		alpha <- rho * alpha
		count <- count + 1
	}
	Lambda.tmp <- Lambda.tmp + alpha * s.dir
	mid.Lambda <- solve(Lambda.tmp) + H2
	G1 <- H1 - H2 %*% solve(mid.Lambda, H1)         #K x m2
	G2 <- H2 - H2 %*% solve(mid.Lambda, H2)         #K x K
}


mll.lambda <- function(lambda, m2, H1, H2) {
	mid.lambda <- solve(lambda) + H2
	return( m2/2 * ( log(det(lambda)) + log(det(mid.lambda)) ) - 1/2 * sum( solve(mid.Lambda, H1) * H1 ) )
}

trace <- function(x) {
	return(sum(diag(x)))
}

fun.bfgs <- function(R.vec, H1, H2, m2) {    #minus log-likelihood
	R <- matrix(R.vec, nrow=nrow(H2), ncol=ncol(H2), byrow=F)
	lambda <- R %*% t(R)
	mid.lambda <- solve(lambda) + H2
	return( m2/2 * ( log(det(lambda)) + log(det(mid.lambda)) ) - 1/2 * sum( solve(mid.Lambda, H1) * H1 ) )
}

grad.bfgs <- function(R.vec, H1, H2, m2) {
	R <- matrix(R.vec, nrow=nrow(H2), ncol=ncol(H2), byrow=F)
	lambda <- R %*% t(R)
	mid.lambda <- solve(lambda) + H2
	g1 <- H1 - H2 %*% solve(mid.lambda, H1)
	ssg1 <- g1 %*% t(g1)
	g2 <- H2 - H2 %*% solve(mid.lambda, H2)
	
	grad.mat <- m2/2 * (g2 %*% R + t(R) %*% g2) + 1/2 * (ssg1 %*% R + t(R) %*% ssg1)
	return(c(grad.mat))
}

tmp.eig <- eigen(Lambda.tmp)
R.vec.start <- c( tmp.eig$vectors %*% diag(sqrt(tmp.eig$values)) %*% t(tmp.eig$vectors) )
out.optim <- optim(R.vec.start, fn=fun.bfgs, gr=grad.bfgs, method="BFGS", H1=H1, H2=H2, m2=m2)




#####    #####

Proj.F <- t(F.mat) %*% solve(F.mat %*% t(F.mat)) %*% F.mat
dir.v <- t(svd(Mean.est %*% Proj.F)$v[,1]) %*% t(F.mat) %*% solve(F.mat %*% t(F.mat))


#Test T, L updates for convergence in the simplest case#
L <- cbind(v.0.1)
Gamma <- Gamma.1
	SinvL <- L / Sigma
	SinvG <- Gamma / Sigma
	SinvA <- A / Sigma
	SinvZ <- Z2 / Sigma
	
	LtSinvZ1 <- t(SinvL) %*% Z1
	LtSinvZ <- t(SinvL) %*% Z2
	GtSinvZ1 <- t(SinvG) %*% Z1
	GtSinvZ <- t(SinvG) %*% Z2
	LtSinvG <- t(SinvL) %*% Gamma
	GtSinvG <- t(Gamma) %*% SinvG
	LtSinvL <- t(L) %*% SinvL
	AtSinvG <- t(A) %*% SinvG
	mid.r <- diag(r) + GtSinvG
	
	##Compute H matrices##
	H1 <- LtSinvZ - LtSinvG %*% solve(mid.r, GtSinvZ)           #s x m2
	H2 <- LtSinvL - LtSinvG %*% solve(mid.r, t(LtSinvG))        #s x s
	H3 <- AtSinvG - AtSinvG %*% solve(mid.r, GtSinvG)         #m1 x r
	H4 <- t(GtSinvZ) - t(GtSinvZ) %*% solve(mid.r, GtSinvG)  #m2 x r
	H5 <- LtSinvG - LtSinvG %*% solve(mid.r, GtSinvG)         #s x r
	H6 <- GtSinvG - GtSinvG %*% solve(mid.r, GtSinvG)       #r x r 
	H7 <- LtSinvZ1 - LtSinvG %*% solve(mid.r, GtSinvZ1)           #s x m2

T.0 <- t(solve(H2, H7) %*% t(F.mat) %*% solve(SFF))
norm.0 <- sqrt(sum(T.0 * T.0))
T.0 <- T.0 / norm.0
L <- L * norm.0

#Follow coordinate 5 and 9

n.iter <- 200
vec.5 <- rep(0,200)
vec.9 <- rep(0,200)
for (i in 2:n.iter) {
	SinvL <- L / Sigma
	LtSinvZ1 <- t(SinvL) %*% Z1
	LtSinvZ <- t(SinvL) %*% Z2
	LtSinvL <- t(L) %*% SinvL
	LtSinvG <- t(SinvL) %*% Gamma
	H1 <- LtSinvZ - LtSinvG %*% solve(mid.r, GtSinvZ)           #s x m2
	H2 <- LtSinvL - LtSinvG %*% solve(mid.r, t(LtSinvG)) 
	H7 <- LtSinvZ1 - LtSinvG %*% solve(mid.r, GtSinvZ1)
	
	T.0 <- t(solve(H2, H7) %*% t(F.mat) %*% solve(SFF))
	norm.0 <- sqrt(sum(T.0*T.0))
	T.0 <- T.0 / norm.0
	new.F <- t(T.0) %*% F.mat
	L <- Z1 %*% t(new.F) %*% solve(new.F %*% t(new.F))
	vec.5[i] <- L[5,1]
	vec.9[i] <- L[9,1]
	
}

out.op <- optim(T.0, fn=func.T, gr=grad.T, method="BFGS", Lambda=Lambda.tmp, F.mat=F.mat, H1=H1, H2=H2, H7=H7, m2=m2.sim)

grad.T <- function(T.vec, F.mat, Lambda, H1, H2, H7, m2) {
	K <- nrow(F.mat); s <- length(T.vec)/K
	T.mat <- matrix(T.vec, nrow=s, ncol=K)
	mid.Lambda <- solve(T.mat %*% Lambda %*% t(T.mat)) + H2
	G1 <- H1 - H2 %*% solve(mid.Lambda, H1)
	G2 <- H2 - H2 %*% solve(mid.Lambda, H2)
	grad.mat <- m2 * G2 %*% T.mat %*% Lambda - G1 %*% t(G1) %*% T.mat %*% Lambda + H2 %*% T.mat %*% F.mat %*% t(F.mat) - H7 %*% t(F.mat)
	return(c(grad.mat))
}

func.T <- function(T.vec, F.mat, Lambda, H1, H2, H7, m2) {
	K <- nrow(F.mat); s <- length(T.vec)/K
	T.mat <- matrix(T.vec, nrow=s, ncol=K)
	Lambda.tilde <- T.mat %*% Lambda %*% t(T.mat)
	mid.Lambda.tilde <- solve(Lambda.tilde) + H2
	F.tilde <- T.mat %*% F.mat
	return( m2/2 * log(det(Lambda.tilde)) + m2/2 * log(det(mid.Lambda.tilde)) - sum( rowSums(H7 * F.tilde) ) + 1/2 * sum( rowSums(H2 * (F.tilde %*% t(F.tilde))) ) - 1/2 * sum( rowSums( (t(H1) %*% solve(mid.Lambda.tilde)) * t(H1) ) ) )
}



#####Can CCA pick up interesting rows of F?#####
p <- 3e5
r <- 4
K <- 3
n <- 140
Sigma <- rep(0.15, p)
Gamma <- matrix(rnorm(4*p), nrow=p, ncol=4) * 0.15/2
F.mat <- matrix(rnorm(K*n), nrow=n, ncol=K)
L <- cbind(rep(0,p), rep(0,p), rnorm(p) * 0.2 * rbinom(p, size=1, prob=0.3))
Z <- L %*% t(F.mat) + Gamma %*% matrix(rnorm(r * n), nrow=r, ncol=n) + matrix(rnorm(n*p), nrow=p, ncol=n) * sqrt(Sigma)

S <- 1/n * t(F.mat) %*% F.mat
S.half <- eigen(S)$vectors %*% diag(sqrt(eigen(S)$values), nrow=K, ncol=K) %*% t(eigen(S)$vectors)
mid.mat <- t(Z) %*% (Z / Sigma) - t(Z) %*% (Gamma / Sigma) %*% solve( diag(r) + t(Gamma) %*% (Gamma / Sigma) ) %*% t(Gamma) %*% (Z / Sigma)

solve(S.half, svd(solve(S.half) %*% t(F.mat) %*% mid.mat %*% F.mat %*% solve(S.half))$v[,1])


#####Using BFGS to solve for Gamma#####
ind.bfgs <- (1:100)
Gamma.1 <- Gamma[ind.bfgs,]
par.start <- c(Gamma.1)
out.bfgs <- optim(par=par.start, fn=Gamma.func.bfgs, gr=Gamma.bfgs.grad, method="BFGS", A=Z1.sim-L%*%F.mat, Z2=Z2.sim, Gamma.rest=Gamma[-ind.bfgs,], L=L, Sigma=Sigma.use, Lambda=Lambda.0, m1=m1.sim, m2=m2.sim, r=ncol(Gamma))


Gamma.bfgs.grad <- function(par, A, Z2, Gamma.rest, L, Sigma, Lambda, m1, m2, r) {
	p <- nrow(L)
	p1 <- length(par) / r
	
	Gamma.1 <- matrix(par, nrow=p1, ncol=r)
	Gamma.tilde.1 <- Gamma.1 / Sigma[1:p1]
	L.tilde.1 <- L[1:p1,] / Sigma[1:p1]
	Gamma <- rbind(Gamma.1, Gamma.rest)
	
	SinvL <- L / Sigma
	SinvG <- Gamma / Sigma
	SinvZ <- Z2 / Sigma
	
	LtSinvZ <- t(SinvL) %*% Z2
	GtSinvZ <- t(SinvG) %*% Z2
	LtSinvG <- t(SinvL) %*% Gamma
	GtSinvG <- t(Gamma) %*% SinvG
	LtSinvL <- t(L) %*% SinvL
	AtSinvG <- t(A) %*% SinvG
	mid.r <- diag(r) + GtSinvG
	mid.L <- solve(Lambda) + LtSinvL
	
	##SGamma and SY##
	S1Gamma <- Gamma.tilde.1 - L.tilde.1 %*% solve( mid.L, LtSinvG )
	S1Z2 <- Z2[1:p1,] / Sigma[1:p1] - L.tilde.1 %*% solve( mid.L, LtSinvZ )

	##Compute H matrices##
	H1 <- LtSinvZ - LtSinvG %*% solve(mid.r, GtSinvZ)           #K x m2
	H2 <- LtSinvL - LtSinvG %*% solve(mid.r, t(LtSinvG))        #K x K
	H3 <- AtSinvG - AtSinvG %*% solve(mid.r, GtSinvG)         #m1 x r
	H4 <- t(GtSinvZ) - t(GtSinvZ) %*% solve(mid.r, GtSinvG)  #m2 x r
	H5 <- LtSinvG - LtSinvG %*% solve(mid.r, GtSinvG)         #K x r	

	##Compute G matrices##
	mid.Lambda <- solve(Lambda) + H2
	G4 <- H4 - t(H1) %*% solve(mid.Lambda, H5)     #m2 x r
	
	##Compute I matrices##
	I1 <- GtSinvG - t(LtSinvG) %*% solve(mid.L, LtSinvG)
	I2 <- GtSinvZ - t(LtSinvG) %*% solve(mid.L, LtSinvZ)
	
	return( -c( -m1 * Gamma.tilde.1 %*% solve(mid.r) + (A[1:p1,] / Sigma[1:p1] - Gamma.tilde.1 %*% solve(mid.r, t(AtSinvG))) %*% H3 - m2 * S1Gamma %*% solve(diag(r) + I1) + (S1Z2 - S1Gamma %*% solve(diag(r) + I1, I2)) %*% G4 ) )
}

Gamma.func.bfgs <- function(par, A, Z2, Gamma.rest, L, Sigma, Lambda, m1, m2, r) {
	p <- nrow(L)
	p1 <- length(par) / r
	
	Gamma <- rbind(matrix(par, nrow=p1, ncol=r), Gamma.rest)

	SinvL <- L / Sigma
	SinvG <- Gamma / Sigma
	SinvA <- A / Sigma
	SinvZ <- Z2 / Sigma
	
	LtSinvZ <- t(SinvL) %*% Z2
	GtSinvZ <- t(SinvG) %*% Z2
	LtSinvG <- t(SinvL) %*% Gamma
	GtSinvG <- t(Gamma) %*% SinvG
	LtSinvL <- t(L) %*% SinvL
	AtSinvG <- t(A) %*% SinvG
	mid.r <- diag(r) + GtSinvG
	
	##Compute H matrices##
	H1 <- LtSinvZ - LtSinvG %*% solve(mid.r, GtSinvZ)           #K x m2
	H2 <- LtSinvL - LtSinvG %*% solve(mid.r, t(LtSinvG))        #K x K
	mid.Lambda <- solve(Lambda) + H2
	
	logdet <- -(m1 + m2) / 2 * log(det(mid.r)) - m2 / 2 * log(det(mid.Lambda))
	tr1.1 <- 1/2 * sum( solve(mid.r, t(AtSinvG)) * t(AtSinvG) )
	tr2.1 <- 1/2 * sum( solve(mid.r, GtSinvZ) * GtSinvZ ) + 1/2 * sum( solve(mid.Lambda, H1) * H1 )
	
	return( (-logdet - tr1.1 - tr2.1)/p )
}







