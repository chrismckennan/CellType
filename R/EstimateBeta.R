##This function estimates beta using GEE
#We assume C.1 is K x n.1, X.i is d x n.i, Omega is K x d
#theta = (beta, l)
#V is n.2*K x n.2*K variance of C.2
#A = - var * X.2' Omega (C.1 P(orthog.X.1) C.1')^(-1) C.1' P(orthog.X.1) Orthog.C
#B.0 = var * X.2' Omega' (C.1 P(orthog.X.1) C.1')^(-1) Omega X.2
#X.tilde = [X.1 Orthog.C.1, X.2], a d x (n.1 - K + n.2) matrix
#z.1 = y.1 Orthog.C.1
#z.2 = y.2 - X.1' Omega' l.hat

Estimator.1 <- function(var, A, B.0, X.tilde, V, l, z.1, z.2) {
	G.l <- Compute.G(var, l, V)
	U <- cbind( rbind( var * diag(length(z.1)), A*var ), rbind(t(A*var), B.0*var + G.l) )
	return( list(beta=solve( X.tilde %*% solve( U, t(X.tilde) ), X.tilde ) %*% solve( U, c(z.1, z.2) ), Var=solve( X.tilde %*% solve(U, t(X.tilde)) )) )
}

Estimator.2 <- function(var, C.1, X.1, X.2, Omega, V, y.1, y.2, theta.0, tol) {
	K <- nrow(C.1)
	n.1 <- ncol(C.1)
	n.2 <- ncol(X.2)
	d <- nrow(X.1)
	
	beta.0 <- theta.0[1:d]
	l.0 <- theta.0[(d+1):(d+K)]
	
	X.tilde.1 <- rbind(X.1, C.1)
	X.tilde.2 <- rbind(X.2, Omega %*% X.2)
	T.0.1 <- 1/var * X.tilde.1 %*% t(X.tilde.1)
	
	diff <- 1 + tol
	count <- 1
	while(diff > tol && count < 1e3) {
		G.0 <- Compute.G(var, l.0, V)
		stand.resid.1.0 <- 1/var * (y.1 - t(X.tilde.1) %*% theta.0)
		stand.resid.2.0 <- solve(G.0, y.2 - t(X.tilde.2) %*% theta.0)
		
		U.0 <- X.tilde.1 %*% stand.resid.1.0 + X.tilde.2 %*% stand.resid.2.0
		diff <- max(abs(U.0))
		T.0 <- T.0.1 + X.tilde.2 %*% solve(G.0, t(X.tilde.2))
		
		theta.0 <- theta.0 + solve(T.0, U.0)
		beta.0 <- theta.0[1:d]
		l.0 <- theta.0[(d+1):(d+K)]
		
		count <- count + 1
	}
	return(list(beta=beta.0, l=l.0, iter=count, grad=diff, T=T.0, naive.T=T.0.1)) 
}



Compute.G <- function(var, l, V) {
	K <- length(l)
	n.2 <- nrow(V)/K
	G <- array(0, dim=c(n.2, n.2))
	for (r in 1:n.2) {
		ind.r <- ((r-1)*K + 1):(r*K)
		for (c in r:n.2) {
			ind.c <- ((c-1)*K + 1):(c*K)
			G[r,c] <- sum( l * V[ind.r, ind.c] %*% l )
			if (c != r) {
				G[c,r] <- G[r,c]
			}
		}
	}
	return(G + var * diag(n.2))
}


#SYY1, SYY2 are p x p
#SY1F is p x K
#SFF is K x K
#L.0 is p x K
#Sigma.0 is a p-vector
ExpMax <- function(L.0, Sigma.0, n.1, n.2, d, SYY1, SYY2, SY1F, SFF, tol) {
	K <- nrow(SFF)
	diff <- 1 + tol
	count <- 1
	while (diff > tol && count < 1e4) {
		H.0 <- solve(L.0 %*% t(L.0) + Sigma.0 * diag(1))
		
		##Update L##
		L.1 <- t( solve( (n.2 - d)*(diag(K) - t(L.0) %*% H.0 %*% L.0) + t(L.0) %*% H.0 %*% SYY2 %*% H.0 %*% L.0 + SFF, t( SYY2 %*% H.0 %*% L.0 + SY1F ) ) )
		H.10 <- solve( L.1 %*% t(L.1) + Sigma.0 * diag(1) )
		
		##Update Sigma##
		Sigma.1 <- as.numeric(2/(n.1+n.2-2*d) * ( 1/2 * diag(SYY2 + SYY1) + 1/2 * diag( L.1 %*% ( (n.2-d) * (diag(K) - t(L.1) %*% H.10 %*% L.1) + t(L.1) %*% H.10 %*% SYY2 %*% H.10 %*% L.1 + SFF ) %*% t(L.1) ) - diag( L.1 %*% ( t(L.1) %*% H.10 %*% SYY2 + t(SY1F) ) ) ))
		
		##Update Parameters##
		diff <- sqrt( norm(t(L.0) - t(L.1), type = "F")^2 + (Sigma.0 - Sigma.1)^2 )
		L.0 <- L.1
		Sigma.0 <- Sigma.1
		count <- count + 1
	}
	return(list(L=L.0, Sigma=Sigma.0, n.iter=count, norm.diff=diff))
}

ExpMax.multL <- function(L.0, Sigma.0, n.1, n.2, d, SYY1, SYY2, SY1F, SFF, tol) {
	K <- nrow(SFF)
	p <- nrow(L.0)
	diff <- 1 + tol
	count <- 1
	while (diff > tol && count < 1e4) {
		Sigma.invL <- L.0
		for (g in 1:p) {
			Sigma.invL[g,] <- L.0[g,]/Sigma.0[g]
		}
		H.0 <- diag(1/Sigma.0) - Sigma.invL %*% solve(diag(K) + t(L.0) %*% Sigma.invL) %*% t(Sigma.invL)
		
		##Update L##
		L.1 <- t( solve( (n.2 - d)*(diag(K) - t(L.0) %*% H.0 %*% L.0) + t(L.0) %*% H.0 %*% SYY2 %*% H.0 %*% L.0 + SFF, t( SYY2 %*% H.0 %*% L.0 + SY1F ) ) )
		Sigma.invL <- L.1
		for (g in 1:p) {
			Sigma.invL[g,] <- L.1[g,]/Sigma.0[g]
		}
		H.10 <- diag(1/Sigma.0) - Sigma.invL %*% solve(diag(K) + t(L.1) %*% Sigma.invL) %*% t(Sigma.invL)
		
		##Update Sigma##
		Sigma.1 <- as.vector(2/(n.1+n.2-2*d) * ( 1/2 * diag(SYY2 + SYY1) + 1/2 * diag( L.1 %*% ( (n.2-d) * (diag(K) - t(L.1) %*% H.10 %*% L.1) + t(L.1) %*% H.10 %*% SYY2 %*% H.10 %*% L.1 + SFF ) %*% t(L.1) ) - diag( L.1 %*% ( t(L.1) %*% H.10 %*% SYY2 + t(SY1F) ) ) ))
		
		##Update Parameters##
		diff <- sqrt( norm(t(L.0) - t(L.1), type = "F")^2 + sum((Sigma.0 - Sigma.1)^2) )
		L.0 <- L.1
		Sigma.0 <- Sigma.1
		count <- count + 1
	}
	return(list(L=L.0, Sigma=Sigma.0, n.iter=count, norm.diff=diff))
}


##The above code works##
#n.sim <- 100000
#p <- 100
#sim.effect <- matrix(rnorm(3*p), nrow=3, ncol=p)
#var.sim <- 1
#sim.effect <- t(sim.effect)
#sim.y <- matrix(rnorm(3*n.sim), nrow=p, ncol=n.sim) + sim.effect %*% matrix(rnorm(3*n.sim), nrow=3, ncol=n.sim)
#tmp.y <- sim.y[1,]
#exp.max <- ExpMax(rbind(sim.effect[1,]), rbind(c(1)), 0, n.sim, 0, as.matrix(0), as.matrix(sum(tmp.y*tmp.y)), rbind(c(0,0,0)), array(0, dim=c(3,3)), 1e-6)



####The below code will perform EM on the FULL model for the confounders + partially observed cell types####

#L is p x k
#Sigma is a p-vector
#Gamma is p x r
#Lambda is k x k
#Y1 is p x (n.1 - d)
#Y2 is p x (n.2 - d)
#F.mat is k x (n.1 - d)
FullEM <- function(L.0, Sigma.0, Gamma.0, Lambda.0, Y1, Y2, F.mat, Update.L, Update.Lambda, Mon.Like, tol, max.iter) {
	n.1 <- ncol(Y1)   #This is actually n.1 - d
	n.2 <- ncol(Y2)  #This is actually n.2 - d
	K <- nrow(F.mat)
	r <- ncol(Gamma.0)
	p <- nrow(Gamma.0)
	SFF <- F.mat %*% t(F.mat)    #A k x k matrix
	Y1Ft <- Y1 %*% t(F.mat)      #A p x K matrix
	diff <- tol + 1
	count <- 1
	if (Mon.Like) {
		l.vec <- rep(NA, max.iter)
	}
	
	while(diff > tol && count < max.iter) {
		###Update for Lambda###
		S.chol.0 <- t(chol(Lambda.0))     #SS' = Lambda.0
		
		SinvL.0 <- sweep(L.0, 1, 1/Sigma.0, "*")    #Sigma.0^{-1} L.0
		SinvG.0 <- sweep(Gamma.0, 1, 1/Sigma.0, "*")    #Sigma.0^{-1} Gamma.0
		
		LtSinvL.0 <- t(L.0) %*% SinvL.0      #L.0' Sigma.0^{-1} L.0, a k x k matrix
		LtSinvG.0 <- t(SinvL.0) %*% Gamma.0    #L' Sigma.0^{-1} Gamma.0
		GtSinvG.0 <- t(SinvG.0) %*% Gamma.0  #Gamma.0' Sigma.0^{-1} Gamma.0
		
		LtSinvY2.0 <- t(SinvL.0) %*% Y2     #L.0' Sigma.0^{-1} Y2
		GtSinvY2.0 <- t(SinvG.0) %*% Y2     #Gamma.0' Sigma.0^{-1} Y2
		
		LtHL.0 <- ( LtSinvL.0 - LtSinvL.0 %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvL.0 ) ) - ( LtSinvG.0 - LtSinvL.0 %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ) ) %*% solve( diag(r) + GtSinvG.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ), t( LtSinvG.0 - LtSinvL.0 %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ) ) )         #L.0' H L.0, a K x K matrix
		
		LtHY2.0 <- ( LtSinvY2.0 - LtSinvL.0 %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvY2.0 ) ) - ( LtSinvG.0 - LtSinvL.0 %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ) ) %*% solve( diag(r) + GtSinvG.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ), GtSinvY2.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvY2.0 ) )        #L.0' H Y2, a K x n.2 matrix

		if (Update.Lambda) {
			Lambda.1 <- ( Lambda.0 - Lambda.0 %*% LtHL.0 %*% Lambda.0 ) + 1/n.2 * ( Lambda.0 %*% LtHY2.0 %*% t(LtHY2.0) %*% Lambda.0 )
			diff.Lambda <- norm(Lambda.1 - Lambda.0, type="F")
			Lambda.0 <- Lambda.1
			S.chol.0 <- t(chol(Lambda.0))     #SS' = Lambda.0. This will not change for the rest of the iteration
		}
		
		###Update L###
		#I DO NOT need to update LtSinvL.0, LtSinvG.0, GtSinvG.0, LtSinvY2.0, GtSinvY2.0
		
		Y1mLF.0 <- Y1 - L.0 %*% F.mat      #Y1 - L.0 F, a p x n.1 matrix. This will also be referred to as 'Z' in the below code
		
		GtSinvZ.0 <- t(SinvG.0) %*% Y1mLF.0        #Gamma.0' Sigma^{-1} (Y1 - L.0 F)
		
		Middle.0 <- diag(r) + GtSinvG.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 )       #See notes
		
		LtHL.0 <- ( LtSinvL.0 - LtSinvL.0 %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvL.0 ) ) - ( LtSinvG.0 - LtSinvL.0 %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ) ) %*% solve( Middle.0, t( LtSinvG.0 - LtSinvL.0 %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ) ) )         #L.0' H L.0, a K x K matrix
		
		LtHY2.0 <- ( LtSinvY2.0 - LtSinvL.0 %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvY2.0 ) ) - ( LtSinvG.0 - LtSinvL.0 %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ) ) %*% solve( Middle.0, GtSinvY2.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvY2.0 ) )        #L.0' H Y2, a K x n.2 matrix
		
		GtGZ.0 <- GtSinvZ.0 - GtSinvG.0 %*% solve( diag(r) + GtSinvG.0, GtSinvZ.0 )     #Gamma.0' G (Y1 - L.0 F), a r x n.1 matrix
		
		GtHY2.0 <- ( GtSinvY2.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvY2.0 ) ) - ( GtSinvG.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ) ) %*% solve( Middle.0, GtSinvY2.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvY2.0 ) )        #Gamma.0' H Y2, a r x n.2 matrix
		
		GtHL.0 <- ( t(LtSinvG.0) - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvL.0 ) ) - ( GtSinvG.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ) ) %*% solve( Middle.0, t(LtSinvG.0) - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvL.0 ) )        #Gamma.0' H L, a r x K matrix
		
		if (Update.L) {
			U.0 <- 	SFF + n.2 * (Lambda.0 - Lambda.0 %*% LtHL.0 %*% Lambda.0) + Lambda.0 %*% LtHY2.0 %*% t(LtHY2.0) %*% Lambda.0         #A K x K matrix
			J.0 <- Y1Ft - Gamma.0 %*% GtGZ.0 %*% t(F.mat) + Y2 %*% t(LtHY2.0) %*% Lambda.0 + Gamma.0 %*% ( n.2 * GtHL.0 %*% Lambda.0 - GtHY2.0 %*% t(LtHY2.0) %*% Lambda.0 )       #A p x K matrix
			L.1 <- J.0 %*% solve(U.0)
			diff.L <- max(abs(L.0 - L.1))
			L.0 <- L.1		
		}
		
		
		###Update Gamma###
		#I DO NOT need to update GtSinvG.0, GtSinvY2.0
		#I NEED to update SinvL.0, LtSinvL.0, LtSinvG.0, LtSinvY2.0, Y1mLF.0, GtSinvZ.0
		
		SinvL.0 <- sweep(L.0, 1, 1/Sigma.0, "*")    #Sigma.0^{-1} L.0
		
		LtSinvL.0 <- t(L.0) %*% SinvL.0      #L.0' Sigma.0^{-1} L.0, a k x k matrix
		LtSinvG.0 <- t(SinvL.0) %*% Gamma.0    #L' Sigma.0^{-1} Gamma.0
		LtSinvY2.0 <- t(SinvL.0) %*% Y2     #L.0' Sigma.0^{-1} Y2
		
		Y1mLF.0 <- Y1 - L.0 %*% F.mat      #Y1 - L.0 F, a p x n.1 matrix. This will also be referred to as 'Z' in the below code
		GtSinvZ.0 <- t(SinvG.0) %*% Y1mLF.0        #Gamma.0' Sigma^{-1} (Y1 - L.0 F)
		
		Middle.0 <- diag(r) + GtSinvG.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 )       #See notes
		
		LtHY2.0 <- ( LtSinvY2.0 - LtSinvL.0 %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvY2.0 ) ) - ( LtSinvG.0 - LtSinvL.0 %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ) ) %*% solve( Middle.0, GtSinvY2.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvY2.0 ) )        #L.0' H Y2, a K x n.2 matrix
		
		GtGZ.0 <- GtSinvZ.0 - GtSinvG.0 %*% solve( diag(r) + GtSinvG.0, GtSinvZ.0 )     #Gamma.0' G (Y1 - L.0 F), a r x n.1 matrix	
		
		GtHY2.0 <- ( GtSinvY2.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvY2.0 ) ) - ( GtSinvG.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ) ) %*% solve( Middle.0, GtSinvY2.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvY2.0 ) )        #Gamma.0' H Y2, a r x n.2 matrix
		
		GtHL.0 <- ( t(LtSinvG.0) - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvL.0 ) ) - ( GtSinvG.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ) ) %*% solve( Middle.0, t(LtSinvG.0) - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvL.0 ) )        #Gamma.0' H L, a r x K matrix
		
		GtGG.0 <- GtSinvG.0 - GtSinvG.0 %*% solve( diag(r) + GtSinvG.0, GtSinvG.0 )          #Gamma.0' G Gamma.0
		
		A.GtHG.0 <- GtSinvG.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve(diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0)       #See notes
		GtHG.0 <- A.GtHG.0 - A.GtHG.0 %*% solve(Middle.0, A.GtHG.0)       #Gamma.0' H Gamma.0
		
		V.0 <- (n.1 + n.2) * diag(r) - n.1 * GtGG.0 - n.2 * GtHG.0 + GtGZ.0 %*% t(GtGZ.0) + GtHY2.0 %*% t(GtHY2.0)        #An r x r matrix
		M.0 <- Y1mLF.0 %*% t(GtGZ.0) + Y2 %*% t(GtHY2.0) - L.0 %*% ( Lambda.0 %*% LtHY2.0 %*% t(GtHY2.0) - n.2 * Lambda.0 %*% t(GtHL.0) )
		Gamma.1 <- M.0 %*% solve(V.0)
		
		#Make Gamma.1 and Sigma.0 identifiable#
		#v.rotate.1 <- svd( t(sweep(Gamma.1, 1, 1/Sigma.0, "*")) %*% Gamma.1 / p )$v
		#Gamma.1 <- Gamma.1 %*% v.rotate.1      #This rotation is such that 1/p * Gamma' Sigma^{-1} Gamma is diagonal with decreasing values along the diagonal (see Hastie paper)
		
		##Compute Log likelihood##
		if (Mon.Like) {
			tmp.Y2.0 <- t(GtSinvY2.0) - t(LtSinvY2.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 )
			l.vec[count] <- -n.1/2 * determinant(GtSinvG.0 + diag(r), logarithm = TRUE)$modulus[1] + 1/2 * sum( rowSums( t(GtSinvZ.0) * t( solve( diag(r) + GtSinvG.0, GtSinvZ.0 ) ) ) ) - n.2/2 * determinant(Middle.0, logarithm = TRUE)$modulus[1] + 1/2 * sum( rowSums( tmp.Y2.0 * t( solve(Middle.0, t(tmp.Y2.0)) ) ) )
		}
		
		diff.Gamma <- max(abs(Gamma.0 - Gamma.1))
		Gamma.0 <- Gamma.1
		
		
		##Update Sigma.0##
		#I need to update everything with a Gamma.0 in it
		
		SinvG.0 <- sweep(Gamma.0, 1, 1/Sigma.0, "*")    #Sigma.0^{-1} Gamma.0
		LtSinvG.0 <- t(SinvL.0) %*% Gamma.0    #L' Sigma.0^{-1} Gamma.0
		GtSinvG.0 <- t(SinvG.0) %*% Gamma.0  #Gamma.0' Sigma.0^{-1} Gamma.0
		GtSinvZ.0 <- t(SinvG.0) %*% Y1mLF.0        #Gamma.0' Sigma^{-1} (Y1 - L.0 F)
		
		Middle.0 <- diag(r) + GtSinvG.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 )       #See notes
		
		A.GtHG.0 <- GtSinvG.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve(diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0)       #See notes
		GtHG.0 <- A.GtHG.0 - A.GtHG.0 %*% solve(Middle.0, A.GtHG.0)       #Gamma.0' H Gamma.0
		
		GtGG.0 <- GtSinvG.0 - GtSinvG.0 %*% solve( diag(r) + GtSinvG.0, GtSinvG.0 )          #Gamma.0' G Gamma.0
		
		GtGZ.0 <- GtSinvZ.0 - GtSinvG.0 %*% solve( diag(r) + GtSinvG.0, GtSinvZ.0 )     #Gamma.0' G (Y1 - L.0 F), a r x n.1 matrix
		
		GtHY2.0 <- ( GtSinvY2.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvY2.0 ) ) - ( GtSinvG.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ) ) %*% solve( Middle.0, GtSinvY2.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvY2.0 ) )        #Gamma.0' H Y2, a r x n.2 matrix
		
		LtHY2.0 <- ( LtSinvY2.0 - LtSinvL.0 %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvY2.0 ) ) - ( LtSinvG.0 - LtSinvL.0 %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ) ) %*% solve( Middle.0, GtSinvY2.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvY2.0 ) )        #L.0' H Y2, a K x n.2 matrix
		
		GtHL.0 <- ( t(LtSinvG.0) - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvL.0 ) ) - ( GtSinvG.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ) ) %*% solve( Middle.0, t(LtSinvG.0) - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvL.0 ) )        #Gamma.0' H L, a r x K matrix
		
		LtHL.0 <- ( LtSinvL.0 - LtSinvL.0 %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvL.0 ) ) - ( LtSinvG.0 - LtSinvL.0 %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ) ) %*% solve( Middle.0, t( LtSinvG.0 - LtSinvL.0 %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ) ) )         #L.0' H L.0, a K x K matrix
		
		C.0 <- Gamma.0 %*% t( (n.1 + n.2) * diag(r) - n.1 * GtGG.0 - n.2 * GtHG.0 + GtGZ.0 %*% t(GtGZ.0) + GtHY2.0 %*% t(GtHY2.0) ) + 2 * L.0 %*% t( GtHY2.0 %*% t(LtHY2.0) %*% Lambda.0 - n.2 * GtHL.0 %*% Lambda.0 ) - 2 * Y1mLF.0 %*% t(GtGZ.0) - 2 * Y2 %*% t(GtHY2.0)      #To be multiplied by Gamma.0 when updating Sigma.0. This is a p x r matrix
		
		D.0 <- L.0 %*% t( n.2 * (Lambda.0 - Lambda.0 %*% LtHL.0 %*% Lambda.0) + Lambda.0 %*% LtHY2.0 %*% t(LtHY2.0) %*% Lambda.0 ) - 2 * Y2 %*% t(Lambda.0 %*% LtHY2.0)      #This is a p x K matrix
		
		##Compute Log likelihood##
		if (Mon.Like) {
			tmp.Y2.0 <- t(GtSinvY2.0) - t(LtSinvY2.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 )
			l.vec[count] <- -l.vec[count] + ( -n.1/2 * determinant(GtSinvG.0 + diag(r), logarithm = TRUE)$modulus[1] + 1/2 * sum( rowSums( t(GtSinvZ.0) * t( solve( diag(r) + GtSinvG.0, GtSinvZ.0 ) ) ) ) - n.2/2 * determinant(Middle.0, logarithm = TRUE)$modulus[1] + 1/2 * sum( rowSums( tmp.Y2.0 * t( solve(Middle.0, t(tmp.Y2.0)) ) ) ) )
		}
		
		Sigma.1 <- 1/(n.1 + n.2) * ( rowSums(Y2 * Y2) + rowSums(Y1mLF.0 * Y1mLF.0) + rowSums(Gamma.0 * C.0) + rowSums(L.0 * D.0) )
		
		diff.Sigma <- max(abs(Sigma.1 - Sigma.0))
		Sigma.0 <- Sigma.1
		
		diff <- diff.Lambda + diff.L + diff.Gamma + diff.Sigma
		count <- count + 1
	}
	return(list( L=L.0, Gamma=Gamma.0, Sigma=Sigma.0, Lambda=Lambda.0, n.iter=count, conv=c(diff.Lambda, diff.L, diff.Sigma, diff.Gamma) ))
}



SingleFactor <- function(Y, Sigma, r, tol, max.iter) {    #It's assumed Sigma is fixed in this code
	Gamma.0 <- fa.pc(t(Y), r)$Gamma + sqrt(0.5) * matrix(rnorm(nrow(Y) * r), nrow=nrow(Y), ncol=r)
	nm1 <- ncol(Y) - 1
	p <- nrow(Y)
	r <- ncol(Gamma.0)
	diff <- tol + 1
	count <- 1

	SinvG.0 <- sweep(Gamma.0, 1, 1/Sigma, "*")
	GtSinvG.0 <- t(Gamma.0) %*% SinvG.0
	YtSinvG.0 <- t(Y) %*% SinvG.0
	GtGY.0 <- t(YtSinvG.0) - GtSinvG.0 %*% solve(diag(r) + GtSinvG.0, t(YtSinvG.0))
	GtGG.0 <- GtSinvG.0 - GtSinvG.0 %*% solve( diag(r) + GtSinvG.0, GtSinvG.0 )
	
	l.vec <- rep(NA, max.iter)

	while(diff > tol && count < max.iter) {
		M.0 <- Y %*% t(GtGY.0)
		U.0 <- nm1 * (diag(r) - GtGG.0) + GtGY.0 %*% t(GtGY.0)
		Gamma.1 <- t( solve(t(U.0), t(M.0)) )
		
		SinvG.1 <- sweep(Gamma.1, 1, 1/Sigma, "*")
		GtSinvG.1 <- t(Gamma.1) %*% SinvG.1
		YtSinvG.1 <- t(Y) %*% SinvG.1
		GtGY.1 <- t(YtSinvG.1) - GtSinvG.1 %*% solve(diag(r) + GtSinvG.1, t(YtSinvG.1))
		GtGG.1 <- GtSinvG.1 - GtSinvG.1 %*% solve( diag(r) + GtSinvG.1, GtSinvG.1 )		
		
		ll.0 <- -1/2/p * determinant(GtSinvG.0 + diag(r), logarithm = TRUE)$modulus[1] - 1/2/p/nm1 * sum( rowSums(YtSinvG.0 * t( solve(diag(r) + GtSinvG.0, t(YtSinvG.0)) )) )
		ll.1 <- -1/2/p * determinant(GtSinvG.1 + diag(r), logarithm = TRUE)$modulus[1] - 1/2/p/nm1 * sum( rowSums(YtSinvG.1 * t( solve(diag(r) + GtSinvG.1, t(YtSinvG.1)) )) )
		
		diff <- abs(ll.0 - ll.1)/abs(ll.0)
		l.vec[count] <- diff
		
		Gamma.0 <- Gamma.1
		SinvG.0 <- SinvG.1
		GtSinvG.0 <- GtSinvG.1
		YtSinvG.0 <- YtSinvG.1
		GtGY.0 <- GtGY.1
		GtGG.0 <- GtGG.1
		
		count <- count + 1
	}
}



#cate.em <- function(Y, r)
#  Y <- t(Y)
#	p <- ncol(Y)
#    n <- nrow(Y)
#    init <- fa.pc(Y, r)
#    Gamma <- init$Gamma  + sqrt(0.5) * matrix(rnorm(ncol(Y) * r), nrow=ncol(Y), ncol=r)
#    invSigma <- 1/init$Sigma
#    llh <- -Inf
#    I <- diag(rep(1, r))
#    sample.var <- colMeans(Y^2)
#    tilde.Gamma <- sqrt(invSigma) * Gamma
#    M <- diag(r) + t(tilde.Gamma) %*% tilde.Gamma
#    eigenM <- eigen(M, symmetric = TRUE)
#    YSG <- Y %*% (invSigma * Gamma)
#    logdetY <- -sum(log(invSigma)) + sum(log(eigenM$values))
#    B <- 1/sqrt(eigenM$values) * t(eigenM$vectors) %*% t(YSG)
#    logtrY <- sum(invSigma * sample.var) - sum(B^2)/n
#    llh <- -logdetY - logtrY
#    converged <- FALSE
#    for (iter in 1:maxiter) {
#        varZ <- eigenM$vectors %*% (1/eigenM$values * t(eigenM$vectors))
#        EZ <- YSG %*% varZ
#        EZZ <- n * varZ + t(EZ) %*% EZ
#        eigenEZZ <- eigen(EZZ, symmetric = TRUE)
#        YEZ <- t(Y) %*% EZ
#        G <- sqrt(eigenEZZ$values) * t(eigenEZZ$vectors)
#        invSigma <- 1/(sample.var - 2/n * rowSums(YEZ * Gamma) + 
#            1/n * rowSums((Gamma %*% t(G))^2))
#        Gamma <- YEZ %*% eigenEZZ$vectors %*% (1/eigenEZZ$values * 
#            t(eigenEZZ$vectors))
#        tilde.Gamma <- sqrt(invSigma) * Gamma
#        M <- diag(r) + t(tilde.Gamma) %*% tilde.Gamma
#        eigenM <- eigen(M, T)
#        YSG <- Y %*% (invSigma * Gamma)
#        old.llh <- llh
#        logdetY <- -sum(log(invSigma)) + sum(log(eigenM$values))
#        B <- 1/sqrt(eigenM$values) * t(eigenM$vectors) %*% t(YSG)
#        logtrY <- sum(invSigma * sample.var) - sum(B^2)/n
#        llh <- -logdetY - logtrY
#        if (abs(llh - old.llh) < tol * abs(llh)) {
#            converged <- TRUE
#            break
#        }
#    }

#My EM for a single factor Gamma matches CATE's updates...It looks like the problem is the starting value for Gamma



##Penalized Huber's Loss proof of principle##
#Assume that d = 1, meaning Y is a p x 1 vector
#Omega is K x d and Gamma is p x K
#k.hub is the point where the second derivative of Huber's loss function is discontinuous
#The penalty is q/2 * (alpha - Omega)' inv.V.Omega (apha - Omega)
#The updates are weighted least squares

Pen.Huber <- function(Y, Sigma, Gamma, Omega, inv.V.Omega, k.hub = 1.345, ind.use=NULL, q = 0, max.iter=1e4, tol=1e-8) {
	K <- ncol(Gamma)
	if (q == 0) {
	  inv.V.Omega <- matrix(0, nrow=K, ncol=K)
	  Omega <- cbind(rep(0, length=K))
	}
	if (! is.null(ind.use)) {
		Y <- Y[ind.use]
		Sigma <- Sigma[ind.use]
		if (K > 1) {
			Gamma <- Gamma[ind.use,]
		} else {
			Gamma <- cbind(Gamma[ind.use])
		}
	}
	Gamma <- Gamma / Sigma
	Y <- Y / Sigma - Gamma %*% Omega
	
	alpha.tilde.0 <- solve( t(Gamma) %*% Gamma + q * inv.V.Omega, t(Gamma) %*% Y)       #Penalized least squares as a starting point

	for (iter in 1:max.iter) {
		psi.0 <- psi.huber(Y - Gamma %*% alpha.tilde.0, k.hub)
		if (K == 1) {
			grad.0 <- q * inv.V.Omega %*% alpha.tilde.0 - sum(Gamma * psi.0)
			diff <- abs(grad.0)
		} else {
			grad.0 <- q * inv.V.Omega %*% alpha.tilde.0 - apply(Gamma * psi.0, 2, sum)
			diff <- max(abs(grad.0))
		}
		if (diff < tol) {
			return(list(alpha=alpha.tilde.0 + Omega, n.iter=iter, out=1, resid=Y - Gamma %*% alpha.tilde.0))
		}
		
		weights.0 <- Weights.huber(Y - Gamma %*% alpha.tilde.0, k.hub)
		tmp.mat <- Gamma * weights.0
		alpha.tilde.0 <- solve( t(tmp.mat) %*% Gamma + q * inv.V.Omega, t(tmp.mat) %*% Y )
		if (K == 1) {
			alpha.tilde.0 <- matrix(alpha.tilde.0, nrow=1, ncol=1)
		}
		
	}
	return(list(alpha=alpha.0, n.iter=iter, out=0))
}



psi.huber <- function(x, k) {
	ind.small <- as.numeric(abs(x) <= k)
	return( ind.small * x/k + (1 - ind.small) * sign.vec(x) )
}

Weights.huber <- function(x, k) {
	ind.small <- as.numeric(abs(x) <= k)
	return( ind.small * 1/k + (1-ind.small) / abs(x) )
}

sign.vec <- function(x) {
	return( as.numeric( x > 0 ) - as.numeric( x < 0 ) )
}

rho.huber <- function(x, k) {
	ind.small <- abs(x) <= k
	return( sum(1/2/k * x[ind.small]^2) + sum(abs(x[!ind.small]) - k/2) )
}

fun.value <- function(x, k, diff, inv.V, q) {
	return( (rho.huber(x, k) + q/2 * t(diff) %*% inv.V %*% diff)/p )
}
 
