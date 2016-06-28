##These functions are to be used as input into TurboEM##
#We need to write EM updates for L with (Sigma, Gamma) fixed and updates for Gamma with (Sigma, L) fixed

###Run optimization for a given simulation###
#This assumes Gamma is 0

Correct.CellType <- function(M, X, C, ind.1, K.use, cov.interest=c(1), B.sim) {     #The covariate of interest is in B.sim, NOT X
	n <- nrow(X)
	K <- ncol(C)
	d <- ncol(X)
	p <- nrow(M)
	Gamma <- cbind(rep(0,p))
	n.1 <- length(ind.1)
	n.2 <- n - n.1
	m.1 <- n.1 - d
	m.2 <- n.2 - d
	
	X.1 <- X[ind.1,]
	X.2 <- X[-ind.1,]
	C.1 <- C[ind.1,]
	
	orthog.X.1 <- diag(n.1) - X.1 %*% solve(t(X.1) %*% X.1, t(X.1))
	orthog.X.2 <- diag(n.2) - X.2 %*% solve(t(X.2) %*% X.2, t(X.2))
	Lambda <- 1/(n.1 - d) * t(C.1) %*% orthog.X.1 %*% C.1  #K x K
	
	Q.X2.orthog <- qr.Q(qr(X.2), complete=T)[,(d+1):n.2]     #n.2 x (n.2 - d)
	
	F.mat <- t(C.1) %*% orthog.X.1      #K x n.1
	Z1 <- M[,ind.1] %*% orthog.X.1
	Z2 <- M[,-ind.1] %*% orthog.X.2
	
	L.0 <- Z1 %*% t(F.mat) %*% solve(F.mat %*% t(F.mat))
	Sigma <- fa.em(t( M %*% qr.Q(qr(X), complete=T)[,(d+1):n] ), r=K.use)$Sigma
	out.turbo <- turboem(par=c(L.0), fixptfn=Update.L, objfn=ll.L, method="squarem", Z1=Z1, Z2=Z2, F.mat=F.mat, Lambda=Lambda, Gamma=Gamma, Sigma=Sigma, m2=m.2, control.run=list(tol=1e-8, convtype="objfn"))
	L <- matrix(out.turbo$pars, nrow=p, ncol=K)
	
	C.2.hat <- solve( t(L) %*% (L / Sigma), t(L) %*% (Z2 / Sigma) ) + t(C.1) %*% X.1 %*% solve(t(X.1) %*% X.1) %*% t(X.2)
	#rm(L)
	cov.total <- cbind(rbind(X.1, X.2), rbind(cbind(C.1), t(C.2.hat)))   #n x (d+K)
	dof <- n - d - K
	M.use <- cbind(M[,ind.1], M[,-ind.1])
	beta.all <- M.use %*% cov.total %*% solve(t(cov.total) %*% cov.total)
	var.beta <- solve(t(cov.total) %*% cov.total)
	orthog.total <- diag(n) - cov.total %*% solve(t(cov.total) %*% cov.total) %*% t(cov.total)
	Sigma <- rowSums((M.use %*% orthog.total) * M.use) / dof
	rm(M.use)
	
	p.values <- cbind(2 - 2 * pt( abs(beta.all[,cov.interest[1]+1]) / sqrt(Sigma) / sqrt(var.beta[cov.interest[1]+1, cov.interest[1]+1]), df=dof ))
	if (length(cov.interest) > 1) {
		for (i in cov.interest[-1]) {
			p.values <- cbind(p.values, 2 - 2 * pt( beta.all[,i+1] / sqrt(Sigma) / sqrt(var.beta[i+1, i+1]), df=dof ))
		}
	}
	
	q.values <- qvalue(p.values[,1])
	fsr.i <- false.sign.results(B.sim[,cov.interest[1]], beta.all[,cov.interest[1]+1], q.values$qvalues)
	return(list(EM=out.turbo, beta=beta.all, p.values=p.values, q.value=q.values, fsr=fsr.i))
}

#The above code works


#####Update L with (Sigma, Gamma) Fixed#####

#turbo.L <- turboem(par=c(tmp.L), fixptfn=Update.L, objfn=ll.L, method="squarem", Z1=Z1.sim, Z2=Z2.sim, F.mat=F.mat, Lambda=Lambda.0, Gamma=Gamma, Sigma=Sigma.use, m2=m2.sim, control.run=list(tol=1e-10, convtype="objfn"))

Update.L <- function(par, Z1, Z2, F.mat, Lambda, Gamma, Sigma, m2) {
	p <- nrow(Z1)
	K <- nrow(F.mat)
	r <- ncol(Gamma)
	
	L <- matrix(par, nrow=p, ncol=K)
	
	##Compute Preliminary Matrices##
	SinvL <- L / Sigma
	SinvG <- Gamma / Sigma
	SinvZ <- Z2 / Sigma
	
	LtSinvZ <- t(SinvL) %*% Z2
	GtSinvZ <- t(SinvG) %*% Z2
	LtSinvG <- t(SinvL) %*% Gamma
	GtSinvG <- t(Gamma) %*% SinvG
	LtSinvL <- t(L) %*% SinvL
	mid.r <- diag(1, nrow=r, ncol=r) + GtSinvG
	
	##Compute H matrices##
	H1 <- LtSinvZ - LtSinvG %*% solve(mid.r, GtSinvZ)           #K x m2
	H2 <- LtSinvL - LtSinvG %*% solve(mid.r, t(LtSinvG))        #K x K
	
	##Compute G matrices##
	mid.Lambda <- solve(Lambda) + H2
	G1 <- H1 - H2 %*% solve(mid.Lambda, H1)         #K x m2
	G2 <- H2 - H2 %*% solve(mid.Lambda, H2)         #K x K

	#####   Update L   #####
	M <- Z1 %*% t(F.mat) + Z2 %*% t(G1) %*% Lambda
	U <- F.mat %*% t(F.mat) + Lambda %*% G1 %*% t(G1) %*% Lambda + m2 * (Lambda - Lambda %*% G2 %*% Lambda)
	return(c(M %*% solve(U)))     #EM update for L
}

ll.L <- function(par, Z1, Z2, F.mat, Lambda, Gamma, Sigma, m2) {
	p <- nrow(Z1)
	K <- nrow(F.mat)
	L <- matrix(par, nrow=p, ncol=K)
	r <- ncol(Gamma)
	
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
	mid.r <- diag(1, nrow=r, ncol=r) + GtSinvG
	
	##Compute H matrices##
	H1 <- LtSinvZ - LtSinvG %*% solve(mid.r, GtSinvZ)           #K x m2
	H2 <- LtSinvL - LtSinvG %*% solve(mid.r, t(LtSinvG))        #K x K
	mid.Lambda <- solve(Lambda) + H2
	
	logdet <- -m2/2 * log(det(mid.Lambda))
	tr1 <- -1/2 * sum( SinvA * A ) + 1/2 * sum( solve(mid.r, t(AtSinvG)) * t(AtSinvG) )
	tr2 <- 1/2 * sum( solve(mid.Lambda, H1) * H1 )
	return((logdet + tr1 + tr2)/p)
}

convfn.user.objfn <- function(old, new) {
	return( abs(new - old)/abs(old) < tol )
}


#####Update Gamma with (Sigma, L) Fixed#####

#turbo.Gamma <- turboem(par=c(Gamma), fixptfn=Update.Gamma, objfn=ll.Gamma, method="squarem", A=Z1.sim - L %*% F.mat, Z2=Z2.sim, Lambda=Lambda.0, L=L, Sigma=Sigma.use, m1=m1.sim, m2=m2.sim, control.run=list(tol=1e-9, convtype="objfn"))

Update.Gamma <- function(par, A, Z2, Lambda, L, Sigma, m1, m2) {    #The update for Gamma depends on Z1 and F.mat ONLY through A, which is known
	p <- nrow(Z1)
	r <- length(par) / p
	Gamma <- matrix(par, nrow=p, ncol=r)

	##Compute Preliminary Matrices##
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
	G3 <- H3															      #m1 x r
	G4 <- H4 - t(H1) %*% solve(mid.Lambda, H5)     #m2 x r
	G5 <- H5 - H2 %*% solve(mid.Lambda, H5)        #K x r
	G6 <- H6																  #r x r
	G7 <- H6 - t(H5) %*% solve(mid.Lambda, H5)     #r x r
	
	#Update Gamma#
	M <- A %*% G3 + Z2 %*% G4 - L %*% ( Lambda %*% G1 %*% G4 - m2 * Lambda %*% G5 )
	U <- t(G3) %*% G3 + m1 * (diag(r) - G6) + t(G4) %*% G4 + m2 * (diag(r) - G7)
	return(c(M %*% solve(U)))    #update for Gamma
}

ll.Gamma <- function(par, A, Z2, Lambda, L, Sigma, m1, m2) {
	p <- nrow(Z1)
	r <- length(par) / p
	Gamma <- matrix(par, nrow=p, ncol=r)
	
	##Compute Preliminary Matrices##
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
	
	logdet <- -m1/2 * log(det(mid.r)) - m2/2 * log(det(mid.Lambda))
	tr1 <- 1/2 * sum( solve(mid.r, t(AtSinvG)) * t(AtSinvG) )
	tr2 <- 1/2 * sum( solve(mid.r, GtSinvZ) * GtSinvZ ) + 1/2 * sum( solve(mid.Lambda, H1) * H1 )
	return((logdet + tr1 + tr2)/p)
}



