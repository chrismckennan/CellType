####This code attempts to perform ML estimation with partially observed covariates with the EM algorithm####


#F.mat is K x m1
PEM <- function(Z1, Z2, F.mat, L, Lambda, Gamma, Sigma, m1, m2, update.Gamma=T, update.Sigma=F, update.Lambda=F, max.iter=1e4, tol=1e-6) {
	if (!is.null(Gamma)) {
		Gamma.old <- Gamma
	} else {
		Gamma <- matrix(rep(0,p), nrow=p, ncol=1)
		Gamma.old <- Gamma
	}
	
	p <- nrow(Z1)
	K <- nrow(F.mat)
	r <- ncol(Gamma)
	Z1Ft <- Z1 %*% t(F.mat)      #Z1 F'
	SFF <- F.mat %*% t(F.mat)  #FF'
	tmp.eig <- eigen(Lambda)
	D <- tmp.eig$vectors %*% diag(sqrt(tmp.eig$values)) %*% t(tmp.eig$vectors)
	
	
	##Compute Preliminary Matrices##
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
	mid.Lambda <- D %*% solve(diag(K) + t(D) %*% H2 %*% D, t(D))     #D(I_K + D'H2D)^{-1}D'
	G1 <- H1 - H2 %*% mid.Lambda %*% H1            #K x m2
	G2 <- H2 - H2 %*% mid.Lambda %*% H2            #K x K
	G3 <- H3															      #m1 x r
	G4 <- H4 - t(H1) %*% mid.Lambda %*% H5        #m2 x r
	G5 <- H5 - H2 %*% mid.Lambda %*% H5           #K x r
	G6 <- H6																  #r x r
	G7 <- H6 - t(H5) %*% mid.Lambda %*% H5        #r x r
	
	##Compute Initial Log-likelihood##
	logdet.0 <- -(m1 + m2) / 2 * ( log(det(mid.r)) + sum(log(Sigma)) ) - m2 / 2 * ( log(det(diag(K) + Lambda %*% H2)) )
	tr1.0 <- -1/2 * sum( SinvA * A ) + 1/2 * sum( solve(mid.r, t(AtSinvG)) * t(AtSinvG) )
	tr2.0 <- -1/2 * sum( SinvZ * Z2 ) + 1/2 * sum( solve(mid.r, GtSinvZ) * GtSinvZ ) + 1/2 * sum( (mid.Lambda %*% H1) * H1 )
	ll.0 <- logdet.0 + tr1.0 + tr2.0
	
	ll.vec <- rep(NA, max.iter+1)
	ll.vec[1] <- ll.0
	
	for (iter in 1:max.iter) {
		#####   Update L   #####
		M <- Z1Ft + Z2 %*% t(G1) %*% Lambda
		U <- SFF + Lambda %*% G1 %*% t(G1) %*% Lambda + m2 * (Lambda - Lambda %*% G2 %*% Lambda)
		L.old <- L
		L <- M %*% solve(U)     #EM update for L
		
		#Update Matrices#
		A <- Z1 - L %*% F.mat
		SinvL <- L / Sigma		
		
		LtSinvZ <- t(SinvL) %*% Z2
		LtSinvG <- t(SinvL) %*% Gamma
		LtSinvL <- t(L) %*% SinvL
		AtSinvG <- t(A) %*% SinvG
		
		#Update H#
		H1 <- LtSinvZ - LtSinvG %*% solve(mid.r, GtSinvZ)
		H2 <- LtSinvL - LtSinvG %*% solve(mid.r, t(LtSinvG))
		H3 <- AtSinvG - AtSinvG %*% solve(mid.r, GtSinvG)         
		H5 <- LtSinvG - LtSinvG %*% solve(mid.r, GtSinvG)         
		
		#Update G#
		mid.Lambda <- D %*% solve(diag(K) + t(D) %*% H2 %*% D, t(D))     #D(I_K + D'H2D)^{-1}D'
		G1 <- H1 - H2 %*% mid.Lambda %*% H1           
		G2 <- H2 - H2 %*% mid.Lambda %*% H2            
		G3 <- H3															      
		G4 <- H4 - t(H1) %*% mid.Lambda %*% H5  
		G5 <- H5 - H2 %*% mid.Lambda %*% H5
		G7 <- H6 - t(H5) %*% mid.Lambda %*% H5
		
		
		if (update.Lambda) {    #Only update Lambda if the estimated covariance has full numerical rank. Otherwise the solution will NOT converge.
			#####   Update Lambda   #####
			tmp.eig <- eigen(Lambda)
			R.vec.start <- c( tmp.eig$vectors %*% diag(sqrt(tmp.eig$values)) %*% t(tmp.eig$vectors) )      #Make starting matrix symmetric
			D <- matrix(optim(R.vec.start, fn=fun.bfgs, gr=grad.bfgs, method="BFGS", H1=H1, H2=H2, m2=m2)$par, nrow=K, ncol=K)
			Lambda <- D %*% t(D)        #D should be symmetric. This is a precaution
			
			#Update G#
			mid.Lambda <- D %*% solve(diag(K) + t(D) %*% H2 %*% D, t(D))
			G1 <- H1 - H2 %*% mid.Lambda %*% H1           
			G2 <- H2 - H2 %*% mid.Lambda %*% H2            												      
			G4 <- H4 - t(H1) %*% mid.Lambda %*% H5  
			G5 <- H5 - H2 %*% mid.Lambda %*% H5
			G7 <- H6 - t(H5) %*% mid.Lambda %*% H5
		}	
		
		
		#####   Update Gamma   #####
		if (update.Gamma && !is.null(Gamma)) {
			M <- A %*% G3 + Z2 %*% G4 - L %*% ( Lambda %*% G1 %*% G4 - m2 * Lambda %*% G5 )
			U <- t(G3) %*% G3 + m1 * (diag(r) - G6) + t(G4) %*% G4 + m2 * (diag(r) - G7)
			Gamma.old <- Gamma
			Gamma <- M %*% solve(U)
			
			#Update Matrices#
			SinvG <- Gamma / Sigma
			
			GtSinvZ <- t(SinvG) %*% Z2
			LtSinvG <- t(SinvL) %*% Gamma
			GtSinvG <- t(Gamma) %*% SinvG
			AtSinvG <- t(A) %*% SinvG
			mid.r <- diag(r) + GtSinvG	
	
			#Update H matrices#
			H1 <- LtSinvZ - LtSinvG %*% solve(mid.r, GtSinvZ)
			H2 <- LtSinvL - LtSinvG %*% solve(mid.r, t(LtSinvG))
			H3 <- AtSinvG - AtSinvG %*% solve(mid.r, GtSinvG)
			H4 <- t(GtSinvZ) - t(GtSinvZ) %*% solve(mid.r, GtSinvG)
			H5 <- LtSinvG - LtSinvG %*% solve(mid.r, GtSinvG)
			H6 <- GtSinvG - GtSinvG %*% solve(mid.r, GtSinvG)
	
			#Update G matrices#
			mid.Lambda <- D %*% solve(diag(K) + t(D) %*% H2 %*% D, t(D))
			G1 <- H1 - H2 %*% mid.Lambda %*% H1 
			G2 <- H2 - H2 %*% mid.Lambda %*% H2            
			G3 <- H3															      
			G4 <- H4 - t(H1) %*% mid.Lambda %*% H5        
			G5 <- H5 - H2 %*% mid.Lambda %*% H5           
			G6 <- H6																  
			G7 <- H6 - t(H5) %*% mid.Lambda %*% H5        
		}
		
		
		if (update.Sigma) {
			#####   Update Sigma   #####
			#There are 4 parts to Sigma (see LaTeX document)
			s1 <- rowSums(A * A)
			s2 <- rowSums(Z2 * Z2)
			mat3.g <- t(G3) %*% G3 + m1 * (diag(r) - G6) + t(G4) %*% G4 + m2 * (diag(r) - G7)
			mat3.l <- Lambda %*% G1 %*% G4 - m2 * Lambda %*% G5
			mat.left.3 <- Gamma %*% mat3.g + 2 * L %*% mat3.l - 2 * A %*% G3 - 2 * Z2 %*% G4
			s3 <- rowSums( mat.left.3 * Gamma )
			mat4.l <- Lambda %*% G1 %*% t(G1) %*% Lambda + m2 * (Lambda - Lambda %*% G2 %*% Lambda)
			mat.left.4 <- L %*% mat4.l - 2 * Z2 %*% t(G1) %*% Lambda
			s4 <- rowSums(mat.left.4 * L)
			Sigma.old <- Sigma
			Sigma <- 1/(m1 + m2) * (s1 + s2 + s3 + s4)
			
			#Update Matrices#
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
			
			#Compute H matrices#
			H1 <- LtSinvZ - LtSinvG %*% solve(mid.r, GtSinvZ)
			H2 <- LtSinvL - LtSinvG %*% solve(mid.r, t(LtSinvG))
			H3 <- AtSinvG - AtSinvG %*% solve(mid.r, GtSinvG)
			H4 <- t(GtSinvZ) - t(GtSinvZ) %*% solve(mid.r, GtSinvG)
			H5 <- LtSinvG - LtSinvG %*% solve(mid.r, GtSinvG)
			H6 <- GtSinvG - GtSinvG %*% solve(mid.r, GtSinvG)
			
			#Compute G matrices#
			mid.Lambda <- D %*% solve(diag(K) + t(D) %*% H2 %*% D, t(D))
			G1 <- H1 - H2 %*% mid.Lambda %*% H1 
			G2 <- H2 - H2 %*% mid.Lambda %*% H2            
			G3 <- H3															      
			G4 <- H4 - t(H1) %*% mid.Lambda %*% H5        
			G5 <- H5 - H2 %*% mid.Lambda %*% H5           
			G6 <- H6																  
			G7 <- H6 - t(H5) %*% mid.Lambda %*% H5  
		}

		##Compute Log-likelihood##
		logdet.1 <- -(m1 + m2) / 2 * ( log(det(mid.r)) + sum(log(Sigma)) ) - m2 / 2 * ( log(det(diag(K) + Lambda %*% H2)) )
		tr1.1 <- -1/2 * sum( SinvA * A ) + 1/2 * sum( solve(mid.r, t(AtSinvG)) * t(AtSinvG) )
		tr2.1 <- -1/2 * sum( SinvZ * Z2 ) + 1/2 * sum( solve(mid.r, GtSinvZ) * GtSinvZ ) + 1/2 * sum( (mid.Lambda %*% H1) * H1 )
		ll.1 <- logdet.1 + tr1.1 + tr2.1
		ll.vec[iter + 1] <- ll.1
		
		if (ll.1 < ll.0) {
			#iter; break
			#tmp.matrices <- Compute.Matrices(L.old, Lambda.old, Gamma.old, Sigma.old, Z1, Z2, F.mat)
			#return( list( L=L.old, Lambda=Lambda.old, Gamma=Gamma.old, Sigma=Sigma.old, n.iter=iter-1, ll=ll.vec, diff=rel.diff, out=1/2, blup.V=Lambda.old %*% tmp.matrices$G1, blup.W1=t(tmp.matrices$G3), blup.W2=t(tmp.matrices$G4) ) )
			return( list( L=L.old, Lambda=Lambda.old, Gamma=Gamma.old, Sigma=Sigma.old, n.iter=iter-1, ll=ll.vec, diff=rel.diff, out=1/2, blup.V=NULL, blup.W1=NULL, blup.W2=NULL ) )
		}
		rel.diff <- abs(ll.1 - ll.0) / abs(ll.0)
		if (rel.diff < tol) {
			#Calculate estimates for latent factors#
			
			#W1#
			GtSinvZ1 <- t(Gamma) %*% ( Z1 / Sigma )
			orthog.F.mat <- diag(ncol(Z1)) - t(F.mat) %*% solve(SFF, F.mat)
			W1.hat <- solve(GtSinvG, GtSinvZ1 %*% orthog.F.mat)
			
			#V and W2#
			tmp.mat <- rbind( cbind( LtSinvL, LtSinvG ), cbind( t(LtSinvG), GtSinvG ) )
			VandW2 <- solve(tmp.mat, rbind( LtSinvZ, GtSinvZ ))
			
			return( list( L=L, Lambda=Lambda, Gamma=Gamma, Sigma=Sigma, n.iter=iter, ll=ll.vec, diff=rel.diff, out=1, blup.V=VandW2[1:K,], blup.W1=W1.hat, blup.W2=VandW2[(K+1):(r+K),] ) )
		}
		ll.0 <- ll.1	
		print(as.character(iter))
	}

return( list( L=L, Lambda=Lambda, Gamma=Gamma, Sigma=Sigma, n.iter=iter, ll=ll.vec, diff=rel.diff, out=0, blup.V=Lambda %*% G1, blup.W1=t(G3), blup.W2=t(G4) ) )	
}

#T.mat is s x K, where it is assumed s < K
PEM.reduced <- function(Z1, Z2, F.mat, T.mat=NULL, L, Lambda.0, Gamma, Sigma, m1, m2, update.Gamma=T, update.Sigma=F, update.Lambda=F, update.T=F, blups=T, max.iter=1e4, tol=1e-6) {
	if (!is.null(Gamma)) {
		Gamma.old <- Gamma
	} else {
		Gamma <- matrix(rep(0,p), nrow=p, ncol=1)
		Gamma.old <- Gamma
	}
	
	K <- nrow(F.mat)
	if (is.null(T.mat)) {
		T.mat <- diag(1, K, K)
	}
	Lambda <- T.mat %*% Lambda.0 %*% t(T.mat)    #This is Lambda.tilde in the Latex document
	
	s <- nrow(T.mat)
	p <- nrow(Z1)
	r <- ncol(Gamma)
	Z1Ft <- Z1 %*% t(F.mat)      #Z1F'
	SFF <- F.mat %*% t(F.mat)  #FF'
	
	
	##Compute Preliminary Matrices##
	A <- Z1 - L %*% T.mat %*% F.mat
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
	
	##Compute G matrices##
	mid.Lambda <- solve(Lambda) + H2
	G1 <- H1 - H2 %*% solve(mid.Lambda, H1)         #s x m2
	G2 <- H2 - H2 %*% solve(mid.Lambda, H2)         #s x s
	G3 <- H3															      #m1 x r
	G4 <- H4 - t(H1) %*% solve(mid.Lambda, H5)     #m2 x r
	G5 <- H5 - H2 %*% solve(mid.Lambda, H5)        #s x r
	G6 <- H6																  #r x r
	G7 <- H6 - t(H5) %*% solve(mid.Lambda, H5)     #r x r
	
	##Compute Initial Log-likelihood##
	logdet.0 <- -(m1 + m2) / 2 * ( log(det(mid.r)) + sum(log(Sigma)) ) - m2 / 2 * ( log(det(Lambda)) + log(det(mid.Lambda)) )
	tr1.0 <- -1/2 * sum( SinvA * A ) + 1/2 * sum( solve(mid.r, t(AtSinvG)) * t(AtSinvG) )
	tr2.0 <- -1/2 * sum( SinvZ * Z2 ) + 1/2 * sum( solve(mid.r, GtSinvZ) * GtSinvZ ) + 1/2 * sum( solve(mid.Lambda, H1) * H1 )
	ll.0 <- logdet.0 + tr1.0 + tr2.0
	
	ll.vec <- rep(NA, max.iter+1)
	ll.vec[1] <- ll.0
	
	for (iter in 1:max.iter) {
		#####   Update L   #####
		M <- Z1Ft %*% t(T.mat) + Z2 %*% t(G1) %*% Lambda
		U <- T.mat %*% SFF %*% t(T.mat) + Lambda %*% G1 %*% t(G1) %*% Lambda + m2 * (Lambda - Lambda %*% G2 %*% Lambda)
		L <- M %*% solve(U)     #EM update for L
		
		#Update Matrices#
		SinvL <- L / Sigma		
		
		LtSinvZ1 <- t(SinvL) %*% Z1
		LtSinvZ <- t(SinvL) %*% Z2
		LtSinvG <- t(SinvL) %*% Gamma
		LtSinvL <- t(L) %*% SinvL
		
		#Update H#
		H1 <- LtSinvZ - LtSinvG %*% solve(mid.r, GtSinvZ)
		H2 <- LtSinvL - LtSinvG %*% solve(mid.r, t(LtSinvG))
		H7 <- LtSinvZ1 - LtSinvG %*% solve(mid.r, GtSinvZ1)
		
		########  Update T  ########		
		T.mat <- matrix(optim(T.mat, fn=Tfun.bfgs, gr=Tgrad.bfgs, method="BFGS", H1=H1, H2=H2, H7=H7, F.mat=F.mat, SFF=SFF, Lambda.0=Lambda.0, m2=m2)$par, nrow=s, ncol=K)
		
		#Update A and Lambda#
		A <- Z1 - L %*% T.mat %*% F.mat
		Lambda <- T.mat %*% Lambda.0 %*% t(T.mat)
		
		#Update Matrices inolving A#
		SinvA <- A / Sigma
		AtSinvG <- t(A) %*% SinvG
		
		#Update H#
		H3 <- AtSinvG - AtSinvG %*% solve(mid.r, GtSinvG)         
		H5 <- LtSinvG - LtSinvG %*% solve(mid.r, GtSinvG)         
		
		#Update G#
		mid.Lambda <- solve(Lambda) + H2
		G1 <- H1 - H2 %*% solve(mid.Lambda, H1)
		G2 <- H2 - H2 %*% solve(mid.Lambda, H2)
		G3 <- H3
		G4 <- H4 - t(H1) %*% solve(mid.Lambda, H5)
		G5 <- H5 - H2 %*% solve(mid.Lambda, H5)
		G7 <- H6 - t(H5) %*% solve(mid.Lambda, H5)
		
		
		if (update.Lambda) {    #Only update Lambda if the estimated covariance has full numerical rank. Otherwise the solution will NOT converge.
			#####   Update Lambda   #####
			Lambda.old <- Lambda
			tmp.eig <- eigen(Lambda)
			R.vec.start <- c( tmp.eig$vectors %*% diag(sqrt(tmp.eig$values), nrow=s, ncol=s) %*% t(tmp.eig$vectors) )
			R.vec.out <- matrix(optim(R.vec.start, fn=fun.bfgs, gr=grad.bfgs, method="BFGS", H1=H1, H2=H2, m2=m2)$par, nrow=s, ncol=s)
			Lambda <- R.vec.out %*% t(R.vec.out)
			
			#Update G#
			mid.Lambda <- solve(Lambda) + H2
			G1 <- H1 - H2 %*% solve(mid.Lambda, H1)
			G2 <- H2 - H2 %*% solve(mid.Lambda, H2)
			G4 <- H4 - t(H1) %*% solve(mid.Lambda, H5)
			G5 <- H5 - H2 %*% solve(mid.Lambda, H5)
			G7 <- H6 - t(H5) %*% solve(mid.Lambda, H5)	
		}	
		
		
		#####   Update Gamma   #####
		if (update.Gamma && !is.null(Gamma.0)) {
			M <- A %*% G3 + Z2 %*% G4 - L %*% ( Lambda %*% G1 %*% G4 - m2 * Lambda %*% G5 )
			U <- t(G3) %*% G3 + m1 * (diag(r) - G6) + t(G4) %*% G4 + m2 * (diag(r) - G7)
			Gamma.old <- Gamma
			Gamma <- M %*% solve(U)
			
			#Update Matrices#
			SinvG <- Gamma / Sigma
			
			GtSinvZ1 <- t(SinvG) %*% Z1
			GtSinvZ <- t(SinvG) %*% Z2
			LtSinvG <- t(SinvL) %*% Gamma
			GtSinvG <- t(Gamma) %*% SinvG
			AtSinvG <- t(A) %*% SinvG
			mid.r <- diag(r) + GtSinvG	
	
			#Update H matrices#
			H1 <- LtSinvZ - LtSinvG %*% solve(mid.r, GtSinvZ)
			H2 <- LtSinvL - LtSinvG %*% solve(mid.r, t(LtSinvG))
			H3 <- AtSinvG - AtSinvG %*% solve(mid.r, GtSinvG)
			H4 <- t(GtSinvZ) - t(GtSinvZ) %*% solve(mid.r, GtSinvG)
			H5 <- LtSinvG - LtSinvG %*% solve(mid.r, GtSinvG)
			H6 <- GtSinvG - GtSinvG %*% solve(mid.r, GtSinvG)
			H7 <- LtSinvZ1 - LtSinvG %*% solve(mid.r, GtSinvZ1)
	
			#Update G matrices#
			mid.Lambda <- solve(Lambda) + H2
			G1 <- H1 - H2 %*% solve(mid.Lambda, H1)
			G2 <- H2 - H2 %*% solve(mid.Lambda, H2)
			G3 <- H3
			G4 <- H4 - t(H1) %*% solve(mid.Lambda, H5)
			G5 <- H5 - H2 %*% solve(mid.Lambda, H5)
			G6 <- H6
			G7 <- H6 - t(H5) %*% solve(mid.Lambda, H5)
		}
		
		
		if (update.Sigma) {
			#####   Update Sigma   #####
			#There are 4 parts to Sigma (see LaTeX document)
			s1 <- rowSums(A * A)
			s2 <- rowSums(Z2 * Z2)
			mat3.g <- t(G3) %*% G3 + m1 * (diag(r) - G6) + t(G4) %*% G4 + m2 * (diag(r) - G7)
			mat3.l <- Lambda %*% G1 %*% G4 - m2 * Lambda %*% G5
			mat.left.3 <- Gamma %*% mat3.g + 2 * L %*% mat3.l - 2 * A %*% G3 - 2 * Z2 %*% G4
			s3 <- rowSums( mat.left.3 * Gamma )
			mat4.l <- Lambda %*% G1 %*% t(G1) %*% Lambda + m2 * (Lambda - Lambda %*% G2 %*% Lambda)
			mat.left.4 <- L %*% mat4.l - 2 * Z2 %*% t(G1) %*% Lambda
			s4 <- rowSums(mat.left.4 * L)
			Sigma.old <- Sigma
			Sigma <- 1/(m1 + m2) * (s1 + s2 + s3 + s4)
			
			#Update Matrices#
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
			
			#Compute H matrices#
			H1 <- LtSinvZ - LtSinvG %*% solve(mid.r, GtSinvZ)
			H2 <- LtSinvL - LtSinvG %*% solve(mid.r, t(LtSinvG))
			H3 <- AtSinvG - AtSinvG %*% solve(mid.r, GtSinvG)
			H4 <- t(GtSinvZ) - t(GtSinvZ) %*% solve(mid.r, GtSinvG)
			H5 <- LtSinvG - LtSinvG %*% solve(mid.r, GtSinvG)
			H6 <- GtSinvG - GtSinvG %*% solve(mid.r, GtSinvG)
			H7 <- LtSinvZ1 - LtSinvG %*% solve(mid.r, GtSinvZ1)
			
			#Compute G matrices#
			mid.Lambda <- solve(Lambda) + H2
			G1 <- H1 - H2 %*% solve(mid.Lambda, H1)
			G2 <- H2 - H2 %*% solve(mid.Lambda, H2)
			G3 <- H3
			G4 <- H4 - t(H1) %*% solve(mid.Lambda, H5)
			G5 <- H5 - H2 %*% solve(mid.Lambda, H5)
			G6 <- H6
			G7 <- H6 - t(H5) %*% solve(mid.Lambda, H5)
		}

		##Compute Log-likelihood##
		logdet.1 <- -(m1 + m2) / 2 * ( log(det(mid.r)) + sum(log(Sigma)) ) - m2 / 2 * ( log(det(Lambda)) + log(det(mid.Lambda)) )
		tr1.1 <- -1/2 * sum( SinvA * A ) + 1/2 * sum( solve(mid.r, t(AtSinvG)) * t(AtSinvG) )
		tr2.1 <- -1/2 * sum( SinvZ * Z2 ) + 1/2 * sum( solve(mid.r, GtSinvZ) * GtSinvZ ) + 1/2 * sum( solve(mid.Lambda, H1) * H1 )
		ll.1 <- logdet.1 + tr1.1 + tr2.1
		ll.vec[iter + 1] <- ll.1
		
		if (ll.1 < ll.0) {
			#iter; break
			#tmp.matrices <- Compute.Matrices(L.old, Lambda.old, Gamma.old, Sigma.old, Z1, Z2, F.mat)
			#return( list( L=L.old, Lambda=Lambda.old, Gamma=Gamma.old, Sigma=Sigma.old, n.iter=iter-1, ll=ll.vec, diff=rel.diff, out=1/2, blup.V=Lambda.old %*% tmp.matrices$G1, blup.W1=t(tmp.matrices$G3), blup.W2=t(tmp.matrices$G4) ) )
			return( list( L=L.old, Lambda=Lambda.old, Gamma=Gamma.old, Sigma=Sigma.old, n.iter=iter-1, ll=ll.vec, diff=rel.diff, out=1/2, blup.V=NULL, blup.W1=NULL, blup.W2=NULL ) )
		}
		rel.diff <- abs(ll.1 - ll.0) / abs(ll.0)
		if (rel.diff < tol) {
			if (blups) {
				#Calculate estimates for latent factors#
				
				#W1#
				GtSinvZ1 <- t(Gamma) %*% ( Z1 / Sigma )
				orthog.F.mat <- diag(ncol(Z1)) - t(T.mat %*% F.mat) %*% solve(T.mat %*% SFF %*% t(T.mat), T.mat %*% F.mat)
				W1.hat <- solve(GtSinvG, GtSinvZ1 %*% orthog.F.mat)
				
				#V and W2#
				tmp.mat <- rbind( cbind( LtSinvL, LtSinvG ), cbind( t(LtSinvG), GtSinvG ) )
				VandW2 <- solve(tmp.mat, rbind( LtSinvZ, GtSinvZ ))
				
				return( list( L=L, Lambda=Lambda, Gamma=Gamma, Sigma=Sigma, n.iter=iter, ll=ll.vec, diff=rel.diff, out=1, blup.V=VandW2[1:K,], blup.W1=W1.hat, blup.W2=VandW2[(K+1):(r+K),] ) )			
			} else {
				return( list( L=L, Lambda=Lambda, Gamma=Gamma, Sigma=Sigma, n.iter=iter, ll=ll.vec[1:(iter+1)], diff=rel.diff, out=1 ) )
			}
		}
		ll.0 <- ll.1	
		
	}

	if (blups) {
		return( list( L=L, Lambda=Lambda, Gamma=Gamma, Sigma=Sigma, n.iter=iter, ll=ll.vec, diff=rel.diff, out=0, blup.V=Lambda %*% G1, blup.W1=t(G3), blup.W2=t(G4) ) )
	}	else {
		return( list( L=L, Lambda=Lambda, Gamma=Gamma, Sigma=Sigma, n.iter=iter, ll=ll.vec[1:(iter+1)], diff=rel.diff, out=0 ) )
	}
}

Starting.Points <- function(Z1, Z2, F.mat, Sigma, L.0=NULL, Gamma.0=NULL, r.nocell, r.cell, update.T=T) {
	K <- nrow(F.mat)
	m1 <- ncol(Z1)
	m2 <- ncol(Z2)
	Lambda.0 <- 1/m1 * F.mat %*% t(F.mat)
	
	orthog.F <- diag(m1) - t(F.mat) %*% solve(F.mat %*% t(F.mat)) %*% F.mat
	if (is.null(Gamma.0)) {
		svd.Z1.tilde <- svd(Z1 %*% orthog.F)
		Gamma.0 <- svd.Z1.tilde$u[,1:r.cell] %*% diag( svd.Z1.tilde$d[1:r.cell], r.cell, r.cell ) / sqrt(m1 - K)      #Initial estimate for Gamma	
	}
	
	#Compute starting point for T.mat#
	if (update.T) {
		s <- max(r.nocell - r.cell, 1)
		mid.mat <- t(Z1) %*% (Z1 / Sigma) - t(Z1) %*% (Gamma.0 / Sigma) %*% solve( diag(r.cell) + t(Gamma.0) %*% (Gamma.0 / Sigma) ) %*% t(Gamma.0) %*% (Z1 / Sigma)
		tmp.eig <- eigen(Lambda.0)
		S.half <- tmp.eig$vectors %*% diag(sqrt(tmp.eig$values), nrow=K, ncol=K) %*% t(tmp.eig$vectors)
		tmp.T.mat <- matrix(solve(S.half, svd(solve(S.half) %*% F.mat %*% mid.mat %*% t(F.mat) %*% solve(S.half) / p)$v[,1:s]), nrow=s, ncol=K)
		T.mat.0 <- t( qr.Q(qr(t(tmp.T.mat))) )
		F.tilde <- T.mat.0 %*% F.mat
	} else {
		s <- K
		T.mat.0 <- diag(K)
		F.tilde <- F.mat
	}
	Lambda.tilde <- T.mat.0 %*% Lambda.0 %*% t(T.mat.0)
	
	#Use T.mat to compute L.0#
	if (is.null(L.0)) {
		L.0 <- Z1 %*% t(F.tilde) %*% solve(F.tilde %*% t(F.tilde))
	}
	
	#Compute matrices#
	
	SinvL <- L.0 / Sigma
	SinvG <- Gamma.0 / Sigma
	SinvZ <- Z2 / Sigma
	
	LtSinvZ1 <- t(SinvL) %*% Z1
	LtSinvZ <- t(SinvL) %*% Z2
	GtSinvZ1 <- t(SinvG) %*% Z1
	GtSinvZ <- t(SinvG) %*% Z2
	LtSinvG <- t(SinvL) %*% Gamma.0
	GtSinvG <- t(Gamma.0) %*% SinvG
	LtSinvL <- t(L.0) %*% SinvL
	mid.r <- diag(r.cell) + GtSinvG
	
	H1 <- LtSinvZ - LtSinvG %*% solve( mid.r,  GtSinvZ)
	H2 <- LtSinvL - LtSinvG %*% solve( mid.r, t(LtSinvG) )
	#H7 <- LtSinvZ1 - LtSinvG %*% solve(mid.r, GtSinvZ1)
	
	F2.hat <- solve(H2, H1)    #Estimated cell type composition for the second set of individuals
	
	mat.cell <- rbind(t(F.tilde), t(F2.hat))
	orthog.cell <- diag(m1 + m2) - mat.cell %*% solve(t(mat.cell) %*% mat.cell) %*% t(mat.cell)
	svd.Z.tilde <- svd(cbind(Z1, Z2) %*% orthog.cell)
	Gamma.0 <- svd.Z.tilde$u[,1:r.cell] %*% diag( svd.Z.tilde$d[1:r.cell], r.cell, r.cell ) / sqrt(m1 + m2 - s)
		
	return(list( Gamma.0=Gamma.0, Lambda.0=Lambda.0, L.0=L.0, T.mat.0=T.mat.0 ))
}


fun.bfgs <- function(par, H1, H2, m2) {    #minus log-likelihood for the components ONLY containing Lambda
	R <- matrix(par, nrow=nrow(H2), ncol=ncol(H2), byrow=F)
	K <- nrow(R)
	lambda <- R %*% t(R)
	mid.lambda <- R %*% solve(diag(1, nrow=K, ncol=K) + t(R) %*% H2 %*% R, t(R))
	return( m2 / 2 * ( log(det(diag(1, nrow=K, ncol=K) + t(R) %*% H2 %*% R)) ) - 1/2 * sum( (mid.lambda %*% H1) * H1 ) )
}

grad.bfgs <- function(par, H1, H2, m2) {   #Gradient of the above minus ll
	R <- matrix(par, nrow=nrow(H2), ncol=ncol(H2), byrow=F)
	K <- nrow(R)
	lambda <- R %*% t(R)
	mid.lambda <- R %*% solve(diag(1, nrow=K, ncol=K) + t(R) %*% H2 %*% R, t(R))
	g1 <- H1 - H2 %*% mid.lambda %*% H1
	ssg1 <- g1 %*% t(g1)
	g2 <- H2 - H2 %*% mid.lambda %*% H2
	
	grad.mat <- m2 * G2 %*% R - ssg1 %*% R
	return(c(grad.mat))
}



user.func <- function(par, H1, H2, m2) {
	return( max(abs( grad.bfgs(par, H1, H2, m2) )) )
}

Tfun.bfgs <- function(T.vec, F.mat, Lambda.0, H1, H2, H7, m2) {
	K <- nrow(Lambda.0)
	s <- length(T.vec) / K
	T.mat <- matrix(T.vec, nrow=s, ncol=K)
	Lambda.tilde <- T.mat %*% Lambda.0 %*% t(T.mat)
	mid.Lambda.tilde <- solve(Lambda.tilde) + H2
	F.tilde <- T.mat %*% F.mat
	return( m2/2 * log(det(Lambda.tilde)) + m2/2 * log(det(mid.Lambda.tilde)) - sum( rowSums(H7 * F.tilde) ) + 1/2 * sum( rowSums(H2 * (F.tilde %*% t(F.tilde))) ) - 1/2 * sum( rowSums( (t(H1) %*% solve(mid.Lambda.tilde)) * t(H1) ) ) )
}

Tgrad.bfgs <- function(T.vec, F.mat, SFF, Lambda.0, H1, H2, H7, m2) {
	K <- nrow(Lambda.0)
	s <- length(T.vec) / K
	T.mat <- matrix(T.vec, nrow=s, ncol=K)
	mid.Lambda.tilde <- solve(T.mat %*% Lambda.0 %*% t(T.mat)) + H2
	G1 <- H1 - H2 %*% solve(mid.Lambda.tilde, H1)
	G2 <- H2 - H2 %*% solve(mid.Lambda.tilde, H2)
	return( c( m2 * G2 %*% T.mat %*% Lambda.0 - G1 %*% t(G1) %*% T.mat %*% Lambda.0 + H2 %*% T.mat %*% SFF - H7 %*% t(F.mat) ) )	
}


Power.Law <- function(Gamma.1.tilde, Gamma.1) {       #Find the maximum eigenvalue and eigenvector of Gamma.1.tilde Gamma.1.tilde' - Gamma.1 Gamma.1' = L.tilde L.tilde'. We can then compute the direction in which to project L onto as v'L
	p <- nrow(Gamma.1)
	v.0 <- rnorm(p)
	v.0 <- v.0 / sqrt(sum(v.0 * v.0))
	v.1 <- Gamma.1.tilde %*% ( t(Gamma.1.tilde) %*% v.0 ) - Gamma.1 %*% ( t(Gamma.1) %*% v.0 )
	norm.1 <- sqrt(sum(v.1 * v.1))
	v.1 <- v.1 / norm.1
	diff <- max(abs(v.1 - v.0))
	v.0 <- v.1
	count <- 1
	while(diff > 1e-12 && count < 1e2) {
		v.1 <- Gamma.1.tilde %*% ( t(Gamma.1.tilde) %*% v.0 ) - Gamma.1 %*% ( t(Gamma.1) %*% v.0 )
		norm.1 <- sqrt(sum(v.1 * v.1))
		v.1 <- v.1 / norm.1
		diff <- max(abs(v.1 - v.0))
		v.0 <- v.1
		count <- count + 1	
	}
	u.1.tilde <- t(Gamma.1.tilde) %*% v.0
	u.1 <- t(Gamma.1) %*% v.0
	max.eig <- sum(u.1.tilde * u.1.tilde) - sum(u.1 * u.1)
	return(list(vec=v.0, eig=max.eig))
}
