##Do we underestimate Omega with GLS?

CompareEstimate <- function(X, alpha, Gamma, Sigma) {     #It is assumed X is n x d where the first column is a vector of 1's
	r <- ncol(Gamma)
	p <- nrow(Gamma)
	d <- ncol(X)
	n <- nrow(X)
	
	n.x <- 1/solve(t(X) %*% X)[2,2]
	
	Q <- qr.Q(qr(X), complete=T)
	R <- qr.R(qr(X))
	Q1 <- Q[,1:d]
	Q2 <- Q[,(d+1):n]
	W <- matrix(rnorm(r*n), nrow=r, ncol=n)
	cov.total <- cbind(X, t(W))
	orthog.total <- diag(n) - cov.total %*% solve(t(cov.total) %*% cov.total, t(cov.total))
	W1 <- W %*% Q1     #r x d
	E <- matrix(rnorm(n*p), nrow=p, ncol=n) * sqrt(Sigma)
	Sigma.hat.ols <- rowSums((E %*% orthog.total) * E) / (n - d)
	E1 <- E %*% Q1
	
	Y <- Gamma %*% (W + alpha %*% t(X[,2:d])) + E
	Y.tilde <- Y %*% Q2
	rm(E)
	out.em <- fa.em( t( Y.tilde ), r=r )
	Gamma.hat <- out.em$Gamma
	W.pca <- svd( t(Y.tilde) %*% Y.tilde )$v[,1:r]    #n-d x r
	orthog.W.pca <- diag(n - d) - W.pca %*% t(W.pca)
	Gamma.hat.pca <- Y.tilde %*% W.pca
	Sigma.hat.pca <- rowSums((Y.tilde %*% orthog.W.pca) * Y.tilde) / (n - d)
	Sigma.hat <- out.em$Sigma
	Info <- t(Gamma.hat / Sigma.hat) %*% Gamma.hat / p
	Info.pca <- svd( t(Gamma.hat.pca / Sigma.hat.pca) %*% Gamma.hat.pca / p )$d
	Rest <- t(t((Gamma) / Sigma.hat) %*% Gamma.hat / p)    #We are interested in the eigenvalues of Rest'Info^{-2}Rest and if they are all strictly less than 1, which assures alpha is UNDERESTIMATED
	Inter <- t( (Gamma.hat - Gamma) / Sigma.hat ) %*% Gamma.hat / p
	n.iter.em <- out.em$niter
	rm(out.em)
	
	alpha.hat <- solve(t(Gamma.hat / Sigma.hat) %*% Gamma.hat, t(Gamma.hat / Sigma.hat) %*% Y %*% Q1 %*% solve(t(R)))[,2:d]
	alpha.true <- alpha + (W1 %*% solve(t(R)))[,2:d]
	
	cate.i <- cate.fit(X.primary = cbind(X[,2]), X.nuis = cbind(X[,1]), Y = t(Y), r=r, fa.method = "ml", calibrate = F)
	V.i <- svd(t(cate.i$Gamma / cate.i$Sigma) %*% cate.i$Gamma)$v
	alpha.rr <- t(V.i) %*% cate.i$alpha
	
	#lad <- rq(Gamma %*% cbind(alpha.true) ~ Gamma.hat, weights=1/sqrt(Sigma))
	
	return(list(alpha.true=alpha.true, alpha.hat=alpha.hat, alpha.rr=alpha.rr, resid=solve(t(Gamma.hat / Sigma.hat) %*% Gamma.hat, t(Gamma.hat / Sigma.hat) %*% E1 %*% solve(t(R)))[,2:d], Info=Info, Rest=Rest, Inter=Inter, Sigma.hat=Sigma.hat, Sigma.ols=Sigma.hat.ols, niter=n.iter.em, Info.pca=Info.pca))
}


CompareEstimate2 <- function(X, alpha, E, Gamma, Sigma) {     #It is assumed X is n x d where the first column is a vector of 1's
                                                       #The only difference between this and the function above is that this takes as input E
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
  E1 <- E %*% Q1
  
  Y <- Gamma %*% (W + alpha %*% t(X[,2:d])) + E
  rm(E)
  out.em <- fa.em( t( Y %*% Q2 ), r=r )
  Gamma.hat <- out.em$Gamma
  Sigma.hat <- out.em$Sigma
  Info <- t(Gamma.hat / Sigma.hat) %*% Gamma.hat / p
  Rest <- t(t((Gamma) / Sigma.hat) %*% Gamma.hat / p)    #We are interested in the eigenvalues of Rest'Info^{-2}Rest and if they are all strictly less than 1, which assures alpha is UNDERESTIMATED
  Inter <- t( (Gamma.hat - Gamma) / Sigma.hat ) %*% Gamma.hat / p
  rm(out.em)
  
  alpha.hat <- solve(t(Gamma.hat / Sigma.hat) %*% Gamma.hat, t(Gamma.hat / Sigma.hat) %*% Y %*% Q1 %*% solve(t(R)))[,2:d]
  alpha.hat.2 <- solve(t(Gamma.hat / Sigma) %*% Gamma.hat, t(Gamma.hat / Sigma) %*% Y %*% Q1 %*% solve(t(R)))[,2:d]
  alpha.true <- alpha + (W1 %*% solve(t(R)))[,2:d]
  return(list(alpha.true=alpha.true, alpha.hat=alpha.hat, alpha.hat.2=alpha.hat.2, resid=solve(t(Gamma.hat / Sigma.hat) %*% Gamma.hat, t(Gamma.hat / Sigma.hat) %*% E1 %*% solve(t(R)))[,2:d], Info=Info, Rest=Rest, Inter=Inter, Sigma.hat=Sigma.hat))
}

###This code will simulate data and analyze it using cate (without cell type info). The results will be used in the Fdr figure###

Sim.Analyze.Fdr <- function(X, C, B, L, Sigma) {   #It is assumed that d = 2 and the first row of X is the intercept
	p <- nrow(L)
	n <- ncol(X)
	r.sim <- ncol(L)
	Y <- B %*% X + L %*% C + matrix(rnorm(n*p), nrow=p, ncol=n) * sqrt(Sigma)
	
	data.X <- t(X)
	colnames(data.X) <- c("Intercept", "Cov1")
	data.X <- data.frame(data.X)
	
	cate.sim <- cate(~Cov1, X.data=data.X, Y=t(Y), r=r.sim, fa.method="ml", adj.method="rr", calibrate=F)
	Z.cate <- cate.sim$Z
	X.tot <- cbind(t(X), Z.cate)
	orthog.X.tot <- diag(n) - X.tot %*% solve(t(X.tot) %*% X.tot, t(X.tot))
	dof <- n - ncol(X.tot)
	B.hat <- (Y %*% X.tot %*% solve(t(X.tot) %*% X.tot))[,2]
	Sigma.hat <- rowSums((Y %*% orthog.X.tot) * Y) / dof
	var.cov1 <- Sigma.hat * solve(t(X.tot) %*% X.tot)[2,2]
	z.cov1 <- B.hat / sqrt(var.cov1)
	p.cov1 <- 2 - 2 * pt(abs(z.cov1), df=dof)
	q.cov1 <- qvalue(p.cov1)
	fdp.sim <- false.sign.results(B[,2], B.hat, q.cov1$qvalue)$fdr
	
	all.data <- cbind(t(X), t(C))
	B.hat.cell <- (Y %*% all.data %*% solve(t(all.data) %*% all.data))[,2]
	orthog.all.data <- diag(n) - all.data %*% solve(t(all.data) %*% all.data, t(all.data))
	Sigma.hat.cell <- rowSums((Y %*% orthog.all.data) * Y) / dof
	var.cov1.cell <- Sigma.hat.cell * solve(t(all.data) %*% all.data)[2,2]
	p.cov1.cell <- 2 - 2 * pt(abs(B.hat.cell / sqrt(var.cov1.cell)), df=dof)
	q.cov1.cell <- qvalue(p.cov1.cell)
	fdp.sim.cell <- false.sign.results(B[,2], B.hat.cell, q.cov1.cell$qvalue)$fdr
	
	return(list(fdp.nocell=fdp.sim, q.nocell=q.cov1, alpha.cate=cate.sim$alpha, fdp.cell=fdp.sim.cell, q.cell=q.cov1.cell))
}

CompareEstimate3 <- function(X, alpha, Gamma, Sigma) {     #It is assumed X is n x d where the first column is a vector of 1's
  r <- ncol(Gamma)
  p <- nrow(Gamma)
  d <- ncol(X)
  n <- nrow(X)
  
  Q <- qr.Q(qr(X), complete=T)
  R <- qr.R(qr(X))
  Q1 <- Q[,1:d]
  Q2 <- Q[,(d+1):n]
  W <- matrix(rnorm(r*n), nrow=r, ncol=n)
  cov.total <- cbind(X, t(W))
  orthog.total <- diag(n) - cov.total %*% solve(t(cov.total) %*% cov.total, t(cov.total))
  W1 <- W %*% Q1     #r x d
  E <- matrix(rnorm(n*p), nrow=p, ncol=n) * sqrt(Sigma)
  Sigma.hat.ols <- rowSums((E %*% orthog.total) * E) / (n - d)
  E1 <- E %*% Q1
  
  Y <- Gamma %*% (W + alpha %*% t(X[,2:d])) + E
  cate.out <- cate.fit(X.primary = cbind(X[,2]), X.nuis = cbind(X[,1]), Y = t(Y), r = r, fa.method = "pc", calibrate = F)
  
  Y.tilde <- Y %*% Q2
  rm(E)
  out.em <- fa.em( t( Y.tilde ), r=r )
  Gamma.hat <- out.em$Gamma
  W.pca <- svd( t(Y.tilde) %*% Y.tilde )$v[,1:r]    #n-d x r
  orthog.W.pca <- diag(n - d) - W.pca %*% t(W.pca)
  Gamma.hat.pca <- Y.tilde %*% W.pca
  Sigma.hat.pca <- rowSums((Y.tilde %*% orthog.W.pca) * Y.tilde) / (n - d)
  Sigma.hat <- out.em$Sigma
  Info <- t(Gamma.hat / Sigma.hat) %*% Gamma.hat / p
  Info.pca <- svd( t(Gamma.hat.pca / Sigma.hat.pca) %*% Gamma.hat.pca / p )$d
  Rest <- t(t((Gamma) / Sigma.hat) %*% Gamma.hat / p)    #We are interested in the eigenvalues of Rest'Info^{-2}Rest and if they are all strictly less than 1, which assures alpha is UNDERESTIMATED
  Inter <- t( (Gamma.hat - Gamma) / Sigma.hat ) %*% Gamma.hat / p
  n.iter.em <- out.em$niter
  rm(out.em)
  
  alpha.hat <- solve(t(Gamma.hat / Sigma.hat) %*% Gamma.hat, t(Gamma.hat / Sigma.hat) %*% Y %*% Q1 %*% solve(t(R)))[,2:d]
  alpha.true <- alpha + (W1 %*% solve(t(R)))[,2:d]
  return(list(alpha.true=alpha.true, alpha.hat=alpha.hat, resid=solve(t(Gamma.hat / Sigma.hat) %*% Gamma.hat, t(Gamma.hat / Sigma.hat) %*% E1 %*% solve(t(R)))[,2:d], Info=Info, Rest=Rest, Inter=Inter, Sigma.hat=Sigma.hat, Sigma.ols=Sigma.hat.ols, niter=n.iter.em, Info.pca=Info.pca, cate.out=cate.out))
}
