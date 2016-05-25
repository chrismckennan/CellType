library('gtools')
##Generate Simulated data
#Omega.cell is K x (d+1)
#sigma.confound is an r-vector of effect sizes for the confounder matrix
#This code simulates data with Cell type distributed as Dirichlet
Sim.Data <- function(n, r, p, pi0.B, pi0.L, ind.BL, sigma.B, sigma.L, sigma.confound, Omega.cell, alpha.cell, alpha.gam, beta.gam) {
	K <- nrow(Omega.cell)
	d <- ncol(Omega.cell) - 1
	Sigma <- rgamma(p, shape=alpha.gam, rate=beta.gam)   #Residual Covariance
	Gamma <- t(sigma.confound * matrix(rnorm(r * p), nrow=r, ncol=p))
	X <- cbind(rep(1, n), matrix( rbinom(n*d, 1, 1/2), nrow=n, ncol=d ))   #An n x (d+1) matrix
  if (ind.BL) {
    B <- matrix( rbinom(p*d, 1, 1-pi0.B), nrow=p, ncol=d, byrow=T ) * t(matrix( rnorm(p*d), nrow=d, ncol=p ) * sigma.B)     #p x d
    L <- matrix( rbinom(p*K, 1, 1-pi0.L), nrow=p, ncol=K, byrow=T ) * t(matrix( rnorm(p*K), nrow=K, ncol=p ) * sigma.L)      #p x K
  } else {
    n.nonzero.B <- floor(p * (1-pi0.B[2]))
    n.nonzero.L <- floor(p * (1-pi0.L[1]))
    B <- cbind( matrix( rbinom(p*(d-1), 1, 1-pi0.B[-2]), nrow=p, ncol=d-1, byrow=T ) * t(matrix( rnorm(p*(d-1)), nrow=d-1, ncol=p ) * sigma.B[-2]), c(sigma.B[2] * rnorm(n.nonzero.B), rep(0, p-n.nonzero.B)) )
    L <- cbind( c(sigma.L[1] * rnorm(n.nonzero.L), rep(0, p-n.nonzero.L)), matrix( rbinom(p*(K-1), 1, 1-pi0.L[-1]), nrow=p, ncol=K-1, byrow=T ) * t(matrix( rnorm(p*(K-1)), nrow=K-1, ncol=p ) * sigma.L[-1]) )
  }
	mean.cell <- Omega.cell %*% t(X)
	C <- matrix(0, nrow=n, ncol=K)     #n x K
	for (i in 1:n) {
		c.i <- mean.cell[,i]
		C[i,] <- rdirichlet(1, alpha.cell * c(c.i, 1-sum(c.i)))[1:K]
	}
	Y <- B %*% t(X[,2:(d+1)]) + L %*% t(C) + Gamma %*% matrix(rnorm(r*n), nrow=r, ncol=n) + matrix(rnorm(n*p), nrow=p, ncol=n) * sqrt(Sigma)
	return(list(M.sim=Y, C.sim=C, X.sim=X, Gamma.sim=Gamma, B.sim=B, L.sim=L, Sigma.sim=Sigma))
}


##The below coded simulates data exactly as above, but draws the cell type confounding from a normal distribution##

Sim.Data_Normal <- function(n, r, p, pi0.B, pi0.L, ind.BL, sigma.B, sigma.L, sigma.confound, Omega.cell, alpha.cell, alpha.gam, beta.gam) {
  K <- nrow(Omega.cell)
  d <- ncol(Omega.cell) - 1
  Sigma <- rgamma(p, shape=alpha.gam, rate=beta.gam)   #Residual Covariance
  if (r > 0) {
    Gamma <- t(sigma.confound * matrix(rnorm(r * p), nrow=r, ncol=p))
  } else {
    Gamma <- 0
  }
  X <- cbind(rep(1, n), matrix( rbinom(n*d, 1, 1/2), nrow=n, ncol=d ))   #An n x (d+1) matrix
  if (ind.BL) {
    B <- matrix( rbinom(p*d, 1, 1-pi0.B), nrow=p, ncol=d, byrow=T ) * t(matrix( rnorm(p*d), nrow=d, ncol=p ) * sigma.B)     #p x d
    L <- matrix( rbinom(p*K, 1, 1-pi0.L), nrow=p, ncol=K, byrow=T ) * t(matrix( rnorm(p*K), nrow=K, ncol=p ) * sigma.L)      #p x K
  } else {
    n.nonzero.B <- floor(p * (1-pi0.B[2]))
    n.nonzero.L <- floor(p * (1-pi0.L[1]))
    B <- cbind( matrix( rbinom(p*(d-1), 1, 1-pi0.B[-2]), nrow=p, ncol=d-1, byrow=T ) * t(matrix( rnorm(p*(d-1)), nrow=d-1, ncol=p ) * sigma.B[-2]), c(sigma.B[2] * rnorm(n.nonzero.B), rep(0, p-n.nonzero.B)) )
    L <- cbind( c(sigma.L[1] * rnorm(n.nonzero.L), rep(0, p-n.nonzero.L)), matrix( rbinom(p*(K-1), 1, 1-pi0.L[-1]), nrow=p, ncol=K-1, byrow=T ) * t(matrix( rnorm(p*(K-1)), nrow=K-1, ncol=p ) * sigma.L[-1]) )
  }
  
  mean.cell <- Omega.cell %*% t(X)
  Var.cell <- 1/(alpha.cell + 1) * ( diag(mean.cell[,1], nrow=K, ncol=K) - cbind(mean.cell[,1]) %*% rbind(mean.cell[,1]) )
  R.cell <- chol(Var.cell)    #R'R = Var.cell
  C <- t(mean.cell) + matrix(rnorm(n*K), nrow=n, ncol=K) %*% R.cell   #n x K

  if (r > 0) {
    Y <- B %*% t(X[,2:(d+1)]) + L %*% t(C) + Gamma %*% matrix(rnorm(r*n), nrow=r, ncol=n) + matrix(rnorm(n*p), nrow=p, ncol=n) * sqrt(Sigma)
  } else {
    Y <- B %*% t(X[,2:(d+1)]) + L %*% t(C) + matrix(rnorm(n*p), nrow=p, ncol=n) * sqrt(Sigma)
  }
  return(list(M.sim=Y, C.sim=C, X.sim=X, Gamma.sim=Gamma, B.sim=B, L.sim=L, Sigma.sim=Sigma))
}


AnalyzeData <- function(M.sim, data.sim, r.cell, r.nocell, B.sim, L.sim) {
	cate.nocell.sim <- cate(~Cov1 + Cov2, X.data=data.sim, Y=t(M.sim), r=r.nocell, fa.method="ml", adj.method="rr", calibrate=F)
	n.sim <- nrow(data.sim)
	Z.nocell.sim <- cate.nocell.sim$Z
	X.nocell.sim <- cbind(rep(1,n.sim), data.sim$Cov1, data.sim$Cov2, Z.nocell.sim)
	d.nocell.sim <- ncol(X.nocell.sim)
	beta.op.nocell.sim <- solve(t(X.nocell.sim) %*% X.nocell.sim) %*% t(X.nocell.sim)
	orthog.X.nocell.sim <- diag(n.sim) - X.nocell.sim %*% solve(t(X.nocell.sim) %*% X.nocell.sim) %*% t(X.nocell.sim)
	B.est.nocell.sim <- M.sim %*% t(beta.op.nocell.sim)    #These will give you the EXACT same beta's. I recommend using the estimated design matrix to calculate sigma^2 for each gene, since you will be able to account for the degrees of freedom used to calculate it
	Sigma.nocell.sim <- 1/(n.sim-d.nocell.sim) * rowSums( M.sim * (M.sim %*% orthog.X.nocell.sim) )
	Var.row.nocell.sim <- solve(t(X.nocell.sim) %*% X.nocell.sim)
	Zscores.nocell.sim <- (B.est.nocell.sim[,2:3] %*% diag(1/sqrt(diag(Var.row.nocell.sim[2:3,2:3]))))/sqrt(Sigma.nocell.sim)
	P.values.nocell.sim <- 2 - 2 * pt( abs(Zscores.nocell.sim), df=n.sim-d.nocell.sim )
	q.nocell.2.sim <- qvalue(P.values.nocell.sim[,2])
	q.nocell.1.sim <- qvalue(P.values.nocell.sim[,1])
	fsr.nocell.2 <- false.sign.results(B.sim[,2], B.est.nocell.sim[,3], q.nocell.2.sim$qvalue)
	#plot(sort(q.nocell.2.sim$qvalue), fsr.nocell.2$fdr, xlab="Est. Q-value", ylab="True FDR", type="l", main="FDR Plot WITHOUT Cell Type for Cov2")
	#plot(fsr.nocell.2$fdr, fsr.nocell.2$power, type="l")
	
	
	#With cell type info#
	#conf.cell.sim <- est.confounder.num(~Cov1 + Cov2 | Tcells + Neutro + Eos, X.data=data.sim, Y=t(M.sim), method="bcv", bcv.plot = T, rmax=10); r.cell.sim <- conf.cell.sim$r
	cate.cell.sim <- cate(~Cov1 + Cov2 | Tcells + Neutro + Eos, X.data=data.sim, Y=t(M.sim), r=r.cell, fa.method="ml", adj.method="rr", calibrate=F)
	Z.cell.sim <- cate.cell.sim$Z
	X.cell.sim <- cbind(rep(1,n.sim), data.sim$Cov1, data.sim$Cov2, data.sim$Tcells, data.sim$Neutro, data.sim$Eos, Z.cell.sim)
	d.cell.sim <- ncol(X.cell.sim)
	beta.op.cell.sim <- solve(t(X.cell.sim) %*% X.cell.sim) %*% t(X.cell.sim)
	orthog.X.cell.sim <- diag(n.sim) - X.cell.sim %*% solve(t(X.cell.sim) %*% X.cell.sim) %*% t(X.cell.sim)
	B.est.cell.sim <- M.sim %*% t(beta.op.cell.sim)    #These will give you the EXACT same beta's. I recommend using the estimated design matrix to calculate sigma^2 for each gene, since you will be able to account for the degrees of freedom used to calculate it
	Sigma.cell.sim <- 1/(n.sim-d.cell.sim) * rowSums( M.sim * (M.sim %*% orthog.X.cell.sim) )
	Var.row.cell.sim <- solve(t(X.cell.sim) %*% X.cell.sim)
	Zscores.cell.sim <- (B.est.cell.sim[,2:3] %*% diag(1/sqrt(diag(Var.row.cell.sim[2:3,2:3]))))/sqrt(Sigma.cell.sim)
	P.values.cell.sim <- 2 - 2 * pt( abs(Zscores.cell.sim), df=n.sim-d.cell.sim )
	q.cell.2.sim <- qvalue(P.values.cell.sim[,2])
	q.cell.1.sim <- qvalue(P.values.cell.sim[,1])
	fsr.cell.2 <- false.sign.results(B.sim[,2], B.est.cell.sim[,3], q.cell.2.sim$qvalue)
	#plot(sort(q.cell.2.sim$qvalue), fsr.cell.2$fdr, xlab="Est. Q-value", ylab="True FDR", type="l", main="FDR Plot WITHOUT Cell Type for Cov2")
	#plot(fsr.cell.2$fdr, fsr.cell.2$power, type="l")
  Zscores.cell.types.sim <- (B.est.cell.sim[,4:6] %*% diag(1/sqrt(diag(Var.row.cell.sim[4:6,4:6]))))/sqrt(Sigma.cell.sim)
  P.values.cell.types.sim <- 2 - 2 * pt( abs(Zscores.cell.types.sim), df=n.sim-d.cell.sim )
  q.cell.neutro.sim <- qvalue(P.values.cell.types.sim[,2])
	
	##Summarize Results
	plot(sort(q.cell.2.sim$qvalue), fsr.cell.2$fdr, xlab="Est. Q-value", ylab="True FDR", main=paste0("FDR Plot With and Without Cell Type, sigma.L = ", paste(as.character(round(100*sigma.L.sim)/100), collapse=",")), type="l")
	lines(sort(q.nocell.2.sim$qvalue), fsr.nocell.2$fdr, col="red")
	legend("bottomright", legend=c("With Cell", "Without Cell"), fill=c("black", "red"))
	abline(a=0,b=1, col="violet")
	
	##Compute risk for true non-zero B[,2]'s and non-zero L[,1]'s##
	ind.nonzero.B <- which(B.sim[,2] != 0)
	ind.nonzero.L <- which(L.sim[,1] != 0)
	out.mat <- rbind(c( sum((B.est.cell.sim[ind.nonzero.B, 3] - B.sim[ind.nonzero.B, 2])^2), sum((B.est.nocell.sim[ind.nonzero.B, 3] - B.sim[ind.nonzero.B, 2])^2), sum((B.est.cell.sim[ind.nonzero.L, 3] - B.sim[ind.nonzero.L, 2])^2),  sum((B.est.nocell.sim[ind.nonzero.L, 3] - B.sim[ind.nonzero.L, 2])^2), length(ind.nonzero.B), length(ind.nonzero.L)))
	colnames(out.mat) <- c("B2 non-zero Cell Risk", "B2 non-zero noCell Risk", "B2 L1non-zero Cell Risk", "B2 L1non-zero noCell Risk", "n_nonzeroB2", "n_nonzeroL1")
	return( list(risk.mat=out.mat, Effects.cell=B.est.cell.sim, Effects.nocell=B.est.nocell.sim, fsr.cell=fsr.cell.2, fsr.nocell=fsr.nocell.2, q.cell=q.cell.2.sim, q.nocell=q.nocell.2.sim, cate.nocell=cate.nocell.sim, cate.cell=cate.cell.sim, q.neutro=q.cell.neutro.sim) )
}


##This data analyzes data from a SINGLE covariate of interest##

AnalyzeData2 <- function(M.sim, data.sim, r.cell, r.nocell, B.sim, L.sim) {
	cate.nocell.sim <- cate(~Cov1, X.data=data.sim, Y=t(M.sim), r=r.nocell, fa.method="ml", adj.method="rr", calibrate=F)
	n.sim <- nrow(data.sim)
	Z.nocell.sim <- cate.nocell.sim$Z
	X.nocell.sim <- cbind(rep(1,n.sim), data.sim$Cov1, Z.nocell.sim)
	d.nocell.sim <- ncol(X.nocell.sim)
	beta.op.nocell.sim <- solve(t(X.nocell.sim) %*% X.nocell.sim) %*% t(X.nocell.sim)
	orthog.X.nocell.sim <- diag(n.sim) - X.nocell.sim %*% solve(t(X.nocell.sim) %*% X.nocell.sim) %*% t(X.nocell.sim)
	B.est.nocell.sim <- M.sim %*% t(beta.op.nocell.sim)    #These will give you the EXACT same beta's. I recommend using the estimated design matrix to calculate sigma^2 for each gene, since you will be able to account for the degrees of freedom used to calculate it
	Sigma.nocell.sim <- 1/(n.sim-d.nocell.sim) * rowSums( M.sim * (M.sim %*% orthog.X.nocell.sim) )
	Var.row.nocell.sim <- solve(t(X.nocell.sim) %*% X.nocell.sim)
	Zscores.nocell.sim <- B.est.nocell.sim[,2:2] * 1/sqrt(Var.row.nocell.sim[2:2,2:2])/sqrt(Sigma.nocell.sim)
	P.values.nocell.sim <- 2 - 2 * pt( abs(Zscores.nocell.sim), df=n.sim-d.nocell.sim )
	#q.nocell.2.sim <- qvalue(P.values.nocell.sim[,2])
	q.nocell.1.sim <- qvalue(P.values.nocell.sim)
	fsr.nocell.2 <- false.sign.results(B.sim[,1], B.est.nocell.sim[,2], q.nocell.1.sim$qvalue)
	#plot(sort(q.nocell.2.sim$qvalue), fsr.nocell.2$fdr, xlab="Est. Q-value", ylab="True FDR", type="l", main="FDR Plot WITHOUT Cell Type for Cov2")
	#plot(fsr.nocell.2$fdr, fsr.nocell.2$power, type="l")
	
	
	#With cell type info#
	#conf.cell.sim <- est.confounder.num(~Cov1 + Cov2 | Tcells + Neutro + Eos, X.data=data.sim, Y=t(M.sim), method="bcv", bcv.plot = T, rmax=10); r.cell.sim <- conf.cell.sim$r
	cate.cell.sim <- cate(~Cov1| Tcells + Neutro + Eos, X.data=data.sim, Y=t(M.sim), r=r.cell, fa.method="ml", adj.method="rr", calibrate=F)
	Z.cell.sim <- cate.cell.sim$Z
	X.cell.sim <- cbind(rep(1,n.sim), data.sim$Cov1, data.sim$Tcells, data.sim$Neutro, data.sim$Eos, Z.cell.sim)
	d.cell.sim <- ncol(X.cell.sim)
	beta.op.cell.sim <- solve(t(X.cell.sim) %*% X.cell.sim) %*% t(X.cell.sim)
	orthog.X.cell.sim <- diag(n.sim) - X.cell.sim %*% solve(t(X.cell.sim) %*% X.cell.sim) %*% t(X.cell.sim)
	B.est.cell.sim <- M.sim %*% t(beta.op.cell.sim)    #These will give you the EXACT same beta's. I recommend using the estimated design matrix to calculate sigma^2 for each gene, since you will be able to account for the degrees of freedom used to calculate it
	Sigma.cell.sim <- 1/(n.sim-d.cell.sim) * rowSums( M.sim * (M.sim %*% orthog.X.cell.sim) )
	Var.row.cell.sim <- solve(t(X.cell.sim) %*% X.cell.sim)
	Zscores.cell.sim <- B.est.cell.sim[,2:2] * 1/sqrt(Var.row.cell.sim[2:2,2:2])/sqrt(Sigma.cell.sim)
	P.values.cell.sim <- 2 - 2 * pt( abs(Zscores.cell.sim), df=n.sim-d.cell.sim )
	#q.cell.2.sim <- qvalue(P.values.cell.sim[,2])
	q.cell.1.sim <- qvalue(P.values.cell.sim)
	fsr.cell.2 <- false.sign.results(B.sim[,1], B.est.cell.sim[,2], q.cell.1.sim$qvalue)
	#plot(sort(q.cell.2.sim$qvalue), fsr.cell.2$fdr, xlab="Est. Q-value", ylab="True FDR", type="l", main="FDR Plot WITHOUT Cell Type for Cov2")
	#plot(fsr.cell.2$fdr, fsr.cell.2$power, type="l")
	
	##Summarize Results
	plot(sort(q.cell.1.sim$qvalue), fsr.cell.2$fdr, xlab="Est. Q-value", ylab="True FDR", main=paste0("FDR Plot With and Without Cell Type, sigma.L = ", paste(as.character(round(100*sigma.L.sim)/100), collapse=",")), type="l")
	lines(sort(q.nocell.1.sim$qvalue), fsr.nocell.2$fdr, col="red")
	legend("bottomright", legend=c("With Cell", "Without Cell"), fill=c("black", "red"))
	abline(a=0,b=1, col="violet")
	
	return( list(Effects.cell=B.est.cell.sim, Effects.nocell=B.est.nocell.sim, fsr.cell=fsr.cell.2, fsr.nocell=fsr.nocell.2, q.cell=q.cell.1.sim, q.nocell=q.nocell.1.sim, cate.nocell=cate.nocell.sim, cate.cell=cate.cell.sim) )
}
