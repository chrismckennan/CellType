###This R document holds code that fit models that take into account cell type
#The methods may be unsupervised (i.e. case-control) or with cell type training data

##Case-Control##

#Fit L and Lambda using EM
#Yi is p x ni, L.0 is p x r, Lambda.0 is r x r, Sigma.0 is a p-vector
EM.CC <- function(Y1, Y2, L.0, Lambda.0, Sigma.0, max.iter, tol) {
	n1 <- ncol(Y1)
	n2 <- ncol(Y2)
	p <- nrow(Y1)
	r <- ncol(L)
	diag.y <- rowSums(Y1 * Y1) + rowSums(Y2 * Y2)
	
	SinvL.0 <- L.0 * 1/Sigma.0
	SinvY1 <- Y1 * 1/Sigma.0
	SinvY2 <- Y2 * 1/Sigma.0
	
	LtSinvL.0 <- t(L.0) %*% SinvL.0
	LtSinvY1.0 <- t(L.0) %*% SinvY1
	LtSinvY2.0 <- t(L.0) %*% SinvY2
	Y1tSinvY1 <- t(Y1) %*% SinvY1
	Y2tSinvY2 <- t(Y2) %*% SinvY2
	
	m.ident <- diag(r) + LtSinvL.0
	m.lambda <- solve(Lambda.0) + LtSinvL.0
	
	ll.0 <- -(n1 + n2)/2 * sum(log(Sigma.0)) - n1/2 * log(det( m.ident )) - n2/2 * log(det( m.lambda )) - n2/2 * log(det(Lambda.0)) - 1/2 * trace( Y1tSinvY1 + Y2tSinvY2 - t(LtSinvY1.0) %*% solve( m.ident, LtSinvY1.0 ) - t(LtSinvY2.0) %*% solve( m.lambda, LtSinvY2.0 ) )
	
	diff <- 1 + tol
	count <- 1
	
}


trace <- function(X) {
	return(sum(diag(X)))
}