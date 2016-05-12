###This creates a randomized matrix and performs SVD.
###The goal here is to learn the number of additional confounders we have in this data

###The goal here is to get a crude estimate for the latent dimension r

Rand.SVD <- function(M) {
	p <- nrow(M)
	n <- ncol(M)
	M.rand <- M
	for (g in 1:p) {
		rand.index <- sample((1:n), n, replace=F)
		M.rand[g,] <- M[g,rand.index]
	}
	SMM.rand <- t(M.rand) %*% M.rand / p
	svd.rand <- svd(SMM.rand)
	return(svd.rand$d/sum(svd.rand$d))
}

Find.r <- function(M, n.ind, n.iter) {
	p <- nrow(M)
	n <- ncol(M)
	SMM <- t(M) %*% M / p
	svd.M <- svd(SMM)
	svd.d <- svd.M$d/sum(svd.M$d)
	svd.d <- svd.d[1:n.ind]
	p.values <- rep(0, n.ind)
	for (i in 1:n.iter) {
		rand.d <- Rand.SVD(M)[1:n.ind]
		p.values <- p.values + 1/n.iter * as.numeric(rand.d > svd.d)
	}
	return(p.values)
}