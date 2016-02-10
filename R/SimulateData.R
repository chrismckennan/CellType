###This script contains functions for various RMD files, including SimulateCellTypeConfounding.Rmd,
###160209

#Sample elements of rows of a matrix independently (to be used in bootstrap procedures to esimate the latent dimension)

PermMat <- function(X) {
	p <- nrow(X)
	n <- ncol(X)
	X.out <- array(0, dim=c(p,n))
	for (r in 1:p) {
		X.out[r,] <- sample(X[p,], n, replace=T)
	}
	return(X.out)
}

#Take the logit transform of all individual elements in a matrix
logit <- function(X) {
	return( log(X / (1 - X)) )
}