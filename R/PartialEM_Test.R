p <- 100000
n <- 200
n.1 <- 150
r <- 3
K <- 2


Gamma.true <- matrix(rnorm(p*r), nrow=p, ncol=r)
F.mat.full <- matrix(rnorm(K*n), nrow=K, ncol=n)    #Lambda is the identity matrix
F.mat <- F.mat.full[,1:n.1]
F.mat <- t( scale( t(F.mat), center=T, scale=F ) )
L.true <- matrix(rnorm(p*K), nrow=p, ncol=K)/sqrt(n)

Y <- L.true %*% F.mat.full + Gamma.true %*% matrix(rnorm(r*n), nrow=r, ncol=n) + matrix(rnorm(n*p), nrow=p, ncol=n)
Y1 <- Y[,1:n.1]
Y2 <- Y[,(n.1+1):n]

Z1 <- t( scale(t(Y1), center=T, scale=F) )
Z2 <- t( scale(t(Y2), center=T, scale=F) )

m1 <- n.1 - 1
m2 <- (n - n.1) - 1

L.0 <- Z1 %*% t(F.mat) %*% solve(F.mat %*% t(F.mat))
Lambda.0 <- F.mat %*% t(F.mat) / m1
Sigma.0 <- rep(1, p)
update.Gamma <- T
Gamma.0 <- Gamma.true

#The partial maximum likelihood works for a GOOD STARTING POINT FOR GAMMA. With a poor starting point, we have little hope of ever reaching the maximum likelihood.