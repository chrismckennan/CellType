n <- 100
p <- 10000
r <- 1
dof.t <- 4.1
dof.t.w <- 2.1e0
Sigma <- rep(1, p)
alpha <- rep(1,r)
L <- matrix(rnorm(r * p), nrow=p, ncol=r) * sqrt(1/n)
L <- L %*% svd(t(L  / Sigma) %*% L)$v
X <- rbind(rep(1,n), rbinom(n,size=1,prob=1/2))
orthog.X <- diag(n) - t(X) %*% solve(X %*% t(X), X)
Q.X <- qr.Q(qr(t(X)), complete=T)[,(3:n)]
W <- matrix(rt(n*r, df=dof.t.w), nrow=r, ncol=n) %*% orthog.X * sqrt((dof.t.w - 2) / dof.t.w) + cbind(rep(0,r), alpha) %*% X
Y.t <-  L %*% W + matrix(rt(n*p, df=dof.t), nrow=p, ncol=n) * sqrt((dof.t - 2) / dof.t) * sqrt(1 / 1)
Y.norm <- L %*% W + matrix(rnorm(n*p), nrow=p, ncol=n) * sqrt(1 / 1)
tmp.fa.t <- fa.em(t(Y.t %*% Q.X), r=r)
tmp.fa.norm <- fa.em(t(Y.norm %*% Q.X), r=r)

info.norm <- svd(t(tmp.fa.norm$Gamma / tmp.fa.norm$Sigma) %*% tmp.fa.norm$Gamma / p)$d
V.norm <- svd(t(tmp.fa.norm$Gamma / tmp.fa.norm$Sigma) %*% tmp.fa.norm$Gamma / p)$v
V.t <- svd(t(tmp.fa.t$Gamma / tmp.fa.t$Sigma) %*% tmp.fa.t$Gamma / p)$v
info.t <- svd(t(tmp.fa.t$Gamma / tmp.fa.t$Sigma) %*% tmp.fa.t$Gamma / p)$d
alpha.norm <- t(V.norm) %*% solve(t(tmp.fa.norm$Gamma / tmp.fa.norm$Sigma) %*% tmp.fa.norm$Gamma, t(tmp.fa.norm$Gamma / tmp.fa.norm$Sigma) %*%Y.norm %*% t(X) %*% solve(X %*% t(X)))
alpha.t <- t(V.t) %*% solve(t(tmp.fa.t$Gamma / tmp.fa.t$Sigma) %*% tmp.fa.t$Gamma, t(tmp.fa.t$Gamma / tmp.fa.t$Sigma) %*%Y.t %*% t(X) %*% solve(X %*% t(X)))
alpha.t

alpha.true <- W %*% t(X) %*% solve(X %*% t(X))
norm.true <- (t(alpha.true) %*% alpha.true)[2,2]
norm.norm <- (t(alpha.norm) %*% alpha.norm)[2,2]
norm.t <- (t(alpha.t) %*% alpha.t)[2,2]
norm.t / norm.true


n <- 500
p <- 5e4
r <- 1
C <- 2     #Variances exist in the set [C^(-2), C^2]
Gamma <- rbind(matrix(rnorm(r * p * (1/n)), nrow=p*(1/n), ncol=r), matrix(0, nrow=(1 - 1/n)*p, ncol=r)) * sqrt(n)
Sigma <- rep(1, p)
Gamma <- Gamma %*% svd(t(Gamma / Sigma) %*% Gamma)$v
Gamma <- Gamma %*% diag(sqrt(c(1))/sqrt(svd(t(Gamma / Sigma) %*% Gamma / p)$d), nrow=r, ncol=r)
#Sigma <- runif(p) * (C^2 - 1/C^2) + 1/C^2
W <- matrix(rnorm(r*n), nrow=r, ncol=n) * sqrt(1/n); #W <- t(svd(t(W) %*% W)$v[,1:r])
#Gamma <- matrix( rnorm(p*r), nrow=p, ncol=r ) * sqrt(Sigma)
#Gamma <- Gamma %*% svd(t(Gamma / Sigma) %*% Gamma)$v
#Gamma[,1] <- Gamma[,1] * sqrt(2)
Y <- Gamma %*% W + matrix(rnorm(n*p), nrow=p, ncol=n) * sqrt(Sigma)
#Gamma.hat <- fa.em(t(Y), r=r)$Gamma
svd.Y <- svd(Y)
L.svd <- cbind(svd.Y$u[,1:r]) %*% diag(svd.Y$d[1:r], nrow=r, ncol=r) / sqrt(n)
sd.hat <- min( sd(sqrt(n) * (L.hat[1:(p/n),] + Gamma[1:(p/n),]/sqrt(n))), sd(sqrt(n) * (L.hat[1:(p/n),] - Gamma[1:(p/n),]/sqrt(n))) )
U.hat <- svd.Y$u[,(r+1):n]
mat.test <- t(L.hat / Sigma) %*% U.hat 
Sigma.hat <- fa.em(t(Y), r=r)$Sigma

W.hat <- t(svd.Y$v[,1:r])
W.hat %*% t(W)
W.hat %*% t(Y / Sigma) %*% Y %*% t(W.hat) / p
W.hat %*% t(Y / Sigma.hat) %*% Y %*% t(W.hat) / p

n <- 100
p <- 100000
Sigma <- runif(p) * (C^2 - 1/C^2) + 1/C^2
for (i in 1:10) {
	E <- matrix(rnorm(n*p), nrow=p, ncol=n) * sqrt(Sigma)
	R <- 1/p * t(E) %*% E - sum(Sigma)/p * diag(n)
	svd(R)$d[1]
}

Er <- matrix(rnorm(n), nrow=r, ncol=n)
svd(Er)$d[1]/sqrt(p)


#EM updates#
L.hat <- L.svd
for (i in 1:100) {
	YtSm1L <- t(Y / Sigma) %*% L.hat
	LtSm1L <- t(L.hat / Sigma) %*% L.hat
	G1 <- YtSm1L - YtSm1L %*% solve( diag(1, nrow=r, ncol=r) + LtSm1L, LtSm1L )
	G2 <- LtSm1L %*% solve(diag(1, nrow=r, ncol=r) + LtSm1L)
	L.hat <- (1/n * Y %*% G1) %*% solve( 1/n * t(G1) %*% G1 + diag(1, nrow=r, ncol=r) - G2 )
}




fa.em.Sigma <- function (Y, r, tol = 1e-08, maxiter = 1000, Sigma) {
    p <- ncol(Y)
    n <- nrow(Y)
    init <- fa.pc(Y, r)
    Gamma <- init$Gamma
    invSigma <- 1/Sigma
    llh <- -Inf
    I <- diag(rep(1, r))
    sample.var <- colMeans(Y^2)
    tilde.Gamma <- sqrt(invSigma) * Gamma
    M <- diag(r) + t(tilde.Gamma) %*% tilde.Gamma
    eigenM <- eigen(M, symmetric = TRUE)
    YSG <- Y %*% (invSigma * Gamma)
    logdetY <- -sum(log(invSigma)) + sum(log(eigenM$values))
    B <- 1/sqrt(eigenM$values) * t(eigenM$vectors) %*% t(YSG)
    logtrY <- sum(invSigma * sample.var) - sum(B^2)/n
    llh <- -logdetY - logtrY
    converged <- FALSE
    for (iter in 1:maxiter) {
        varZ <- eigenM$vectors %*% (1/eigenM$values * t(eigenM$vectors))
        EZ <- YSG %*% varZ
        EZZ <- n * varZ + t(EZ) %*% EZ
        eigenEZZ <- eigen(EZZ, symmetric = TRUE)
        YEZ <- t(Y) %*% EZ
        G <- sqrt(eigenEZZ$values) * t(eigenEZZ$vectors)
        #invSigma <- 1/(sample.var - 2/n * rowSums(YEZ * Gamma) + 
            #1/n * rowSums((Gamma %*% t(G))^2))
        Gamma <- YEZ %*% eigenEZZ$vectors %*% (1/eigenEZZ$values * 
            t(eigenEZZ$vectors))
        tilde.Gamma <- sqrt(invSigma) * Gamma
        M <- diag(r) + t(tilde.Gamma) %*% tilde.Gamma
        eigenM <- eigen(M, T)
        YSG <- Y %*% (invSigma * Gamma)
        old.llh <- llh
        logdetY <- -sum(log(invSigma)) + sum(log(eigenM$values))
        B <- 1/sqrt(eigenM$values) * t(eigenM$vectors) %*% t(YSG)
        logtrY <- sum(invSigma * sample.var) - sum(B^2)/n
        llh <- -logdetY - logtrY
        if (abs(llh - old.llh) < tol * abs(llh)) {
            converged <- TRUE
            break
        }
    }
    svd.H <- svd(t(Gamma) %*% (invSigma * Gamma))
    Z <- Y %*% (invSigma * Gamma) %*% (svd.H$u %*% (1/svd.H$d * 
        t(svd.H$v)))
    return(list(Gamma = Gamma, Sigma = 1/invSigma, Z = Z, niter = iter, 
        converged = converged))
}

p <- 1e5
n <- 500
E1 <- matrix(rnorm(K * p), nrow=p, ncol=K)
E2 <- matrix(rnorm(n*p), nrow=p, ncol=n)
Q <- qr.Q(qr(matrix(rnorm(n*K), nrow=n, ncol=K)))
test <- t(E1) %*% E2 %*% t(E2) %*% (E2 %*% Q)
test.1 <- t(E1) %*% E2 %*% Q %*% t(E2 %*% Q) %*% (E2 %*% Q)

p <- 1e4
n <- 1e3
k <- 1
E1 <- matrix(rnorm(k*p), nrow=k, ncol=p)
E2 <- matrix(rnorm(n*p), nrow=p, ncol=n)
E3 <- matrix(rnorm(p*k), nrow=p, ncol=k)
#Q <- qr.Q(qr(matrix(rnorm(n*k), nrow=n, ncol=k)))

1/p^2 * E1 %*% E2 %*% t(E2) %*% E3 * p/n

vec <- rep(NA, 20)
for (i in 1:length(vec)) {
	E1 <- matrix(rnorm(k*p), nrow=k, ncol=p)
	E2 <- matrix(rnorm(n*p), nrow=p, ncol=n)
	#1/p^2 * E1 %*% E2 %*% t(E2) %*% t(E1)
	vec[i] <- sqrt(n) * sqrt(sum((diag(1, nrow=k, ncol=k) - 1/p^2 * E1 %*% E2 %*% t(E2) %*% t(E1) * p/n)^2))/2
}


n <- 300
p <- 4e4
r <- 2
Gamma <- rbind(matrix(rnorm(r * floor(p * (1/n))), nrow=floor(p/n), ncol=r), matrix(0, nrow=(p - floor(p/n)), ncol=r)) * sqrt(n)
Sigma <- rep(1, p)
Gamma <- Gamma %*% svd(t(Gamma / Sigma) %*% Gamma)$v
Gamma <- Gamma %*% diag(sqrt(c(1))/sqrt(svd(t(Gamma / Sigma) %*% Gamma / p)$d), nrow=r, ncol=r)
W <- matrix(rnorm(r*n), nrow=r, ncol=n) * sqrt(1/n); W <- t(svd(t(W) %*% W)$v[,1:r])
Y <- Gamma %*% W + matrix(rnorm(n*p), nrow=p, ncol=n) * sqrt(Sigma)
svd.Y <- svd(Y)
L.svd <- cbind(svd.Y$u[,1:r]) %*% diag(svd.Y$d[1:r], nrow=r, ncol=r) / sqrt(n)
t(L.svd) %*% L.svd / p * n - diag(2, nrow=r, ncol=r)
W.hat <- svd.Y$v[,1:r]

E <- matrix(rnorm(n*p), nrow=p, ncol=n)
sd(svd(t(E) %*% E / p)$d)

#########Test of Partial SVD method#########
n <- 300
n.1 <- floor(5*n/10)
n.2 <- n - n.1
p <- 3.5e5
K <- 1

X <- cbind(rep(1, n), rbinom(n, size=1, prob=0.5))  #n x d
d <- ncol(X)
C <- matrix(rnorm(n*K), nrow=n, ncol=K) #n x K
Sigma <- rep(1, p)
L <- rbind(matrix(rnorm(K * floor(p * (1/n))), nrow=p*(1/n), ncol=K), matrix(0, nrow=p - floor(p * (1/n)), ncol=K)) * sqrt(2)
#L <- matrix(rnorm(p*K), nrow=p, ncol=K) * sqrt(2) * sqrt(1/n)
L <- L %*% svd(t(L / Sigma) %*% L)$v
ind.1 <- sort(sample((1:n), size=n.1, replace = F))

orthog.1 <- diag(n.1) - X[ind.1,] %*% solve(t(X[ind.1,]) %*% X[ind.1,], t(X[ind.1,]))
Q.X2 <- qr.Q(qr(X[-ind.1,]), complete=T)[,(d+1):(n-n.1)]
Q.X1 <- qr.Q(qr(X[ind.1,]), complete=T)[,(d+1):(n.1)]
orthog.x <- diag(n) - X %*% solve(t(X) %*% X, t(X))
tmp.C.1 <- C[ind.1,] - X[ind.1,] %*% solve(t(X[ind.1,]) %*% X[ind.1,], t(X[ind.1,]) %*% C[ind.1,])
tmp.C.2 <- C[-ind.1,] - X[-ind.1,] %*% solve(t(X[-ind.1,]) %*% X[-ind.1,], t(X[-ind.1,]) %*% C[-ind.1,])
tmp.svd <- svd(t(tmp.C.1) %*% tmp.C.1 / (n.1 - d)); tmp.Lambda <- tmp.svd$v %*% diag(sqrt(tmp.svd$d), nrow=K, ncol=K) %*% t(tmp.svd$v)
C <-  C %*% solve(tmp.Lambda)

M <- L %*% t(C) + matrix(rnorm(p*n), nrow=p, ncol=n) * sqrt(Sigma)
m <- n - 2*d
m.1 <- n.1 - d
m.2 <- n.2 - d
Y <- cbind(M[,ind.1] %*% Q.X1, M[,-ind.1] %*% Q.X2)
C1 <- t(Q.X1) %*% C[ind.1,]
C2 <- t(Q.X2) %*% C[-ind.1,]
C.rotate <- rbind(C1, C2)
svd.g <- svd(1/(n-2*d) * t(Y) %*% Y)
W <- t(svd.g$v[,1:K])
W1 <- rbind(W[,1:m.1])
W2 <- rbind(W[,(m.1+1):(n - 2*d)])
svd.W1 <- svd(W1)
B <- svd.W1$v
Lambda <- diag(svd.W1$d, nrow=K, ncol=K)
A <- svd.W1$u
A.tilde.0 <- svd(t(B) %*% C1)$v %*% t(svd(t(B) %*% C1)$u)
A.tilde <- A.tilde.0 %*% solve(Lambda, t(A))
C2.hat <- t(A.tilde %*% W2) * sqrt(m.2)
C1.hat <- t(A.tilde %*% W1) * sqrt(m.1)
diff.svd.hat <- sqrt(m) * (t(C2.hat) %*% C2.hat / m.2 - t(C2) %*% C2 / m.2)

A.tilde.1 <- t(C1) %*% t(W1) %*% solve(W1 %*% t(W1)) / sqrt(m)
A.tilde.all <- t(C.rotate) %*% t(W) / sqrt(m)
C2.hat.all <- t(A.tilde.all %*% W2) * sqrt(m)
C2.hat.1 <- t(A.tilde.1 %*% W2) * sqrt(m)
sqrt(m) * (t(C2.hat.1) %*% C2.hat.1 / m.2 - t(C2) %*% C2 / m.2)
sqrt(m) * (t(C2.hat.all) %*% C2.hat.all / m.2 - t(C2) %*% C2 / m.2)

