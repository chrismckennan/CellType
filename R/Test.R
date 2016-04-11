p <- nrow(beta.values)
n <- ncol(beta.values)
model.mat <- model.matrix(~sex + as.factor(chip) + as.factor(Treatment) + Age + asthma + DNA_concE1 + DNA_concE2, data=data.cov)
P.orthog <- diag(n) - model.mat %*% solve(t(model.mat) %*% model.mat, t(model.mat))
M <- log(beta.values/(1-beta.values))
M.center <- M %*% P.orthog
SMM <- t(M.center) %*% M.center / p
svd.SMM <- svd(SMM)
eigs <- svd.SMM$d
v <- svd.SMM$v
u <- M.center %*% v[,1:40]
plot(eigs, pch=".")

corr.pvalues.1 <- rep(NA, p)
for (i in 1:p) {
	y.i <- M.center[i,]
	u.i <- u[i,1]
	var <- 1/(n-1) * (sum(y.i * y.i) - u.i^2)
	t.i <- abs(u.i/sqrt(var))
	corr.pvalues.1[i] <- 2*pt(-t.i, df=n-1)
}
hist(corr.pvalues.1, breaks=90)
pi.0 <- 0.1/0.9 * length(corr.pvalues.1[corr.pvalues.1 > 0.2])/p

k <- 30
n <- 400
n.1 <- 60
n.2 <- n - n.1
p <- 100000
W <- matrix(0.1*rnorm(k*n), nrow=k, ncol=n)
W.1 <- W[,1:n.1]
W.2 <- W[,(n.1+1):n]
D <- matrix(rnorm(p*k), nrow=p, ncol=k)
E <- D %*% W
E.1 <- E[,1:n.1]
E.2 <- E[,(n.1+1):n]

SEE.1 <- t(E.1) %*% E.1
svd.SEE.1 <- svd(SEE.1)
v.1 <- svd.SEE.1$v[,1:k]
u.1 <- E.1 %*% v.1

SEE.2 <- t(E.2) %*% E.2
svd.SEE.2 <- svd(SEE.2)
v.2 <- svd.SEE.2$v[,1:k]
u.2 <- E.2 %*% v.2

guess.v2 <- solve(t(u.1)%*%u.1, t(u.1)) %*% E.2
proj.guess.v2 <- t(guess.v2) %*% solve(guess.v2%*%t(guess.v2)) %*% guess.v2
proj.W2 <- t(W.2) %*% solve(W.2 %*% t(W.2)) %*% W.2

norm(proj.W2 - proj.guess.v2, "2")


X <- rdirichlet(50, 15*c(0.2, 0.5, 0.3))
H <- X %*% solve(t(X)%*%X, t(X))
beta <- c(0.01, 0.02, 0.03)
Y <- X %*% beta + 0.001*rnorm(50)
Y[Y < 0] <- -Y[Y<0]

Y.mean <- Y - mean(Y)
SSY <- sum(Y.mean^2)
SSreg <- sum((H%*%Y.mean)^2)
R2 <- SSreg/SSY

Y.logit <- log(Y/(1-Y))
Y.logit.mean <- Y.logit - mean(Y.logit)
SSY.logit <- sum(Y.logit.mean^2)
SSreg.logit <- sum((H%*%Y.logit.mean)^2)
R2.logit <- SSreg.logit/SSY.logit






