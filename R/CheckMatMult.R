####Check matrix multiplications in EM####
n.1 <- 10
n.2 <- 10
r <- 5
K <- 3
p <- 100

L.0 <- matrix(rnorm(K*p), nrow=p, ncol=K)
Gamma.0 <- matrix(rnorm(r*p), nrow=p, ncol=r)
Y1 <- matrix(rnorm(n.1*p), nrow=p, ncol=n.1)
Y2 <- matrix(rnorm(n.2*p), nrow=p, ncol=n.2)
F.mat <- matrix(rnorm(n.1*K), nrow=K, ncol=n.1)
Sigma.0 <- abs(rnorm(p))
tmp.S <- matrix(rnorm(K^2), nrow=K, ncol=K)
Lambda.0 <- tmp.S %*% t(tmp.S)
SFF <- F.mat %*% t(F.mat)    #A k x k matrix
Y1Ft <- Y1 %*% t(F.mat)      #A p x K matrix

S.chol.0 <- t(chol(Lambda.0))     #SS' = Lambda.0
		
SinvL.0 <- sweep(L.0, 1, 1/Sigma.0, "*")    #Sigma.0^{-1} L.0
SinvG.0 <- sweep(Gamma.0, 1, 1/Sigma.0, "*")    #Sigma.0^{-1} Gamma.0
		
LtSinvL.0 <- t(L.0) %*% SinvL.0      #L.0' Sigma.0^{-1} L.0, a k x k matrix
LtSinvG.0 <- t(SinvL.0) %*% Gamma.0    #L' Sigma.0^{-1} Gamma.0
GtSinvG.0 <- t(SinvG.0) %*% Gamma.0  #Gamma.0' Sigma.0^{-1} Gamma.0
		
LtSinvY2.0 <- t(SinvL.0) %*% Y2     #L.0' Sigma.0^{-1} Y2
GtSinvY2.0 <- t(SinvG.0) %*% Y2     #Gamma.0' Sigma.0^{-1} Y2

Y1mLF.0 <- Y1 - L.0 %*% F.mat      #Y1 - L.0 F, a p x n.1 matrix. This will also be referred to as 'Z' in the below code
		
GtSinvZ.0 <- t(SinvG.0) %*% Y1mLF.0        #Gamma.0' Sigma^{-1} (Y1 - L.0 F)
		
Middle.0 <- diag(r) + GtSinvG.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 )       #See notes

H <- solve(diag(Sigma.0) + L.0 %*% Lambda.0 %*% t(L.0) + Gamma.0 %*% t(Gamma.0))
G <- solve(diag(Sigma.0) + Gamma.0 %*% t(Gamma.0))

##L.0' H L.0, a K x K matrix##
LtHL.0 <- ( LtSinvL.0 - LtSinvL.0 %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvL.0 ) ) - ( LtSinvG.0 - LtSinvL.0 %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ) ) %*% solve( diag(r) + GtSinvG.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ), t( LtSinvG.0 - LtSinvL.0 %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ) ) )         #L.0' H L.0, a K x K matrix

LtHL.true <- t(L.0) %*% H %*% L.0

#This works

##L.0' H Y2, a K x n.2 matrix##
LtHY2.0 <- ( LtSinvY2.0 - LtSinvL.0 %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvY2.0 ) ) - ( LtSinvG.0 - LtSinvL.0 %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ) ) %*% solve( diag(r) + GtSinvG.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ), GtSinvY2.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvY2.0 ) ) 

LtHY2.true <- t(L.0) %*% H %*% Y2

#This works

##Gamma.0' G (Y1 - L.0 F)##
GtGZ.0 <- GtSinvZ.0 - GtSinvG.0 %*% solve( diag(r) + GtSinvG.0, GtSinvZ.0 )

GtGZ.true <- t(Gamma.0) %*% G %*% (Y1 - L.0 %*% F.mat)

#This works

###Gamma.0' H Y2, a r x n.2 matrix###
GtHY2.0 <- ( GtSinvY2.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvY2.0 ) ) - ( GtSinvG.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ) ) %*% solve( Middle.0, GtSinvY2.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvY2.0 ) )

GtHY2.true <- t(Gamma.0) %*% H %*% Y2

#This works

###Gamma.0' H L, a r x K matrix###
GtHL.0 <- ( t(LtSinvG.0) - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvL.0 ) ) - ( GtSinvG.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0 ) ) %*% solve( Middle.0, t(LtSinvG.0) - t(LtSinvG.0) %*% S.chol.0 %*% solve( diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvL.0 ) )

GtHL.true <- t(Gamma.0) %*% H %*% L.0


###Gamma.0' G Gamma.0###
GtGG.0 <- GtSinvG.0 - GtSinvG.0 %*% solve( diag(r) + GtSinvG.0, GtSinvG.0 )

GtGG.true <- t(Gamma.0) %*% G %*% Gamma.0

#This works

###Gamma.0' H Gamma.0###
A.GtHG.0 <- GtSinvG.0 - t(LtSinvG.0) %*% S.chol.0 %*% solve(diag(K) + t(S.chol.0) %*% LtSinvL.0 %*% S.chol.0, t(S.chol.0) %*% LtSinvG.0)       #See notes
GtHG.0 <- A.GtHG.0 - A.GtHG.0 %*% solve(Middle.0, A.GtHG.0)

GtHG.true <- t(Gamma.0) %*% H %*% Gamma.0

#This works

##All Woodbury identities are coded correctly


####Check to see if EM algorithm is working correctly via simulation###

n.1 <- 100
n.2 <- 100
p <- 1000
K <- 3
r <- 5
v <- 0.1^2

L <- matrix(0, nrow=p, ncol=K)
tmp.S <- matrix(rnorm(K^2), nrow=K, ncol=K)*1/10
Lambda <- tmp.S %*% t(tmp.S)
Lambda <- diag(K)
F.mat <- tmp.S %*% matrix( rnorm(n.1*K), nrow=K, ncol=n.1 )
Gamma <- matrix(rnorm(r*p), nrow=p, ncol=r)
Sigma <- v * rep(1,p)

Y1 <- L %*% F.mat + Gamma %*% matrix( rnorm(r*n.1), nrow=r, ncol=n.1 ) + sqrt(v) * matrix(rnorm(n.1*p), nrow=p, ncol=n.1)
Y2 <- L %*% tmp.S %*% matrix( rnorm(n.2*K), nrow=K, ncol=n.2 ) + Gamma %*% matrix( rnorm(r*n.2), nrow=r, ncol=n.2 ) + sqrt(v) * matrix(rnorm(n.2*p), nrow=p, ncol=n.2)

L.0 <- L
Gamma.0 <- Gamma + sqrt(0.5)*matrix(rnorm(p*r), nrow=p, ncol=r)
Sigma.0 <- Sigma
Lambda.0 <- Lambda


#It looks like we are overestimating Gamma by a factor of the noise of the starting value...



##Check update for Gamma##
#If L = 0, update is Gamma = ((Z1Z1' + Z2Z2') G Gamma.0) ((n.1 + n.2)(I_r - Gamma.0' G Gamma.0) + Gamma.0' G (Z1Z1' + Z2Z2') G Gamma.0)^{-1}
G <- solve(diag(Sigma.0) + Gamma.0 %*% t(Gamma.0))
Update.true <- ( ( Y1 %*% t(Y1) + Y2 %*% t(Y2) ) %*% G %*% Gamma.0 ) %*% solve( (n.1 + n.2)*(diag(r) - t(Gamma.0) %*% G %*% Gamma.0) + t(Gamma.0) %*% G %*% (Y1 %*% t(Y1) + Y2 %*% t(Y2)) %*% G %*% Gamma.0 )

#Updates are correct. I think I need to spend more time updating Gamma, since Sigma and L appear to converge quickly. Also the initialized values should be good (or at least the values from CATE).

##Test my EM code with CATE's starting value for Gamma
n.1 <- 300
p <- 1000
r <- 5
v <- 0.1^2
tol <- 1e-6


Gamma.start <- matrix(rnorm(r*p), nrow=p, ncol=r)
Sigma <- v * rep(1,p)
Y1 <- Gamma.start %*% matrix( rnorm(r*n.1), nrow=r, ncol=n.1 ) + sqrt(v) * matrix(rnorm(n.1*p), nrow=p, ncol=n.1)
Y <- Y1 - apply(Y1, 1, mean)




