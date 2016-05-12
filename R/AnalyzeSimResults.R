#This script looks at the true false discovery proportion, false sign proportion and power
false.sign.results <- function(beta.true, beta.est, qvalue.est) {
	p <- length(beta.true)
	n.nonzero <- max(length(which(beta.true != 0)), 1)
	false.sign.rate <- rep(0,p)
	fdr <- rep(0,p)
	power <- rep(0,p)
	order.q <- order(qvalue.est)
	qvalue.est <- sort(qvalue.est)
	beta.true <- beta.true[order.q]
	beta.est <- beta.est[order.q]
	for (i in 1:p) {
		b.true <- beta.true[i]
		b.est <- beta.est[i]
		if (b.true != 0) {
			sgn.true <- b.true/abs(b.true)
		} else {
			sgn.true <- 0
		}
		if (b.est != 0) {
			sgn.est <- b.est/abs(b.est)
		} else {
			sgn.est <- 0
		}		
		if (i == 1) {
			false.sign.rate[i] <- as.numeric( (sgn.true >= 0 && sgn.est < 0) || (sgn.true <= 0 && sgn.est > 0) )
			fdr[i] <- as.numeric( sgn.true == 0 && sgn.est != 0 )
			power[i] <- (1 - abs(abs(sgn.est) - abs(sgn.true)))/n.nonzero
		} else {		
			false.sign.rate[i] <- (i-1)/i * false.sign.rate[i-1] + 1/i * as.numeric( (sgn.true >= 0 && sgn.est < 0) || (sgn.true <= 0 && sgn.est > 0) )	
			fdr[i] <- (i-1)/i * fdr[i-1] + 1/i * as.numeric( sgn.true == 0 && sgn.est != 0 )
			power[i] <- power[i-1] + (1 - abs(abs(sgn.est) - abs(sgn.true)))/n.nonzero
		}
	}
	return(list(fsr=false.sign.rate, fdr=fdr, qorder=order.q, power=power))
}