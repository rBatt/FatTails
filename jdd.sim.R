


# jdd.sim <- function(N, j.sigma=0.05, mu=10E-4, sigma=1E-2, lambda=10/350){
# 	# http://www.stat.berkeley.edu/~aldous/Research/Ugrad/ZY3.pdf
# 	S0 <- rnorm(1, 1, sd=0.05)
# 	S <- rep(NA, N)
# 	n.jump <- rpois(N, lambda=lambda)
# 	y.jump <- mapply(function(x)sum(rnorm(x, sd=j.sigma)), n.jump)
# 	# y.jump <- rnorm(N, sd=j.sigma)
# 	
# 	# dW <- sigma*rnorm(N)
# 	# W <- cumsum(dW)
# 	W <- cumsum(rnorm(N))*sigma
# 	
# 	dJ <- n.jump*y.jump
# 	J <- (1-y.jump)^cumsum(n.jump)
# 	
#  	# S0*exp((mu-0.5*sigma^2*(1:N)) + W)*J
#  	S0*exp((mu-0.5*sigma^2)*(1:N) + W*J)
# }



jdd.sim <- function(N, S0=0, dt=1, mu.star=0, sigma=0.1, lambda=0.1, mu.y=0, sigma.y=0.1){
	# N is the number of times steps
	# S0 is the starting price
	# dt is the size of the time step (less than or equal to 1)
	# mu.star is the interest rate
	# sigma is the standard deviation of the geometric brownian motion
	# lambda is the intensity of the jump process – i.e., expected number of jumps at each time step (well, dt*lambda)
	# mu.y is the mean of a jump
	# sigma.y is the standard deviation of a jump
	# http://arxiv.org/pdf/0812.4210.pdf
	m.sim <- c()
	for(i in 1:N){
		jumpnb <- rpois(1, lambda*dt)
		jump <- rnorm(1, (mu.y*(jumpnb-lambda*dt)), (sqrt(jumpnb)*sigma.y))
		m.sim[i] <- mu.star*dt + sigma*sqrt(dt)*rnorm(1) + jump
	}
	
	# cumsum(m.sim)+S0
	cumprod(c(S0,exp(m.sim*dt)))[-1]
}



# sample.size <- trunc(seq(2, 1E4, length.out=1E2))
# test.lognorm <- matrix(NA, nrow=length(sample.size), ncol=2, dimnames=list(NULL, c("max", "mean")))
# for(i in 1:length(sample.size)){t.sim <- rlnorm(sample.size[i]); test.lognorm[i-1,] <- c(max(t.sim), mean(t.sim))}
# 
# test.jdd <- matrix(NA, nrow=length(sample.size), ncol=2, dimnames=list(NULL, c("max", "mean")))
# for(i in 1:length(sample.size)){t.sim <- jdd.sim(sample.size[i]); test.jdd[i-1,] <- c(max(t.sim), mean(t.sim))}
# 
# dev.new()
# par(mfrow=c(3,2))
# plot(sample.size, log10(test.lognorm[,1]), xlab="sample size", ylab="log10 maximum")
# plot(sample.size, log10(test.jdd[,1]), xlab="sample size", ylab="log10 maximum")
# plot(sample.size, log10(test.lognorm[,2]), xlab="sample size", ylab="log10 mean")
# plot(sample.size, log10(test.jdd[,2]), xlab="sample size", ylab="log10 mean")
# plot(sample.size, cumsum(log10(test.lognorm[,1]))/(sample.size), xlab="sample size", ylab="log10 cummulative mean of maxima")
# plot(sample.size, cumsum(log10(test.jdd[,1]))/(sample.size), xlab="sample size", ylab="log10 cummulative mean of maxima")






