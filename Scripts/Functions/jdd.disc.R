
jdd.disc <- function(N, S0=1, dt=1, mu=0, sigma=0.1, lambda=0.1, sigma.y=0.1){
	# http://www0.gsb.columbia.edu/faculty/mjohannes/PDFpapers/JF_Diss.pdf
	# Equation 1
	
	# N is the number of times steps
	# S0 is the starting price
	# dt is the size of the time step (less than or equal to 1)
	# mu is the interest rate (drift)
	# sigma is the standard deviation of the geometric brownian motion
	# lambda is the intensity of the jump process – i.e., expected number of jumps at each time step (well, dt*lambda)
	# sigma.y is the standard deviation of a jump
	
	r <- c(S0, rep(NA, N-1))
	
	for(i in 2:(N+1)){
		n.jump <- rpois(1, lambda=lambda*dt)
		r[i] <- r[i-1] + mu*r[i-1]*dt + sigma*r[i-1]*rnorm(1)*sqrt(dt) + sum(rnorm(n.jump, sd=sigma.y))
	}
	
	r[-1]

}


n.maxima <- 50
n.obs.per.block <- 11:100
jd.sims <- mapply(function(x)apply(mapply(jdd.disc, rep(x, n.maxima)), 2, max), n.obs.per.block)
jd.xi <- apply(jd.sims, 2, function(x)gev.fit(x)$mle["sh_0"])
plot(n.obs.per.block, jd.xi); abline(lm(jd.xi~n.obs.per.block))
