# RDB 30-June-2014
# Simulate Jump Drift Diffusion time series from 2 models
# Models may not be different, but I've coded them differently
# They come from 2 papers



# ==========================
# = Model 1: Johannes Eq 1 =
# ==========================
jdd.johan <- function(N, S0=1, dt=1, mu=0, sigma=0.1, lambda=0.1, sigma.y=0.1){
	# Programmed by me from Equation 1, 
	# I did not refer to any code or pseudo code
	
	# http://www0.gsb.columbia.edu/faculty/mjohannes/PDFpapers/JF_Diss.pdf

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



# ================================
# = Model 2: Brigo Matlab Code 6 =
# ================================
jdd.brigo <- function(N, S0=1, dt=1, mu.star=0, sigma=0.1, lambda=0.1, mu.y=0, sigma.y=0.1){
	# Translated from Matlab CODE 6 in Brigo et al. 2007 arxiv paper
	# Code is between Eq 30 and 31, but I did not refer to those equations 
	
	# http://arxiv.org/pdf/0812.4210.pdf
	
	# N is the number of times steps
	# S0 is the starting price
	# dt is the size of the time step (less than or equal to 1)
	# mu.star is the interest rate
	# sigma is the standard deviation of the geometric brownian motion
	# lambda is the intensity of the jump process – i.e., expected number of jumps at each time step (well, dt*lambda)
	# mu.y is the mean of a jump
	# sigma.y is the standard deviation of a jump
	
	m.sim <- c()
	for(i in 1:N){
		jumpnb <- rpois(1, lambda*dt)
		jump <- rnorm(1, (mu.y*(jumpnb-lambda*dt)), (sqrt(jumpnb)*sigma.y))
		m.sim[i] <- mu.star*dt + sigma*sqrt(dt)*rnorm(1) + jump
	}
	cumprod(c(S0,exp(m.sim*dt)))[-1]
}




