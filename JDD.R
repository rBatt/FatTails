# RDB 30-June-2014
# Simulate Jump Drift Diffusion time series from 2 models
# Models may not be different, but I've coded them differently
# They come from 2 papers

# Contents:
# 1) Line 15: Function to simulate JDD from Johannes 2004 The Journal of Finance
# 2) Line 45: Function to simulate JDD from Brigo 2007 arxiv.org (parts published in Jrnl Risk Mngmnt Fincl Insts, 2009)
# 3) Line 79: Function to fit GEV parameters from a vector of values (code based on Katz 2005 Ecology)
# 4) Line 200: Simulate JDD time series from both models
# 5) Line 226: Plot maxima of simulations, and regressions of xi vs. # obs per maximum



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








# =================================
# = Core function for fitting GEV =
# =================================	
# This is the function that I use to fit the GEV distribution

gev.fit <- function (xdat, ydat = NULL, mul = NULL, sigl = NULL, shl = NULL, mulink = identity, siglink = identity, shlink = identity, show = FALSE, method = "Nelder-Mead", maxit = 10000, ...){

    z <- list() #these are going to be aspects of the output, more later
    npmu <- length(mul) + 1 #number of parameters to be fitted for mu.  if mul is NULL, then length(mul) == 0, and npmu is 1
    npsc <- length(sigl) + 1 #number of params for sigma
    npsh <- length(shl) + 1 # n param for shape
    z$trans <- FALSE #transformation?

	Xi_init <- 0.1
	g1 <- gamma(1+Xi_init)
	g2 <- gamma(1+2*Xi_init)
	sigma_init <- sqrt((var(xdat)/(g2-g1^2))*Xi_init^2)
	mu_init <- mean(xdat) + sigma_init/Xi_init - sigma_init/(Xi_init*g1)
	
    # in2 <- sqrt(6 * var(xdat))/pi #the initial guess for sigma to be used in the optimization
    # in1 <- mean(xdat) - 0.57722 * in2 #initial guess for the location

    if(is.null(mul)){
        mumat <- as.matrix(rep(1, length(xdat))) #if mul is NULL, then only 1 parameter, the intercept, is needed for mu, and the predictor variable is 1 
        muinit <- mu_init
    }else{
        z$trans <- TRUE #if we do have non-NULL mul, then the transformation (transfer function?) is TRUE
        mumat <- cbind(rep(1, length(xdat)), ydat[, mul]) #the predictor variable matrix will include a column of 1's (for the intercept), and then all of the columns in ydat that correspond to other desired predictor variables.  E.g., mul could be equal to "Year", in which case our predictor matrix, mumat, would have a column of 1's of length xdat, and a column of the Years, also of length xdat.  The model would be similar to Y[i] = 1*B0 + Year[i]*B1, or Y = mumat%*%params, where Y is actually xdat, xdat has N rows and 1 column, mumat has N rows and 2 columns, and params is a matrix with 2 rows and 1 column
        muinit <- c(mu_init, rep(0, length(mul))) #now we have to make more guesses for the initial values --- we calculate the intercept as the "default" estimate of the parameter, and 0 (or no effect) as the default estimate of the other parameters
    }
    if (is.null(sigl)){ #same idea as for mul, and subsequent estiamtes of initials and matrices
        sigmat <- as.matrix(rep(1, length(xdat)))
        siginit <- sigma_init
    }else{
        z$trans <- TRUE
        sigmat <- cbind(rep(1, length(xdat)), ydat[, sigl])
        siginit <- c(sigma_init, rep(0, length(sigl)))
    }

    if (is.null(shl)){ #same idea as for mu and sigma
        shmat <- as.matrix(rep(1, length(xdat)))
        shinit <- 0.1
    }else{
        z$trans <- TRUE
        shmat <- cbind(rep(1, length(xdat)), ydat[, shl])
        shinit <- c(0.1, rep(0, length(shl)))
    }

    z$model <- list(mul, sigl, shl) #this is a list of 3 elements, where each element is a character vector whose elements in turn specify a column (name) of ydat to be used as a predictor variables for the respective GEV parameter
    z$link <- deparse(substitute(c(mulink, siglink, shlink))) #turns the link functions used into character strings
    init <- c(muinit, siginit, shinit) #Consider log-transforming the siginit, and then inside gev.like exponentiating. This would avoid the problem of sigma being <= 0.
		#however, problem with transforming is that some of the sigma parameters can be negative, so long as none of the sc are negative. 
		#essentially, the sum of some parameters related to sigma can't be negative, but any one of them *may* be allowed negative
		#if we leave the 10^6 solution as-is, the optimization can get "lost" on a massively-elevated NLL plane
		#perhaps a solution would be to slope the plane such that the NLL gets worse and worse as sc gets further below 0. same idea for Y.  
		#an initial jump in NLL when either sc or y reach 0, and then a steady increase afterwards.  maybe I'll try this.

    gev.lik <- function(a){
        mu <- mulink(mumat %*% (a[1:npmu]))#unpack the mu parameters from the a vector
        sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)])) #unpack the sigma parameters from the a vector
        xi <- shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)])) #unpack the shape parameters from the a vector
        y <- (xdat - mu)/sc #this line is part of t(x)
        y <- 1 + xi * y #this line is more of t(x), and is now only missing the ^(-1/Xi) part
        if (any(y <= 0) || any(sc <= 0)){
			if(any(y <= 0)){ypenal <- (1/exp(y[which(y<=0)]))}else{ypenal <- 0}
			if(any(sc <= 0)){scpenal <- (1/exp(sc[which(sc<=0)]))}else{scpenal <- 0}
			BadReturn <- 10^5 * (mean(ypenal) + mean(scpenal))
			return(BadReturn)
            # return(10^6) #if you get an impossible value for y or sc, just return a huge nll, which will discourage the optimizer from picking the current parameters
			#also, this is a shortcut for the support that X is an element of (mu-sigma/Xi to INF) when Xi >0, and X is an element of (-INF to mu-sigma/Xi) when Xi < 0
				#there is only support for those values of X because t_x must be defined... 
				#in the case where Xi!=0, (1+((x-mu)/sigma)*xi), which is defined as "y" above, ends up being in the denominator when Y^(-1/xi)
				# Y^(-1/xi) is t(x)
				#while Y can't be zero, b/c denom, it can't be negative either b/c we take the log of y in the likelihood equation
		}
        sum(log(sc)) + sum(y^(-1/xi)) + sum(log(y)*(1/xi + 1)) #negative log likelihood (from the pdf of the gev)
			#sum(log(sc)) is the -log() of 1/sigma
			#sum(y^(-1/xi)) is the -log() of exp(-t(x))
			#sum(log(y) * (1/xi + 1)) is the -log() of t(x)^(xi+1) 
				#because had (y^-1/xi)^*(xi+1)  
				#==> y^(-1 + -1/xi) 
				#==> log(y) * -(1 + 1/xi) 
				#==> -log(y)*(1/xi + 1) (and then take neg for nll)
    }

    x <- optim(init, gev.lik, hessian = TRUE, method = method, control = list(maxit = maxit, ...))
    z$conv <- x$convergence
    mu <- mulink(mumat %*% (x$par[1:npmu]))
    sc <- siglink(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))
    xi <- shlink(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))
    z$nllh <- x$value
    z$data <- xdat

    if (z$trans){
        z$data <- -log(as.vector((1 + (xi * (xdat - mu))/sc)^(-1/xi)))
    }
    z$mle <- x$par
	names(z$mle) <- c(paste("mu",c(0,mul),sep="_"), paste("sig",c(0,sigl),sep="_"), paste("sh",c(0,shl),sep="_"))
    z$cov <- solve(x$hessian)#I should use the hessian() function in package numDeriv to calculate this, and compare it to the hessian from the optimization algorithm.
    z$se <- sqrt(diag(z$cov))
	names(z$se) <- c(paste("mu",c(0,mul),sep="_"), paste("sig",c(0,sigl),sep="_"), paste("sh",c(0,shl),sep="_"))
    z$vals <- cbind(mu, sc, xi)
    if(show){
        if(z$trans){
            print(z[c(2, 3, 4)])
        }else{print(z[4])}
        if(!z$conv){
            print(z[c(5, 7, 9)])
		}
    }
    invisible(z)
}









# ==============================
# = Simulations of both models =
# ==============================
n.maxima <- 50
n.obs.per.block <- seq(10, 1E3, length.out=10) # I did 199 originally

# Big simulation involving 50 simulations of each of 199 time series that individually range between 10 and 1E3 time steps
jdd.johan.sims <- mapply(function(x)apply(mapply(jdd.johan, rep(x, n.maxima)), 2, max), n.obs.per.block)
jdd.johan.xi <- apply(jdd.johan.sims, 2, function(x)gev.fit(x)$mle["sh_0"])
jdd.johan.xi.log <- apply(jdd.johan.sims, 2, function(x)gev.fit(log(x))$mle["sh_0"])

jdd.brigo.sims <- mapply(function(x)apply(mapply(jdd.brigo, rep(x, n.maxima)), 2, max), n.obs.per.block)
jdd.brigo.xi <- apply(jdd.brigo.sims, 2, function(x)gev.fit(x)$mle["sh_0"])
jdd.brigo.xi.log <- apply(jdd.brigo.sims, 2, function(x)gev.fit(log(x))$mle["sh_0"])

# Small simulation of a single time series (50 years of 500 observations per year)
set.seed(2)
jdd.johan.1ts <- jdd.johan(N=500)

set.seed(4)
jdd.brigo.1ts <- jdd.brigo(N=500)





# ===============================
# = Plotting Simulation Results =
# ===============================
# Create hot to cold color scheme
ramp.cols <- col2rgb(colorRampPalette(c("blue", "cyan", "green", "yellow", "orange", "red"))(length(n.obs.per.block)), alpha=FALSE) # create color ramp
mycols <- rgb(red=ramp.cols["red",], green=ramp.cols["green",], blue=ramp.cols["blue",], alpha=50, maxColorValue=256) # add transparency

dev.new(width=5.5, height=8.5) # new graphical device
# png("~/Desktop/JDD_2mods.png", width=5.5, height=8.5, units="in", res=150) # use this line instead of dev.new() if you want to save the png
par(mfcol=c(4, 2), mar=c(2.25, 3, 0.5, 1.5), oma=c(0.1, 0.1, 10, 0.1), ps=10, mgp=c(0.9, 0.15, 0), tcl=-0.20, family="Times", cex=1) # graphical params

# Plotting Johannes 2004 Eq 1
plot(jdd.johan.1ts, type="l", xlab="time (within 1 year)", ylab="Returns (?) (r)")
abline(h=0)

# Add Johannes equation
mtext(bquote( # begin mtext of equations
	atop({
				{r[t+Delta]-r[t]} == {mu*(r[t])*Delta+sigma*(r[t])*epsilon[t+Delta]*sqrt(Delta)+J[t+Delta]*Z[t+Delta]}
			},
			{
				{J[t+Delta]*Z[t+Delta]} == {sum({N*(0*','*sigma[z])[i]}, i==1, {Pois*(lambda*Delta)})}
			}
	))
, line=3.5)
mtext("Model 1\n(Eq 1 Johannes 2004 )", line=8, font=2) # mtext for model label

# Plot time series of max returns
jdd.johan.ylim <- log(range(jdd.johan.sims)) # calculate ylims for plot of maxima time series
plot(log(jdd.johan.sims[,length(n.obs.per.block)]), type="l", ylim=jdd.johan.ylim, col=mycols[1], xlab="time (across 50 years)", ylab="log Max Returns (r)") # plot first maxima
for(i in (length(n.obs.per.block)-1):1){ # loop through other maxima, plotting with lines()
	lines(log(jdd.johan.sims[,i]), type="l", col=mycols[i])
}
legend("topright", legend=c("N=10", "N=1000"), text.col=c("blue", "red"), cex=1, bty="n") # add color scheme legend


# scatter of xi vs. number of observations used per each maximum (x axis is analagous to # obs per year when calculating xi from annual maxima)
plot(n.obs.per.block, jdd.johan.xi, xlab="N Obs per Maximum", ylab=bquote(xi~~of~~r~~Maxima))
abline(lm(jdd.johan.xi~n.obs.per.block), lwd=3) # add solid black regression line
abline(lm(jdd.johan.xi~n.obs.per.block), lwd=1, col="white")
abline(h=median(jdd.johan.xi), lty="dashed", col="blue", lwd=4) # add dotted blue median line
abline(h=median(jdd.johan.xi), lty="dashed", col="lightblue", lwd=2) # add dotted blue median line

# same as above, but for log(maxima)
plot(n.obs.per.block, jdd.johan.xi.log, xlab="N Obs per Maximum", ylab=bquote(xi~~of~~log[e]*(r)~~Maxima))
abline(lm(jdd.johan.xi.log~n.obs.per.block), lwd=3)
abline(lm(jdd.johan.xi.log~n.obs.per.block), lwd=1, col="white")
abline(h=median(jdd.johan.xi.log), lty="dashed", col="blue", lwd=4)
abline(h=median(jdd.johan.xi.log), lty="dashed", col="lightblue", lwd=2)


# Plotting Brigo 2007 Code 6
# Most of this code is same as above – only commenting on differences
# Plotting Johannes 2004 Eq 1
plot(jdd.brigo.1ts, type="l", xlab="time (within 1 year)", ylab="Price (?) (S)", ylim=c(-0.1, max(jdd.brigo.1ts)))
abline(h=0)

# Add Brigo Equation
mtext(bquote(
	atop(
		{
			{r[t+Delta]-r[t]} == {mu*(r[t])*Delta+sigma*(r[t])*epsilon[t+Delta]*sqrt(Delta)+J[t+Delta]*Z[t+Delta]}
		},
		{
			{J[t+Delta]*Z[t+Delta]} == {sum({N*(0*','*sigma[z]*sqrt(n*(t+Delta)))[i]}, i==1, {n*(t+Delta)})} # note that this equation is slightly different from previous (multiply sigma[z] by sqrt(n(t+Delta)))
		}
		)
	)
, line=3.5)
mtext(bquote( # only difference in plotting is that I need 2 more equations, which are inserted here
	atop(
		{
			{n*(t+Delta)}=={Pois*(lambda*Delta)}
		},
		{
			{S*(T)} == {prod({exp*(r[t]*Delta)}, t==1, T)}
		}
		)
	)
, line=0)
mtext("Model 2\n(Code 6 from Brigo 2007 in arxiv.org)", line=8, font=2)

# Plot max prices
jdd.brigo.ylim <- log(range(jdd.brigo.sims))
plot(log(jdd.brigo.sims[,length(n.obs.per.block)]), type="l", ylim=jdd.brigo.ylim, col=mycols[1], xlab="time (across 50 years)", ylab="log Max Price (S)") # I think that this is the Price, and others are Returns
for(i in (length(n.obs.per.block)-1):1){
	lines(log(jdd.brigo.sims[,i]), type="l", col=mycols[i])
}
legend("topright", legend=c("N=10", "N=1000"), text.col=c("blue", "red"), cex=1, bty="n")

# Scatter plots
plot(n.obs.per.block, jdd.brigo.xi, xlab="N Obs per Maximum", ylab=bquote(xi~~of~~S~~Maxima))
abline(lm(jdd.brigo.xi~n.obs.per.block), lwd=3)
abline(lm(jdd.brigo.xi~n.obs.per.block), lwd=1, col="white")
abline(h=median(jdd.brigo.xi), lty="dashed", col="blue", lwd=4)
abline(h=median(jdd.brigo.xi), lty="dashed", col="lightblue", lwd=2)

plot(n.obs.per.block, jdd.brigo.xi.log, xlab="N Obs per Maximum", ylab=bquote(xi~~of~~log[e]*(S)~~Maxima)) 
abline(lm(jdd.brigo.xi.log~n.obs.per.block), lwd=3)
abline(lm(jdd.brigo.xi.log~n.obs.per.block), lwd=1, col="white")
abline(h=median(jdd.brigo.xi.log), lty="dashed", col="blue", lwd=4)
abline(h=median(jdd.brigo.xi.log), lty="dashed", col="lightblue", lwd=2)
# dev.off() # use this if you want to save the png

