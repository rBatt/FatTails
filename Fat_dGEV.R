#Ryan Batt
#24-April-2013
#Modified from Katz et al. 2005, Coles 2001
	#changed the way that initial values were guessed (they had used the equations for Xi=0, I used equations for Xi < 1 b/c our guess for Xi is 0.1)
	#I changed the likelihood penalty for having y<=0 or sc<=0; originally 10^6.  Changed to increase the penalty as y or sc became more negative.
	# Note to self: use model.matrix(), to get the ydat matrix --- very helpful for turning a column of factors (w/ n levels) into a matrix with n columns that contains 0 or 1
# _v0.2 (02-May-2013): Added names to the estimates in the output
	
gev.fit <- function (xdat, ydat = NULL, mul = NULL, sigl = NULL, shl = NULL, 

    mulink = identity, siglink = identity, shlink = identity, 

    show = TRUE, method = "Nelder-Mead", maxit = 10000, ...) 

{

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

    if (is.null(mul)) {

        mumat <- as.matrix(rep(1, length(xdat))) #if mul is NULL, then only 1 parameter, the intercept, is needed for mu, and the predictor variable is 1 

        muinit <- mu_init

    }

    else {

        z$trans <- TRUE #if we do have non-NULL mul, then the transformation (transfer function?) is TRUE

        mumat <- cbind(rep(1, length(xdat)), ydat[, mul]) #the predictor variable matrix will include a column of 1's (for the intercept), and then all of the columns in ydat that correspond to other desired predictor variables.  E.g., mul could be equal to "Year", in which case our predictor matrix, mumat, would have a column of 1's of length xdat, and a column of the Years, also of length xdat.  The model would be similar to Y[i] = 1*B0 + Year[i]*B1, or Y = mumat%*%params, where Y is actually xdat, xdat has N rows and 1 column, mumat has N rows and 2 columns, and params is a matrix with 2 rows and 1 column

        muinit <- c(mu_init, rep(0, length(mul))) #now we have to make more guesses for the initial values --- we calculate the intercept as the "default" estimate of the parameter, and 0 (or no effect) as the default estimate of the other parameters

    }

    if (is.null(sigl)) { #same idea as for mul, and subsequent estiamtes of initials and matrices

        sigmat <- as.matrix(rep(1, length(xdat)))

        siginit <- sigma_init

    }

    else {

        z$trans <- TRUE

        sigmat <- cbind(rep(1, length(xdat)), ydat[, sigl])

        siginit <- c(sigma_init, rep(0, length(sigl)))

    }

    if (is.null(shl)) { #same idea as for mu and sigma

        shmat <- as.matrix(rep(1, length(xdat)))

        shinit <- 0.1

    }

    else {

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

    x <- optim(init, gev.lik, hessian = TRUE, method = method, 

        control = list(maxit = maxit, ...))

    z$conv <- x$convergence

    mu <- mulink(mumat %*% (x$par[1:npmu]))

    sc <- siglink(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))

    xi <- shlink(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))

    z$nllh <- x$value

    z$data <- xdat

    if (z$trans) {

        z$data <- -log(as.vector((1 + (xi * (xdat - mu))/sc)^(-1/xi)))

    }

    z$mle <- x$par
	names(z$mle) <- c(paste("mu",c(0,mul),sep="_"), paste("sig",c(0,sigl),sep="_"), paste("sh",c(0,shl),sep="_"))

    z$cov <- solve(x$hessian)#I should use the hessian() function in package numDeriv to calculate this, and compare it to the hessian from the optimization algorithm.

    z$se <- sqrt(diag(z$cov))
	names(z$se) <- c(paste("mu",c(0,mul),sep="_"), paste("sig",c(0,sigl),sep="_"), paste("sh",c(0,shl),sep="_"))

    z$vals <- cbind(mu, sc, xi)

    if (show) {

        if (z$trans) 

            print(z[c(2, 3, 4)])

        else print(z[4])

        if (!z$conv) 

            print(z[c(5, 7, 9)])

    }

    invisible(z)

}


