
# =================================
# = Core function for fitting GEV =
# =================================	
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




# ================================
# = Method of independent storms =
# ================================
MIS <- function(x, Thresh=NULL, quant=0.9){
	if(is.null(Thresh)){
		Thresh <- quantile(x, quant)
	}
	Above <- as.integer(x>=Thresh)
	udF <- ifelse(Above[1]==1, 1, 0)
	udL <- ifelse(rev(Above)[1]==1, -1, 0)
	UpDown <- c(udF,diff(Above),udL)
	StartStop <- matrix(sort(c(which(UpDown==1), (which(UpDown==-1)-1))), ncol=2, byrow=TRUE)
	Storms <- rep(NA, nrow(StartStop))
	StormIndex <- Storms
	for(i in 1:nrow(StartStop)){
		ThisInd <- StartStop[i,1]:StartStop[i,2]
		StormIndex[i] <- which.max(x[ThisInd]) + ThisInd[1] - 1
		Storms[i] <- x[StormIndex[i]]

	}
	return(matrix(c(StormIndex, Storms), ncol=2, dimnames=list(NULL, c("Index", "Max"))))
}

# ========================================
# = Find the peak of a cycle via wavelet =
# ========================================
PeakCycle <- function(Data, SearchFrac=0.028){
	if(!library("wmtsa", logical.return=TRUE)){stop("install package 'wmtsa'")}
	# using package "wmtsa"
	#the SearchFrac parameter just controls how much to look to either side 
	#of wavCWTPeaks()'s estimated maxima for a bigger value
	#see dRange
	Wave <- wavCWT(Data)
	WaveTree <- wavCWTTree(Wave)
	WavePeaks <- wavCWTPeaks(WaveTree, snr.min=5)
	WavePeaks_Times <- attr(WavePeaks, which="peaks")[,"iendtime"]
	
	NewPeakTimes <- c()
	dRange <- round(SearchFrac*length(Data))
	for(i in 1:length(WavePeaks_Times)){
		NewRange <- max(c(WavePeaks_Times[i]-dRange, 1)):min(c(WavePeaks_Times[i]+dRange, length(Data)))
		NewPeakTimes[i] <- which.max(Data[NewRange])+NewRange[1]-1
	}
	
	return(matrix(c(NewPeakTimes, Data[NewPeakTimes]), ncol=2, dimnames=list(NULL, c("PeakIndices", "Peaks"))))
}

# ==============================
# = Replace non-finite w/ NA's =
# ==============================
Inf2NA <- function(x) {x[which(x==-Inf | x==Inf, arr.ind=TRUE)] <- NA; x}


# ============================
# = Remove linear time trend =
# ============================
Stationary <- function(x){
	if(is.null(dim(x))){
		x1 <- 0:(length(x)-1)
		y1 <- x
	}else{
		if(dim(x)[2]!=2){stop("dim of x must be 2, or x must be 1 dimensional")}
		x1 <- x[,1] - min(x[,1], na.rm=TRUE)
		y1 <- x[,2]
	}
	
	if(sum(!is.na(y1))<5){return(x)}
	res <- as.vector(residuals(lm(y1~x1, na.action=na.exclude))) + sum(y1, na.rm=TRUE)/length(y1)
	
	if(!is.null(dim(x))){
		x[,2] <- res
		return(x)
	}else{
		return(cbind(x1, res))
	}
}

# ========================================
# = Basic waiting time functions for GEV =
# ========================================
ReturnLevel <- function(RetQ, mu, sc, xi){
	mu-(sc/xi)*(1-(-log(RetQ))^(-xi)) 	
}
ReturnQuant <- function(RetTime, nExts, TS_Duration){
	-TS_Duration/(RetTime*nExts) + 1
}
ReturnTime <- function(RetQ, nExts, TS_Duration){
	1/((1-RetQ)*(nExts/TS_Duration))
}

xYr_Lvl <- function(RetTime, nExts, TS_Duration, mu, sc, xi){
	RetQ <- -TS_Duration/(RetTime*nExts) + 1
	mu-(sc/xi)*(1-(-log(RetQ))^(-xi))
}

# ======================================================
# = Calculate the 1.1*max(x) return level for data.stat =
# ======================================================
calc.level <- function(x, thresh=1.1){
	# calculate the number of observations,
	# the duration of the time series
	# and the level 2 return level
	noNA <- is.finite(x[,"Data"])
	N <- sum(noNA)
	Dur <- diff(range(x[noNA,"year4"]))+1
	luse <- max(x[,"Data"], na.rm=TRUE)*thresh
	return(data.frame("N"=N, "Duration"=Dur, "level"=luse))
		
}

# =============================================================
# = Calculate waiting time for a given return level (for GEV) =
# =============================================================
lvlX_ReturnTime <- function(lvl, a, nExts, TS_Duration){ #Added _v4
	if(any(is.na(c(lvl,a)) | is.nan(c(lvl,a)))){return(NA)}
	
	mu <- a[1]
	sc <- a[2]
	xi <- a[3]
	
	if(xi < 0 & lvl > (mu-sc/xi)){return(NA)}
	
	tx <- (1 + xi * (lvl - mu)/sc)^(-1/xi) #the t(x) from Wikipedia GEV article
	#pdf <- (1/sc)*tx^(xi+1)*exp(-tx) #this is the pdf, if I ever want to use it
	cdf <- exp(-tx) #this is the CDF, which gives the percentile for a certain return level
	1/((1-cdf)*(nExts/TS_Duration))
}

# =================================================
# = wrapper for calculating waiting time from GEV =
# =================================================
lvl_return <- function(x, level, returnFull=FALSE){ # New version to calculate return time (equiv to level 2)
	lvl0 <- x[,"level"] #as.numeric(x[paste("Level",level, "_residual", sep="")])
	a0 <- as.numeric(x[c("mu_0","sig_0","sh_0")])
	nExts0 <- as.numeric(x["N"])
	TS_Duration0 <- as.numeric(x["Duration"])
	result <- lvlX_ReturnTime(lvl=lvl0, a=a0, nExts=nExts0, TS_Duration=TS_Duration0)
	if(returnFull){
		cbind(x, "Level2_time"=result)
	}else{
		result
	}
	
}

# ============================================
# = Convenient changes to reshape() function =
# ============================================
reshape2 <- function(...){
	a <- reshape(...)
	row.names(a) <- NULL
	a <- a[,!names(a)=="id"]
	a
}


# =============================================
# = Estimate mu and sd for normal + lognormal =
# =============================================
musd <- function(x){
	if(!is.null(dim(x))){
		x <- x[,"Data"]
	}
	px <- x>0
	lmu <- mean(log(x[px]), na.rm=TRUE)
	lsd <- sd(log(x[px]), na.rm=TRUE)
	nmu <- mean(x, na.rm=TRUE)
	nsd <- sd(x, na.rm=TRUE)
	data.frame("logMean"=lmu, "logSd"=lsd, "mean"=nmu, "sd"=nsd)
	
}

# =================================================
# = calculate waiting time from level in data.stat =
# =================================================
normTime <- function(x, Level=2){
	RetQ <- pnorm(x[,"level"], x[,"mean"], x[,"sd"])
	data.frame("Level2_normTime"=1/((1-RetQ)*(x[,"N"]/x[,"Duration"])))
}

# =================================================
# = calculate waiting time from level in data.stat =
# =================================================
lnormTime <- function(x, Level=2){
	RetQ <- plnorm(x[,"level"], x[,"logMean"], x[,"logSd"])
	data.frame("Level2_logNormTime"=1/((1-RetQ)*(x[,"N"]/x[,"Duration"])))
}

# ==========================
# = The GEV Wrapper for Me =
# ==========================
GEV <- function(x){
	tryCatch(
		{
			qdat <- x[,"Data"][!is.na(x[,"Data"])]
			qgev <- gev.fit(qdat)
			c(qgev$mle[c("mu_0","sig_0","sh_0")], "se"=qgev$se["sh_0"])
			# rvs <- c(qgev$mle[c("mu_0","sig_0","sh_0")], "se"=qgev$se["sh_0"])
			# data.frame("mu_0"=rvs[1], "sig_0"=rvs[2], "sh_0"=rvs[3], "se.sh_0"=revs[4])
		}, 
		error=function(cond)rep(NA,4)
	)
}


# ========================================================
# = Fit the stable distribution using package stabledist =
# ========================================================
stable.fit <- function(x){
	if(!library("fBasics", logical.return=TRUE)){stop("install package 'fBasics'")}
	sFit <- stableFit(x, type="mle", doplot=FALSE)
	sFit@fit$estimate
}

# ==================================================
# = Calculate waiting time for stable distribution =
# ==================================================
stableTime <- function(lvl, a, nExts, TS_Duration){ #Added _v4
	if(any(is.na(c(lvl,a)) | is.nan(c(lvl,a)))){return(NA)}
	
	alpha <- a[1]
	beta <- a[2]
	cc <- a[3]
	mu <- a[4]
	
	# if(xi < 0 & lvl > (mu-sc/xi)){return(NA)}
	# tx <- (1 + xi * (lvl - mu)/sc)^(-1/xi) #the t(x) from Wikipedia GEV article
	#pdf <- (1/sc)*tx^(xi+1)*exp(-tx) #this is the pdf, if I ever want to use it
	cdf <- pstable(lvl, alpha=alpha, beta=beta, gamma=cc, delta=mu) #this is the CDF, which gives the percentile for a certain return level
	1/((1-cdf)*(nExts/TS_Duration))
}


# ===============================================================
# = The largest Xi that the estimated Xi is signif greater than =
# ===============================================================
fattestSig <- function(mu, se){
	qnorm(p=0.05, mean=mu, sd=se)
}

